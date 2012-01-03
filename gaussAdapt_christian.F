#include "fintrf.h"
*
#if 0
*     
*     gaussAdapt.F
*     .F file needs to be preprocessed to generate .for equivalent
*     
#endif
*     
*     gaussAdapt.F


*#######################################################################
*     gaussAdapt
*#######################################################################

      subroutine gaussAdapt(fitFun, opts, xMin, fMin, countEval, output)
         use userdefined
         use gaAModule
         implicit none


*=======================================================================
*     Interface declaration
*=======================================================================

      interface

*     mvnTruncEigCov
*-----------------------------------------------------------------------
      subroutine mvnTruncEigCov(mu, sigma, eigC, lambdaC, A, b, 
     +   numSamples, burnin, thinInterval, X, Z)
         use userdefined
         implicit none

      type(vector), intent(in) :: mu
      real*8, intent(in) :: sigma
      type(matrix), intent(in) :: eigC
      type(vector), intent(in) :: lambdaC
      type(matrix), intent(in) :: A
      type(vector), intent(in) :: b
      mwSize, intent(in) :: numSamples, burnin, thinInterval
      
      type(matrix), intent(inout) :: X
      type(matrix), intent(inout), optional :: Z

      end subroutine

      end interface

*=======================================================================
*     Variables declaration
*=======================================================================

      mwPointer, intent(in) :: fitFun
      type(gaAOpts), intent(in) :: opts
      type(vector), intent(inout) :: xMin
      real*8, intent(inout) :: fMin
      integer*4, intent(inout) :: countEval
      mwPointer, intent(out) :: output

      integer*4 :: ret, mexEvalString

      type(vector) :: b, mu, lambdaC, work, d, xStart, xOld

      type(matrix) :: A, C, eigC, arx, ary, arz
      real*8 :: r, c_T, f, xsi, ss, sf, N_mu, N_C, N_T

      integer*4 :: nrhs
      integer*4, parameter :: nrhsDisp4 = 4, nrhsDisp3 = 3, nrhsPlot = 2
      mwPointer, allocatable, dimension(:) :: prhs
      mwPointer, dimension(4) :: prhsDisp
      mwPointer, dimension(2) :: prhsPlot
      mwPointer :: ptrX, ptrXplot

      integer*4, parameter :: nlhs = 1, nlhsDisp = 0, nlhsPlot = 0
      mwPointer, dimension(nlhs) :: plhs
      mwPointer, dimension(0) :: NULL

      mwSize, parameter :: numSamples = 1
      mwSize, parameter :: burnin = 0, thinInterval = 1

      real*8 :: condC, detC, scaleC

      type(matrix) :: XDisp, YDisp, FDisp

      type(matrix) :: xRaw, xAcc, xBest, muSave, lambdaSave
      type(vector) :: fRaw, fAcc, fBest, cntAcc
      type(vector) :: rSave, c_TSave, PSave

      integer*4 :: numAcc, rNumAcc, histLen
      real*8 :: rfMin
      type(vector) :: rxMin

      logical :: continueLoop, restartNeeded

      real*8 :: biasVal
      mwSize :: interimSaves
      integer*4 :: totalSaves, accIndex, saveCount
      integer*4, allocatable, dimension(:) :: saveIndex
      real*8, allocatable, dimension(:) :: dSaveIndex

      integer*4, parameter :: figPlot = 43, figDisp = 24

      integer*4, parameter :: structDim = 12
      character(len=30), dimension(structDim) :: fieldnames
      mwPointer, dimension(structDim) :: fields

      integer*4 :: info
      real*8, dimension(8) :: bar 

      character(80) :: foo
      
      integer*4 :: alloc_status, dealloc_status, i, j  
   
      integer*4 :: mexPrintf, mexCallMATLAB
      mwPointer :: mxGetPr, mxCreateDoubleMatrix, mxCreateStructMatrix
      mwPointer :: mxCreateDoubleScalar, mxCreateString
      integer*4 :: mxIsDouble

      mwSize, parameter :: TWO = 2
      mwSize, parameter :: ONE = 1
      integer*4, parameter :: ZERO = 0

      mwIndex :: mwIndexTmp

      intrinsic :: sqrt, dble

#if TIMING_GA
      integer*8 :: now, past, rate
      integer*8, dimension(8) :: time

      time = 0

      call SYSTEM_CLOCK(past, rate)
      call SYSTEM_CLOCK(now)
      write(foo,'(A,D10.3,A)') '0) ',dble(now-past)*1D+3/dble(rate),'ms'
      ret = mexPrintf( foo//achar(10) )
#endif


*=======================================================================
*     Initialize variables
*=======================================================================

      mu%m = xMin%m
      allocate( mu%data(mu%m), stat = alloc_status )
      xStart%m = xMin%m
      allocate( xStart%data(xStart%m), stat = alloc_status )
      xOld%m = xMin%m
      allocate( xOld%data(xOld%m), stat = alloc_status )
      rxMin%m = xMin%m
      allocate( rxMin%data(rxMin%m), stat = alloc_status )

      arx%m = xMin%m
      arx%n = numSamples
      allocate( arx%data(arx%m, arx%n), stat = alloc_status )
      if( .not. opts%bounds%Active ) then
         ary%m = xMin%m
         ary%n = numSamples
         allocate( ary%data(ary%m, ary%n), stat = alloc_status )
         arz%m = xMin%m
         arz%n = numSamples
         allocate( arz%data(arz%m, arz%n), stat = alloc_status )
      endif
      d%m = xMin%m
      allocate( d%data(d%m), stat = alloc_status )

      C%m = xMin%m
      C%n = xMin%m
      allocate( C%data(C%m, C%n), stat = alloc_status )

      lambdaC%m = xMin%m
      allocate( lambdaC%data(lambdaC%m), stat = alloc_status )
      eigC%m = xMin%m
      eigC%n = xMin%m
      allocate( eigC%data(eigC%m, eigC%n), stat = alloc_status )

      if( opts%bounds%Active ) then
         A = opts%bounds%A
         b = opts%bounds%b
      end if

      ptrX = mxCreateDoubleMatrix( xMin%m, ONE, ZERO )
      ptrXplot = mxCreateDoubleMatrix( TWO, ONE, ZERO )

      nrhs = 2 + opts%algorithm%nFunArgs
      allocate( prhs(nrhs), stat = alloc_status )
      prhs(1) = fitFun
      prhs(2) = ptrX
      do i = 3, nrhs
         prhs(i) = opts%algorithm%FunArgs(i-2)
      end do


*     Displaying
*-----------------------------------------------------------------------
      if( opts%output%Display ) then
         XDisp%m = opts%output%length
         XDisp%n = opts%output%length
         allocate( XDisp%data(XDisp%m, XDisp%n), stat = alloc_status )
         YDisp%m = opts%output%length
         YDisp%n = opts%output%length
         allocate( YDisp%data(YDisp%m, YDisp%n), stat = alloc_status )
         FDisp%m = opts%output%length
         FDisp%n = opts%output%length
         allocate( FDisp%data(FDisp%m, FDisp%n), stat = alloc_status )
      end if

*     Saving
*-----------------------------------------------------------------------
      totalSaves = int( opts%stopCrit%MaxIter / 
     +                       opts%output%SavingModulo ) + 2

      interimSaves = max( int( opts%output%VerboseModulo / 
     +      opts%output%SavingModulo ) , 1) + 1

      if( opts%output%Saving ) then
         if( .not. opts%output%Plotting ) then
            interimSaves = 1
         endif
         allocate( saveIndex(interimSaves), stat = alloc_status )
         allocate( dSaveIndex(interimSaves), stat = alloc_status )

         xRaw%m = xMin%m
         xRaw%n = totalSaves
         allocate( xRaw%data(xRaw%m, xRaw%n), stat = alloc_status )
         fRaw%m = totalSaves
         allocate( fRaw%data(fRaw%m), stat = alloc_status )

         xAcc%m = xMin%m
         xAcc%n = totalSaves
         allocate( xAcc%data(xAcc%m, xAcc%n), stat = alloc_status )
         fAcc%m = totalSaves
         allocate( fAcc%data(fAcc%m), stat = alloc_status )
         cntAcc%m = totalSaves
         allocate( cntAcc%data(cntAcc%m), stat = alloc_status )

         xBest%m = xMin%m
         xBest%n = totalSaves
         allocate( xBest%data(xBest%m, xBest%n), stat = alloc_status )
         fBest%m = totalSaves
         allocate( fBest%data(fBest%m), stat = alloc_status )

         muSave%m = xMin%m
         muSave%n = totalSaves
         allocate( muSave%data(muSave%m, muSave%n), stat=alloc_status )
         lambdaSave%m = xMin%m
         lambdaSave%n = totalSaves
         allocate( lambdaSave%data(lambdaSave%m, lambdaSave%n), 
     +              stat = alloc_status )
         rSave%m = totalSaves
         allocate( rSave%data(rSave%m), stat = alloc_status )
         PSave%m = totalSaves
         allocate( PSave%data(PSave%m), stat = alloc_status )
         c_TSave%m = totalSaves
         allocate( c_TSave%data(c_TSave%m), stat = alloc_status )
      end if

#if TIMING_GA
      past = now

      call SYSTEM_CLOCK(now)
      write(foo,'(A,D10.3,A)') '1) ',dble(now-past)*1D+3/dble(rate),'ms'
      ret = mexPrintf( foo//achar(10) )
#endif


*=======================================================================
*     Initial values
*=======================================================================

      mu%data = xMin%data
      xStart%data = xMin%data

      r = opts%algorithm%r
      ss = opts%algorithm%ss
      sf = opts%algorithm%sf
      N_mu = opts%algorithm%N_mu
      N_C = opts%algorithm%N_C
      N_T = opts%algorithm%N_T

      eigC%data = 0D+0
      lambdaC%data = r**2
      C%data = 0D+0
      do i = 1, C%m
         eigC%data(i,i) = 1D+0
         C%data(i,i) = r**2
      end do

      condC = 1D+0
      countEval = 0
      numAcc = 0
      rNumAcc = 0
      histLen = 10*xMin%m

      continueLoop = .TRUE.
      restartNeeded = .FALSE.

      
*     Determine work for eig decomp
*-----------------------------------------------------------------------
      work%m = -1
      call DSYEV( 'V', 'L', eigC%n, eigC%data, eigC%m, lambdaC%data,
     +      bar, work%m, info)
   
      work%m = int( bar(1) )
      allocate( work%data(work%m), stat = alloc_status )

      
*     Determine start fitness value and threshold
*-----------------------------------------------------------------------
      call mxCopyReal8ToPtr( xMin%data, mxGetPr(ptrX), xMin%m )
      if( opts%mode%Design ) then
         c_T = opts%algorithm%c_T
      else
         ret = mexCallMATLAB( nlhs, plhs, nrhs, prhs, 'feval' )
         if( mxIsDouble( plhs(1) ) ) then
            call mxCopyPtrToReal8( mxGetPr(plhs(1)), c_T, ONE )
            countEval = countEval + 1
         else
            call mexErrMsgTxt( 'Return of the fitFun is not a double' )
         end if
         call mxDestroyArray(plhs(1))
      end if
      fMin = c_T
      f = fMin
      rfMin = c_T


*     Saving
*-----------------------------------------------------------------------
      if( opts%output%Saving ) then

         saveIndex = 1
         saveCount = 2
         accIndex = 2

         xRaw%data = 0D+0
         xRaw%data(:,saveIndex(1)) = xMin%data
         fRaw%data = 0D+0
         fRaw%data(saveIndex(1)) = fMin

         xAcc%data = 0D+0
         xAcc%data(:,1) = xMin%data
         fAcc%data = 0D+0
         fAcc%data(1) = fMin

         xBest%data = 0D+0
         xBest%data(:,saveIndex(1)) = xMin%data
         fBest%data = 0D+0
         fBest%data(saveIndex(1)) = fMin

         muSave%data = 0D+0
         muSave%data(:,saveIndex(1)) = mu%data
         rSave%data = 0D+0
         rSave%data(saveIndex(1)) = r
         c_TSave%data = 0D+0
         c_TSave%data(saveIndex(1)) = c_T
         lambdaSave%data = 0D+0
         lambdaSave%data(:,saveIndex(1)) = r**2
      endif

#if TIMING_GA
      past = now

      call SYSTEM_CLOCK(now)
      write(foo,'(A,D10.3,A)') '2) ',dble(now-past)*1D+3/dble(rate),'ms'
      ret = mexPrintf( foo//achar(10) )
#endif
 

*=======================================================================
*     Initialize Matlab stuff for "Display" and "Plotting"
*=======================================================================

      if( opts%output%Plotting ) then

         write(foo,'(A,I2,A)') 'figure(',figPlot,');'
         ret = mexEvalString( foo )

         prhsPlot(1) = mxCreateDoubleMatrix( interimSaves, ONE, ZERO )
         prhsPlot(2) = mxCreateDoubleMatrix( interimSaves, ONE, ZERO )

      end if

      if( opts%output%Display ) then
         write(foo,'(A,I2,A)') 'figure(',figDisp,');hold on;'
         ret = mexEvalString( foo )

         prhs(2) = ptrXplot
         bar(3) = opts%output%UBounds(1) - opts%output%LBounds(1)
         bar(4) = opts%output%UBounds(2) - opts%output%LBounds(2)
         do j = 1, opts%output%length
            do i = 1, opts%output%length
               XDisp%data(i, j) = opts%output%LBounds(1) +
     +            (bar(3) * (i-1D+0) / (opts%output%length-1D+0))
               YDisp%data(i, j) = opts%output%LBounds(2) +
     +            (bar(4) * (j-1D+0) / (opts%output%length-1D+0))

               bar(1) = XDisp%data(i, j)
               bar(2) = YDisp%data(i, j)
 
               call mxCopyReal8ToPtr( 
     +                  bar, mxGetPr(ptrXplot), TWO )
               ret = mexCallMATLAB( nlhs, plhs, nrhs, prhs, 'feval' )
               call mxCopyPtrToReal8( 
     +                  mxGetPr(plhs(1)), bar(5), ONE )
               call mxDestroyArray(plhs(1))
               FDisp%data(i,j) = bar(5)
            enddo
         enddo
         prhs(2) = ptrX

         prhsDisp(1) = mxCreateDoubleMatrix( XDisp%m, XDisp%n, ZERO )
         prhsDisp(2) = mxCreateDoubleMatrix( YDisp%m, YDisp%n, ZERO )
         prhsDisp(3) = mxCreateDoubleMatrix( FDisp%m, FDisp%n, ZERO )
         prhsDisp(4) = mxCreateDoubleScalar( dble(opts%output%length) )

         call mxCopyReal8ToPtr( XDisp%data, mxGetPr(prhsDisp(1)),
     +          XDisp%m*XDisp%n )
         call mxCopyReal8ToPtr( YDisp%data, mxGetPr(prhsDisp(2)),
     +          YDisp%m*YDisp%n )
         call mxCopyReal8ToPtr( FDisp%data, mxGetPr(prhsDisp(3)),
     +          FDisp%m*FDisp%n )

         if(opts%mode%Optimization .or. opts%mode%Sampling) then
            ret = mexCallMATLAB( 
     +              nlhsDisp, NULL, nrhsDisp4, prhsDisp, 'contour' )
         else
            ret = mexCallMATLAB( 
     +              nlhsDisp, NULL, nrhsDisp3, prhsDisp, 'surf' )
            ret = mexEvalString( 'colormap(''autumn'')' )
            ret = mexEvalString( 'shading flat' )
            ret = mexEvalString( 'view(0,90)' )
            ret = mexEvalString( foo )
            write(foo,'(A,F7.4,A,F7.4,A)') 'xlim([',
     +          opts%output%LBounds(1),' ',opts%output%UBounds(1),'])'
            ret = mexEvalString( foo )
            write(foo,'(A,F7.4,A,F7.4,A)') 'ylim([',
     +          opts%output%LBounds(2),' ',opts%output%UBounds(2),'])'
            ret = mexEvalString( foo )
         endif

      end if  

#if TIMING_GA
      past = now

      call SYSTEM_CLOCK(now)
      write(foo,'(A,D10.3,A)') '3) ',dble(now-past)*1D+3/dble(rate),'ms'
      ret = mexPrintf( foo//achar(10)//'GA = ['//achar(10) )
#endif  

*=======================================================================
*     Algorithm: Gaussian Adaptation
*=======================================================================

      do while((countEval < (opts%stopCrit%MaxIter))
     +   .and. (f > opts%stopCrit%StopFitness)
     +   .and. continueLoop)
         countEval = countEval + 1

#if TIMING_GA
      past = now

      call SYSTEM_CLOCK(now)
      time(1) = time(1) + (now-past)
      if(mod(countEval, opts%output%VerboseModulo) == 0) then
         write(foo,'(D10.3)') dble(time(1))*1D+3/dble(rate)
         ret = mexPrintf( foo//'  ' )
      endif
#endif


*     Sample arx ~ N(mu, C)
*-----------------------------------------------------------------------
         if( opts%mode%Sampling ) then
            xOld%data = arx%data(:,1)
         endif
         if( opts%bounds%Active ) then
            call mvnTruncEigCov(mu, r, eigC, lambdaC, A, b, 
     +          numSamples, burnin, thinInterval, arx)
         else
            call RANDOM_NUMBER(arz%data)
            arz%data = (arz%data*2D+0) - 1D+0
            call VDERFINV( arz%m, arz%data, ary%data)

            ary%data(:,1) = r*(sqrt(lambdaC%data)*ary%data(:,1))
            arx%data(:,1) = mu%data

            call DGEMV('N', eigC%m, eigC%n, 1D+0, eigC%data,
     +          eigC%m, ary%data, int(1), 1D+0, arx%data, int(1))
         endif

#if TIMING_GA
      past = now

      call SYSTEM_CLOCK(now)
      time(2) = time(2) + (now-past)
      if(mod(countEval, opts%output%VerboseModulo) == 0) then
         write(foo,'(D10.3)') dble(time(2))*1D+3/dble(rate)
         ret = mexPrintf( foo//'  ' )
      endif
#endif


*     Eval fitFun
*-----------------------------------------------------------------------
         call mxCopyReal8ToPtr( arx%data, mxGetPr(ptrX), arx%m )
         ret = mexCallMATLAB( nlhs, plhs, nrhs, prhs, 'feval' )
         call mxCopyPtrToReal8( mxGetPr(plhs(1)), f, ONE )
         call mxDestroyArray(plhs(1))

#if TIMING_GA
      past = now

      call SYSTEM_CLOCK(now)
      time(3) = time(3) + (now-past)
      if(mod(countEval, opts%output%VerboseModulo) == 0) then
         write(foo,'(D10.3)') dble(time(3))*1D+3/dble(rate)
         ret = mexPrintf( foo//'  ' )
      endif
#endif


*     Check for acceptable solution
*-----------------------------------------------------------------------
         call RANDOM_NUMBER(xsi)
         if( ((f < c_T) .and. opts%mode%Design) .or.
     +       ((f < c_T) .and. opts%mode%Optimization) .or.
     +       ((xsi < (rfMin / f)) .and. opts%mode%Sampling) ) then

*     Save the current best
*-----------------------------------------------------------------------
            if( f < rfMin ) then
               rfMin = f
               rxMin%data = arx%data(:,1)
            end if
         
            if( f < fMin ) then
               fMin = f
               xMin%data = arx%data(:,1)
            end if


*     Adapt parameters
*-----------------------------------------------------------------------
            r = min( ss * r, opts%algorithm%MaxR )
            if(opts%mode%Optimization) then
               c_T = (1-1/N_T)*c_T + (1/N_T)*f
            endif
            mu%data = (1-1/N_mu)*mu%data + (1/N_mu)*arx%data(:,1)

            d%data = arx%data(:,1) - mu%data
            do i = 1, C%n
               C%data(i:C%n, i) = (1-1/N_C)*C%data(i:C%n, i) +
     +            (1/(N_C*(r**2)))*(d%data(i:C%n)*d%data(i))
            end do


*     Eigendecomp of C
*-----------------------------------------------------------------------
         eigC%data = C%data
      
         call DSYEV( 'V', 'L',eigC%n, eigC%data, eigC%m, lambdaC%data,
     +      work%data, work%m, info)
      
         if( info .ne. 0 ) then
            call mexErrMsgTxt( 'Eigenvector decomposition failed' )
         endif


*     Normalizing C
*-----------------------------------------------------------------------
         detC = PRODUCT( lambdaC%data )

         if( (detC .eq. 0) .or. 
     +       (detC .eq. INF) .or. 
     +       (detC .eq. -INF) .or. 
     +       (detC .eq. NAN) ) then
            write(foo,'(A,D10.3)') 'det(C) = ',detC
            call mexErrMsgTxt( foo )
         endif

         scaleC = detC**(1D+0/C%m)

         C%data = C%data / scaleC
         lambdaC%data = lambdaC%data / scaleC

         condC = lambdaC%data(lambdaC%m) / lambdaC%data(1)
         if( condC > opts%algorithm%MaxCond ) then
            do i = 1, C%m
               C%data(i,i) = C%data(i,i) + 1D+0/C%m
            enddo
         endif


*     Adapt moments
*-----------------------------------------------------------------------
            numAcc = numAcc + 1
            rNumAcc = rNumAcc + 1

            if( opts%output%Saving .and.
     +          mod( numAcc, opts%output%SavingModulo ) == 0) then

               accIndex = accIndex + 1

               xAcc%data(:,accIndex) = arx%data(:,1)
               fAcc%data(accIndex) = f
               cntAcc%data(accIndex) = numAcc
            endif


*     Solution no accepted
*-----------------------------------------------------------------------
         else

            r = sf * r

            if( opts%mode%Sampling ) then
               arx%data(:,1) = xOld%data
            endif

         end if

#if TIMING_GA
      past = now

      call SYSTEM_CLOCK(now)
      time(4) = time(4) + (now-past)
      if(mod(countEval, opts%output%VerboseModulo) == 0) then
         write(foo,'(D10.3)') dble(time(4))*1D+3/dble(rate)
         ret = mexPrintf( foo//'  ' )
      endif
#endif


*     Draw error elipse
*-----------------------------------------------------------------------
         if( opts%output%Display 
     +       .and. mod(countEval, 
     +       opts%output%VerboseModulo) == 0 ) then

            write(foo,'(A,I2,A)') 'figure(',figDisp,');hold on;'
            ret = mexEvalString( foo )

            call drawErrorElipse( C%data(1:2, 1:2), mu%data(1:2), 
     +                            r, opts%output )
         endif


*     Plot
*-----------------------------------------------------------------------
         if( opts%output%Plotting .and.
     +       opts%output%Saving .and.
     +       (mod(countEval, opts%output%VerboseModulo) == 0) ) then

            if( opts%stopCrit%StopFitness > -Inf ) then
               biasVal = -opts%stopCrit%StopFitness
            else
               biasVal = 0D+0
            endif

            dSaveIndex = dble( saveIndex * opts%output%SavingModulo )
            call mxCopyReal8ToPtr(dSaveIndex,
     +               mxGetPr(prhsPlot(1)), interimSaves)

            write(foo,'(A,I2,A)') 'figure(',figPlot,');'
            ret = mexEvalString( foo )

            ret = mexEvalString( 'subplot(2,2,1);'  )
            do i = 1, muSave%m
               call mxCopyReal8ToPtr(
     +               muSave%data(i,saveIndex), 
     +               mxGetPr(prhsPlot(2)), interimSaves)

               ret = mexCallMATLAB( 
     +                  nlhsPlot, NULL, nrhsPlot, prhsPlot, 'plot' )
            enddo
            write(foo,'(A,D9.2,A)') 'title([''Current Best: '',''',
     +                 rfMin,''']);'
            ret = mexEvalString( foo )
            ret = mexEvalString( 'grid on;hold on;' )


            ret = mexEvalString( 'subplot(2,2,2);'  )
            call mxCopyReal8ToPtr(
     +               rSave%data(saveIndex), 
     +               mxGetPr(prhsPlot(2)), interimSaves)
            ret = mexCallMATLAB( 
     +                  nlhsPlot, NULL, nrhsPlot, prhsPlot, 'semilogy' )
            write(foo,'(A,D9.2,A)') 'title([''Current r: '',''',
     +               r,''']);'
            ret = mexEvalString( foo )
            ret = mexEvalString( 'grid on;hold on;' )


            ret = mexEvalString( 'subplot(2,2,3);'  )
            call mxCopyReal8ToPtr(
     +               fBest%data(saveIndex) + biasVal, 
     +               mxGetPr(prhsPlot(2)), interimSaves)
            ret = mexCallMATLAB(
     +                  nlhsPlot, NULL, nrhsPlot, prhsPlot, 'semilogy' )
            ret = mexEvalString( 'grid on;hold on;' )

            call mxCopyReal8ToPtr(
     +               c_TSave%data(saveIndex) + biasVal, 
     +               mxGetPr(prhsPlot(2)), interimSaves)
            ret = mexCallMATLAB( 
     +                  nlhsPlot, NULL, nrhsPlot, prhsPlot, 'semilogy' )
            write(foo,'(A,D9.2,A)') 'title([''Current c_T: '',''',
     +               c_T,''']);'
            ret = mexEvalString( foo )

            ret = mexEvalString( 'subplot(2,2,4);'  )
            do i = 1, lambdaSave%m
               call mxCopyReal8ToPtr(
     +               1D+0 / lambdaSave%data(i,saveIndex), 
     +               mxGetPr(prhsPlot(2)), interimSaves)
               ret = mexCallMATLAB( 
     +                  nlhsPlot, NULL, nrhsPlot, prhsPlot, 'semilogy' )
            enddo
            write(foo,'(A,D9.2,A)') 'title([''Condition of C: '',''',
     +               condC,''']);'
            ret = mexEvalString( foo )
            ret = mexEvalString( 'grid on;hold on;' )

            saveIndex(1) = saveIndex(saveCount)
            saveCount = 1


         endif

#if TIMING_GA
      past = now

      call SYSTEM_CLOCK(now)
      time(5) = time(5) + (now-past)
      if(mod(countEval, opts%output%VerboseModulo) == 0) then
         write(foo,'(D10.3)') dble(time(5))*1D+3/dble(rate)
         ret = mexPrintf( foo//'  ' )
      endif
#endif


*     Restart
*-----------------------------------------------------------------------
         if( opts%mode%Optimization ) then
            if( r < opts%stopCrit%TolR ) then
               write(foo,'(A,D9.2)') 'Restart due to: TolR: ',r
               restartNeeded = .TRUE.
            elseif( (abs(c_T-rfMin) < opts%stopCrit%TolCon) .and. 
     +              (f > opts%stopCrit%StopFitness) .and. 
     +              (rNumAcc > histLen)) then
               write(foo,'(A,D9.2)') 'Restart due to: TolCon: ',
     +               abs(c_T-rfMin)
               restartNeeded = .TRUE.
            endif
         endif

         if( restartNeeded ) then
            continueLoop = opts%stopCrit%Restart 
     +          .and. (f >= opts%stopCrit%StopFitness)
            restartNeeded = .FALSE.

            ret = mexPrintf( foo//achar(10) )

            r = opts%algorithm%r

            eigC%data = 0D+0
            lambdaC%data = r**2
            C%data = 0D+0
            do i = 1, C%m
               eigC%data(i,i) = 1D+0
               C%data(i,i) = r**2
            enddo
            condC = 1D+0

            if( opts%bounds%Active ) then
               call mvnTruncEigCov(xStart, r, eigC, lambdaC, A, b, 
     +                  1, burnin, thinInterval, arx)
               mu%data = arx%data(:,1)
            else
               mu%data = xStart%data
            endif

            if( opts%mode%Design ) then
               c_T = opts%algorithm%c_T
            else
               call mxCopyReal8ToPtr( mu%data, mxGetPr(ptrX), mu%m )
               ret = mexCallMATLAB( nlhs, plhs, nrhs, prhs, 'feval' )
               call mxCopyPtrToReal8( 
     +                  mxGetPr(plhs(1)), c_T, ONE )
               call mxDestroyArray(plhs(1))
               countEval = countEval + 1
            end if

            rNumAcc = 0
            N_T = opts%algorithm%inc_T * N_T
            rfMin = c_T
            rxMin%data = mu%data
            
         endif

#if TIMING_GA
      past = now

      call SYSTEM_CLOCK(now)
      time(6) = time(6) + (now-past)
      if(mod(countEval, opts%output%VerboseModulo) == 0) then
         write(foo,'(D10.3)') dble(time(6))*1D+3/dble(rate)
         ret = mexPrintf( foo//'  ' )
      endif
#endif


*     Save data
*-----------------------------------------------------------------------
         if( opts%output%Saving .and.
     +       mod(countEval, opts%output%SavingModulo) == 0) then

            if( opts%output%Plotting ) then
                saveCount = saveCount + 1
                saveIndex(saveCount) = saveIndex(saveCount-1) + 1
            else
                saveCount = 1
                saveIndex(saveCount) = saveIndex(saveCount) + 1
            endif
            rSave%data(saveIndex(saveCount)) = r
            muSave%data(:,saveIndex(saveCount)) = mu%data
            lambdaSave%data(:,saveIndex(saveCount)) = lambdaC%data
            c_TSave%data(saveIndex(saveCount)) = c_T
            PSave%data(saveIndex(saveCount)) = 
     +              dble(numAcc) / dble(countEval)
            
            fRaw%data(saveIndex(saveCount)) = f
            xRaw%data(:,saveIndex(saveCount)) = arx%data(:,1)
            
            fBest%data(saveIndex(saveCount)) = rfMin
            xBest%data(:,saveIndex(saveCount)) = rxMin%data

         endif

#if TIMING_GA
      past = now

      call SYSTEM_CLOCK(now)
      time(7) = time(7) + (now-past)
      if(mod(countEval, opts%output%VerboseModulo) == 0) then
         write(foo,'(D10.3)') dble(time(7))*1D+3/dble(rate)
         ret = mexPrintf( foo//';'//achar(10) )
         time = 0
      endif
#endif
 

*     Verbose output
*-----------------------------------------------------------------------
         if( opts%output%Printing .and.
     +       (mod(countEval, opts%output%VerboseModulo) == 0)) then
            ret = mexPrintf( '----------------------------'//achar(10) )

            write(foo,'(A,I6)') 'Number of iterations: ',countEval
            ret = mexPrintf( foo//achar(10) )

            write(foo,'(A,D9.2)') 'P_acc: ',dble(numAcc)/dble(countEval)
            ret = mexPrintf( foo//achar(10) )

            write(foo,'(A,D9.2)') 'Search radius: ',r
            ret = mexPrintf( foo//achar(10) )

            if( opts%mode%Optimization ) then
               write(foo,'(A,D9.2)') 'Current fitness: ',f
               ret = mexPrintf( foo//achar(10) )

               write(foo,'(A,D9.2)') 'Threshold: ',c_T
               ret = mexPrintf( foo//achar(10) )

               write(foo,'(A,D9.2)') 'TolCon: ',c_T-rfMin
               ret = mexPrintf( foo//achar(10) )
            endif

         endif

      enddo

#if TIMING_GA
         ret = mexPrintf(']'//achar(10) )
#endif
      

*     Create output structure
*-----------------------------------------------------------------------
 

      if( opts%output%Saving ) then

         fieldnames(1) = 'xRaw'
         fieldnames(2) = 'fRaw'
         fieldnames(3) = 'xAcc'
         fieldnames(4) = 'fAcc'
         fieldnames(5) = 'cntAcc'
         fieldnames(6) = 'xBest'
         fieldnames(7) = 'fBest'
         fieldnames(8) = 'muVec'
         fieldnames(9) = 'eigVec'
         fieldnames(10) = 'rVec'
         fieldnames(11) = 'P_empVec'
         fieldnames(12) = 'c_TVec'
         
         output = mxCreateStructMatrix(ONE, ONE, structDim, fieldnames)

         fields(1) = mxCreateDoubleMatrix( xRaw%m, xRaw%n, ZERO )
         call mxCopyReal8ToPtr( 
     +          xRaw%data, mxGetPr(fields(1)), xRaw%m*xRaw%n )

         fields(2) = mxCreateDoubleMatrix( fRaw%m, ONE, ZERO )
         call mxCopyReal8ToPtr( fRaw%data, mxGetPr(fields(2)), fRaw%m )

         fields(3) = mxCreateDoubleMatrix( xAcc%m, xAcc%n, ZERO )
         call mxCopyReal8ToPtr( 
     +          xAcc%data, mxGetPr(fields(3)), xAcc%m*xAcc%n)

         fields(4) = mxCreateDoubleMatrix( fAcc%m, ONE, ZERO )
         call mxCopyReal8ToPtr( fAcc%data, mxGetPr(fields(4)), fAcc%m )

         fields(5) = mxCreateDoubleMatrix( cntAcc%m, ONE, ZERO )
         call mxCopyReal8ToPtr( 
     +          cntAcc%data, mxGetPr(fields(5)), cntAcc%m)

         fields(6) = mxCreateDoubleMatrix( xBest%m, xBest%n, ZERO )
         call mxCopyReal8ToPtr( 
     +          xBest%data, mxGetPr(fields(6)), xBest%m*xBest%n)

         fields(7) = mxCreateDoubleMatrix( fBest%m, ONE, ZERO )
         call mxCopyReal8ToPtr(fBest%data, mxGetPr(fields(7)), fBest%m)

         fields(8) = mxCreateDoubleMatrix( muSave%m, muSave%n, ZERO )
         call mxCopyReal8ToPtr( 
     +          muSave%data, mxGetPr(fields(8)), muSave%m*muSave%n)

         fields(9) = mxCreateDoubleMatrix(
     +                  lambdaSave%m, lambdaSave%n, ZERO )
         call mxCopyReal8ToPtr( lambdaSave%data, 
     +          mxGetPr(fields(9)), lambdaSave%m*lambdaSave%n)

         fields(10) = mxCreateDoubleMatrix( rSave%m, ONE, ZERO )
         call mxCopyReal8ToPtr( 
     +          rSave%data, mxGetPr(fields(10)), rSave%m)

         fields(11) = mxCreateDoubleMatrix( PSave%m, ONE, ZERO )
         call mxCopyReal8ToPtr( 
     +          PSave%data, mxGetPr(fields(11)), PSave%m)

         fields(12) = mxCreateDoubleMatrix( c_TSave%m, ONE, ZERO )
         call mxCopyReal8ToPtr( 
     +          c_TSave%data, mxGetPr(fields(12)), c_TSave%m)
      
      mwIndexTmp = 1
      do i = 1, structDim 
         call mxSetFieldByNumber(output, mwIndexTmp, i, fields(i))
      enddo

      else
         output = mxCreateString('No output')

      endif



*=======================================================================
*     Finalize
*=======================================================================

      call mxDestroyArray(ptrX)
      call mxDestroyArray(ptrXplot)

      deallocate( mu%data, stat = dealloc_status )
      deallocate( xStart%data, stat = dealloc_status )
      deallocate( xOld%data, stat = dealloc_status )
      deallocate( rxMin%data, stat = dealloc_status )
      deallocate( arx%data, stat = dealloc_status )
      if( .not. opts%bounds%Active ) then
         deallocate( ary%data, stat = dealloc_status )
         deallocate( arz%data, stat = dealloc_status )
      endif
      deallocate( d%data, stat = dealloc_status )
      deallocate( C%data, stat = dealloc_status )
      deallocate( lambdaC%data, stat = dealloc_status )
      deallocate( eigC%data, stat = dealloc_status )
      deallocate( prhs, stat = dealloc_status )

      deallocate( work%data, stat = dealloc_status )


      if( opts%output%Plotting ) then
         call mxDestroyArray(prhsPlot(1))
         call mxDestroyArray(prhsPlot(2)) 
      endif

      if( opts%output%Display ) then
         call mxDestroyArray(prhsDisp(1))
         call mxDestroyArray(prhsDisp(2))
         call mxDestroyArray(prhsDisp(3))
         call mxDestroyArray(prhsDisp(4))

         deallocate( XDisp%data, stat = dealloc_status )
         deallocate( YDisp%data, stat = dealloc_status )
         deallocate( FDisp%data, stat = dealloc_status )
      endif

      if( opts%output%Saving ) then
         deallocate( saveIndex, stat = dealloc_status )
         deallocate( dSaveIndex, stat = dealloc_status )
         deallocate( xRaw%data, stat = dealloc_status )
         deallocate( fRaw%data, stat = dealloc_status )
         deallocate( xAcc%data, stat = dealloc_status )
         deallocate( fAcc%data, stat = dealloc_status )
         deallocate( cntAcc%data, stat = dealloc_status )
         deallocate( xBest%data, stat = dealloc_status )
         deallocate( fBest%data, stat = dealloc_status )
         deallocate( muSave%data, stat = dealloc_status )
         deallocate( lambdaSave%data, stat = dealloc_status )
         deallocate( rSave%data, stat = dealloc_status )
         deallocate( PSave%data, stat = dealloc_status )
         deallocate( c_TSave%data, stat = dealloc_status )
      endif

      return
      end