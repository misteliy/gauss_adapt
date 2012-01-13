	program test_gauss
    use determinante
    use random
	implicit none

	integer,parameter                  :: n=10,maxiter=10000 
	character(60)                      :: task
	real(kind=dp)                      :: f,c_t,x(n),eta(n),m(n),eps
	real(kind=dp)                      :: Q(n,n),C(n,n),r,xmin(n),fmin,detQ
	integer                            :: i,counter(2),indx(n)
	logical                            :: iprint


    call init_random_seed()
	task='start'
	eps=1e-9
!	x=1._dp 
    call random_number(x)
    x=(x-0.5_dp)*10._dp
	counter=0
	iprint = .true.
	
	if (iprint) print *,'**********************start************************'
	do i = 0,maxiter
		if (task(1:5).eq.'new_x') call fun(x,n,f)
		if (task.eq.'new_x_accepted') then
			if (iprint) print *,'accepted'
			counter(1)=counter(1)+1
		else
			if (iprint) print *,'rejected'
			counter(2)=counter(2)+1
		endif
		if (iprint) then
			print *,'---------------------------------'
			print *,'loop',i
			write(*,'(A,1X,10E13.6E1)') 'x',x
            print *,'f',f,'ct',c_t
!			print *,'m',m
!			print *,'xmin',xmin,'fmin',fmin
			print *,'r',r
!			print *,'Q',Q
            if (i.ne.0) call DTRM(Q,N,detQ,INDX)
            print *,'detQ',detQ
!			print *,'C',matmul(r*Q,transpose(r*Q))
			print *,'---------------------------------'
            write(111,*) i,xmin,c_t
		endif
		call GaussAdapt(n,x,eta,f,m,C,Q,r,c_t,xmin,fmin,task)
		if (task(1:4).eq.'stop') then
			print *,task 
			exit
		endif
		if (abs(fmin).le.eps) then
			if (iprint) print *,'min found'
			exit
		endif
	enddo
	print *,'=============================================================='
	print *,'final x',xmin,'final f',fmin
	print *,'accepted',counter(1),'rejected',counter(2),'iterations',counter(2)+counter(1)
	end program
!----------------------------objective function---------------------------------
	subroutine fun(x,n,f)
	use random
	implicit none
	integer,intent(in)   :: n
	real(kind=dp) :: x(n),a,f
    f = sum(x**2)
!	f = (x(1)-5._dp)**2 + (x(2)-3._dp)**2
	!a=sqrt(2._dp)
	!f = 100*(x(2)-x(1)**2)**2+(a-x(1))
	end subroutine fun
!-----------------------Gaussian Adaption Algorithm-----------------------------
	subroutine GaussAdapt(n,x,eta,f,m,C,Q,r,c_t,xmin,fmin,task)
	use random
!    use ziggurat
	implicit none
	character(60),intent(inout)        :: task
	integer,intent(in)                 :: n
	real(kind=dp),intent(inout)        :: f,c_t,x(n),eta(n),m(n)
	real(kind=dp),intent(inout)        :: Q(n,n),C(n,n),r
	real(kind=dp)                      :: dC(n,n),Id(n,n)
	real(kind=dp)                      :: wm,wc,wt,beta,p,fe,fc,dx(n)
	real(kind=dp)                      :: fmin,xmin(n),eig_val(n),work(3*n-1)
	real(kind=dp)                      :: detQ,dQ(n,n),test,temp(n,n),delta(n)
	integer                            :: i,j,info
!----------------------------set GauAdapt prameter------------------------------
	Id=0._dp
    do i=1,n; Id(i,i)=1._dp; enddo
    wm=exp(1._dp)*n
!    wm=1._dp
	wc=(n+1._dp)**2/(log(n+1._dp))
	wt=exp(1._dp)*n
!    wt=1._dp
	beta=1._dp/wc
	p=1._dp/exp(1._dp)
    fe=1._dp + beta*(1._dp-p)
!	print *,'f_e',fe
     fc=1._dp-beta*p
!	print *,'f_c',fc
!-------------------------------------------------------------------------------
	if (task(1:5).eq.'new_x') then
		if (f .lt. c_t) then 
            delta = x - m
			task='new_x_accepted'
!----------------------------save current best----------------------------------
			if (f .lt. fmin) then
			fmin = f 
			xmin = x 
			endif
!------------------------------adapt parameters---------------------------------
			r=fe*r
			c_T = (1._dp-1._dp/wt)*c_T+f/wt
			m = (1._dp-1._dp/wm)*m + x/wm
!-----------------------------update convariance--------------------------------

!-----------------------------rank 1 update ------------------------------------
			call DGER(n,n,1._dp/wc,delta,int(1),delta,int(1),(1._dp-1._dp/wc)*C,n)
			dQ=C
!-----------------------------eig.val decomposition-----------------------------
			call DSYEV( 'V', 'l', n, dQ, n, eig_val, WORK, 3*n-1, INFO )
			if (info.ne.0) task='stop_eigendecomp'
            if (any(isnan(eig_val))) task='stop'
            if (any(eig_val.lt.0)) then
                task='stop eig_val'
                return 
            endif
!-----------------------------normalization of Q and C--------------------------
			detQ=product(eig_val)
            detQ=detQ**(1.d0/n)
            C=C/detQ
            eig_val=eig_val/detQ
            eig_val=sqrt(eig_val)
!-----------------------------update Q such that Q=Q*D^1/2------------------------
            do i=1,n; Q(:,i)=dQ(:,i)*sqrt(eig_val(i)); enddo
!--------------------------------update Q---------------------------------------
		else
!---------------rejected lower step size dont adopt cov. + mean-----------------
			task='new_x_rejected'
			r=fc*r
		endif
!------------------------------------sample-------------------------------------
		if (r.lt.1e-9) then
            task='stop_r'
            return
        endif
		do i=1,n; eta(i)=random_normal(); enddo
!		print *,'eta',eta
		x=m
		call DGEMV('N',n,n,r,Q,n,eta,int(1),1._dp,x,int(1))
		if (any(isnan(x))) task='stop_sampling x'
	endif
	
	if (task(1:5).eq.'start') then
		Q = 0._dp; do i = 1,n; Q(i,i) = 1._dp; end do
        C = Q
!		print *,'ini C',C
		r=1._dp
		m = x
		task='new_x'
		call fun(x,n,f)
		c_t=f
		fmin=f
		xmin=x
	endif
	
	end subroutine


SUBROUTINE init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
END SUBROUTINE
