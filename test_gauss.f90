	program test_gauss
	implicit none

	integer,parameter                  :: dp=8,n=2,maxiter=1000000
	character(60)                      :: task
	real(kind=dp)                      :: f,c_t,x(n),eta(n),m(n),eps
	real(kind=dp)                      :: Q(n,n),r,xmin(n),fmin
	integer                            :: i,counter(2)
	logical                            :: iprint

	task='start'
	eps=1e-5
	x=1._dp
	counter=0
	iprint = .true.
	
	if (iprint) print *,'**********************start************************'
	do i = 1,maxiter
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
			print *,'x',x,'f',f,'ct',c_t
			print *,'m',m
			print *,'xmin',xmin,'fmin',fmin
			print *,'r',r
			print *,'Q',Q
            print *,'detQ',q(1,1)*q(2,2)-q(1,2)*q(2,1)
			print *,'C',matmul(r*Q,transpose(r*Q))
			print *,'---------------------------------'
		endif
		call GaussAdapt(n,x,eta,f,m,Q,r,c_t,xmin,fmin,task)
		if (task(1:4).eq.'stop') then
			print *,task 
			exit
		endif
		if (abs(fmin).lt.eps) then
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
	!f = (x(1)-4._dp)**2 + (x(2)-3._dp)**2
	a=sqrt(2._dp)
	f = 100*(x(2)-x(1)**2)**2+(a-x(1))
	end subroutine fun
!-----------------------Gaussian Adaption Algorithm-----------------------------
	subroutine GaussAdapt(n,x,eta,f,m,Q,r,c_t,xmin,fmin,task)
	use random
	implicit none
	character(60),intent(inout)        :: task
	integer,intent(in)                 :: n
	real(kind=dp),intent(inout)        :: f,c_t,x(n),eta(n),m(n)
	real(kind=dp),intent(inout)        :: Q(n,n),r
	real(kind=dp)                      :: dC(n,n),Id(n,n)
	real(kind=dp)                      :: wm,wc,wt,beta,p,fe,fc,dx(n)
	real(kind=dp)                      :: fmin,xmin(n),eig_val(n),work(3*n-1)
	real(kind=dp)                      :: detQ,dQ(n,n),test,temp(n,n)
	integer                            :: i,j,info
!----------------------------set GauAdapt prameter------------------------------
	Id=0._dp
    do i=1,n; Id(i,i)=1._dp; enddo
    wm=exp(1._dp)*n
	wc=(n+1._dp)**2/(log(n+1._dp))
	wt=exp(1._dp)*n
	beta=1._dp/wc
	p=1._dp/exp(1._dp)
     fe=1._dp + beta*(1._dp-p)
!	print *,'f_e',fe
     fc=1._dp-beta*p
!	print *,'f_c',fc
!-------------------------------------------------------------------------------
	if (task(1:5).eq.'new_x') then
		if (f .lt. c_t) then 
			task='new_x_accepted'
!----------------------------save current best----------------------------------
			if (f .lt. fmin) then
			fmin = f 
			xmin = x 
			endif
!------------------------------adapt parameters---------------------------------
			r=fe*r
			c_T = (1-1/wt)*c_T+1/wt*f
			m = (1-1/wm) * m + 1/wm*x
!-----------------------------update convariance--------------------------------
			call DGER(n,n,1._dp,eta,int(1),eta,int(1),dQ,n)
			dC=(1._dp-1._dp/wc)*Id+1._dp/wc*dQ
			dQ=dC
			call DSYEV( 'V', 'l', n, dQ, n, eig_val, WORK, 3*n-1, INFO )
			if (info.ne.0) task='stop_eigendecomp'
            if (any(isnan(eig_val))) task='stop'
            if (any(eig_val.eq.0)) task='stop' 
!            print *,'eig_val',eig_val 
			detQ=product(eig_val)
            detQ=detQ**(1.d0/n)
            eig_val=eig_val/detQ
            eig_val=sqrt(eig_val)
!            print *,'before eigenvalue multi',dQ
            do i=1,n; dQ(:,i)=dQ(:,i)*sqrt(eig_val(i)); enddo
            !test=dq(1,1)*dq(2,2)-dq(1,2)*dq(2,1)
            !print *,test,'test first Qdet'
!            print *,'dQ half',dQ
            !dq= matmul(dQ,transpose(dq))
            !call DGEMM('n','t',n,n,n,1.0d0,dQ,n,dQ,n,0.0d0,temp,n)
            call DSYRK('l','n',n,n,1._dp,dQ,n,0._dp,temp,n)
            dQ=temp
            !do i=1,n
            !   do j=1,i
            !      dq(i,j)=dq(j,i)
            !   enddo
            !enddo
!            print *,'dQ',dQ
            !call DGEMM('n','t',n,n,n,1._dp,Id,n,dQ,n,0._dp,Id,n) 
            !print *,'Id',Id
            !print *,'dQ',dQ
            !dQ=matmul(dQ,Id)
            !call DGEMM('n','n',n,n,n,1._dp,dQ,n,Id,n,0._dp,dQ,n) 
		!	dQ=matmul(dQ,matmul(Id,transpose(dQ)))
!--------------------------------normalize dQ-----------------------------------
			!detQ=product(eig_val)
			!detQ=detQ**(1D+0/n)
			!dQ=dQ/detQ
			!print *,detQ,'detQ'
!            test=dq(1,1)*dq(2,2)-dq(1,2)**2
!            print *,test,'testQdet'
!            print *,'old Q',Q
!--------------------------------update Q---------------------------------------
            call DSYMM('r','l',n,n,1._dp,dQ,n,Q,n,0._dp,temp,n)
            !call DGEMM('n','n',n,n,n,1.0d0,Q,n,dQ,n,0.0d0,temp,n)
			!Q=matmul(Q,dQ)
            Q=temp
!            print *,'final Q',Q
!            test=q(1,1)*q(2,2)-q(1,2)*q(2,1)
!            print *,test,'test final Q det'
		else
!---------------rejected lower step size dont adopt cov. + mean-----------------
            !test=q(1,1)*q(2,2)-q(1,2)*q(2,1)
            !print *,test,'test final Q det'
			task='new_x_rejected'
			r=fc*r
		endif
!------------------------------------sample-------------------------------------
		do i=1,n; eta(i)=random_normal(); enddo
!		print *,'eta',eta
		x=m
		call DGEMV('N',n,n,r,Q,n,eta,int(1),1._dp,x,int(1))
		if (any(isnan(x))) task='stop_sampling'
!		if (r.lt.1e-15) task='stop_r'
	endif
	
	if (task(1:5).eq.'start') then
		Q = 0._dp; do i = 1,n; Q(i,i) = 1._dp; end do
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
