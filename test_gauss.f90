	program test_gauss
	implicit none

	integer,parameter                  :: dp=8,n=2,maxiter=10000
	character(60)                      :: task
	real(kind=dp)                      :: f,c_t,x(n),m(n),eps
	real(kind=dp)                      :: Q(n,n),r,xmin(n),fmin
	integer                            :: i,counter(2)
	logical                            :: iprint

	task='start'
	eps=1e-6
	x=3._dp
	counter=0
	iprint = .true.
	
	if (iprint) print *,'*********************start************************'
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
			print *,'C',r**2*matmul(Q,transpose(Q))
			print *,'---------------------------------'
		endif
		call GaussAdapt(n,x,f,m,Q,r,c_t,xmin,fmin,task)
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
	f = (x(1)-1._dp)**2 + (x(2)-3._dp)**2
	!a=sqrt(2._dp)
	!f = 100*(x(2)-x(1)**2)**2+(a-x(1))
	end subroutine fun
!-----------------------Gaussian Adaption Algorithm-----------------------------
	subroutine GaussAdapt(n,x,f,m,Q,r,c_t,xmin,fmin,task)
	use random
	implicit none
	character(60),intent(inout)        :: task
	integer,intent(in)                 :: n
	real(kind=dp),intent(inout)        :: f,c_t,x(n),m(n)
	real(kind=dp),intent(inout)        :: Q(n,n),r
	real(kind=dp)                      :: Q_old(n,n),dC(n,n),eta(n),Id(n,n)
	real(kind=dp)                      :: wm,wc,wt,beta,p,fe,fc,dx(n)
	real(kind=dp)                      :: fmin,xmin(n),eig_val(n),work(3*n-1)
	real(kind=dp)                      :: detQ,dQ(n,n)
	integer                            :: i,j,info
!----------------------------set GauAdapt prameter------------------------------
	Id = 0._dp; do i = 1,n; Id(i,i) = 1._dp; end do
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
			Q_old=Q
!-----------------------------update convariance--------------------------------
			call DGER(n,n,1._dp,dx,int(1),dx,int(1),dQ,n)
			dC=(1._dp-1._dp/wc)*Id+1._dp/wc*dQ
			dQ=dC
			call DSYEV( 'V', 'L', n, dQ, n, eig_val, WORK, 3*n-1, INFO )
			if (info.ne.0) task='stop_eigendecomp'
!			print *,'dQ',dQ
!			print *,'eig_val',eig_val
			do i=1,n; Id(i,i)=Id(i,i)*sqrt(eig_val(i)); enddo
			Q=matmul(dQ,matmul(Id,transpose(dQ)))
!--------------------------------normalize dQ-----------------------------------
			detQ=product(eig_val)
			detQ=detQ**(1D+0/n)
			dQ=dQ/detQ
!			print *,detQ,'detQ'
!--------------------------------update Q---------------------------------------
			Q=matmul(Q_old,dQ)
		else
!---------------rejected lower step size dont adopt cov. + mean-----------------
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
