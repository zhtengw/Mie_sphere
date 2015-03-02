	program main
	
	
	implicit none
	integer :: nlayer=3	!球壳的层数

	integer i,j,n,nargs,state
	double precision,allocatable :: r0(:)		!各层球壳的外径
	real*4 :: executetime(2),clocktime,etime
	character(len=128) :: buffer
	
	nargs=iargc()
c	nargs保存参数个数
c	按设计，第一个参数是整数的球壳层数，
c	之后的参数是从内到外每一层的外径
c
c	如果参数个数大于等于1，即输入了壳层
c	数，那么将第一个参数读到整型nlayer中
	if (nargs .ge. 1) then
	  call getarg(1,buffer)
	  read(buffer,*,iostat=state)nlayer
	endif
c	
c	如果读入错误或者nlayer不是大于1的整数
c	又或者参数个数小于1，从键盘读入nlayer
	j=nargs
	do while(state .ne. 0 .or. j .lt. 1 .or. nlayer .lt. 1)
	  print *, 'Enter number of layers of the sphere(>=1):'
	  read(*,*,iostat=state)nlayer
	  j=j+1
	enddo
	print*,'nlayer',nlayer

	allocate(r0(nlayer))	!给半径数组r0分配内存
	

c
c	跟读入nlayer时相似。数组r0的维数是nlayer，
c	先尝试从参数表的剩余参数中读入双精度的r,
c	第i个半径由第i+1个参数提供
	do i=1,nlayer
	  if ( i+1 .le. nargs) then
	      call getarg(i+1,buffer)	      
	      read(buffer,*,iostat=state)r0(i)
	      print *,buffer
	  endif
	  !如果读入数据错误，或者输入的半径参数不够，
	  !则从键盘读入半径
	  do while(state .ne. 0 .or. i+1 .gt. nargs)
	      print '(a24,i2,a17)', 'Enter the radius of the ',i,
     &' layer (unit:nm):'
	      read(*,*,iostat=state)r0(i)
	      if(state .eq. 0)exit	!读入半径符合实数时跳出循环
	  enddo
	  print '(a3,i2,a3,d12.4)', 'r0(',i,')  ', r0(i)
	  r0(i)=r0(i)*1.0d-9		!输入数据是nm单位，转换为m
	enddo
	
	call Mie_eff(nlayer,r0)
	
	deallocate(r0)
	
	clocktime=etime(executetime)
	write(*,*) executetime(1)
	end
	
