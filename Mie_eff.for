	subroutine Mie_eff(nlayer, r0)
c ...................................................................................
c ...................................................................................
	implicit none
	integer, intent(in):: nlayer	!球壳的层数，只能在调用Mie_eff时
					!读入，在执行Mie_eff过程中不可改变
	integer kstep,kn
	double precision Lambda,kstart,kint,Qext,Qscat,Qabs
	double precision r0(nlayer)		!各层球壳的外径
	double precision Cw, PE
	complex*16 :: m(nlayer+1)		!各个区域的折射率
	integer i,j, INFO
	integer, allocatable :: IPIV(:)		!可变大小的数组IPIV?
	real*4 :: executetime(2),clocktime,etime

	Cw=8.0655D5	
	!r0(1)=5.5d-9				!redius of CdSe core, m
	!r0(2)=7.5d-9				!radius of silica shell, m
	!r0(3)=10.0d-9				!radius of nanoparticle, m
      	print*,r0
	kstart=2.00d-7				!初始波长200nm
	kint=1.0d-9				!步长1nm
	kstep=2001				!扫601个波长，即算到800nm

	open(20,file='Mie_sphere.csv')
c	char(9)是ascii码横向制表符
	write(20,*) 'wavelength',char(9),'Qext'
     &  ,char(9),'Qscat',char(9),'Qabs'

	do kn=1,kstep
	Lambda=kstart+(kn-1)*kint		!wavelength, m
c ...................................................................
c refractive index of silica core
c ...................................................................

	m(1)=2.5				!refractive index of CdSe core
	m(2)=1.53

!	m(2)=1.44851+3171.95/((lambda*1.0D9)**2)
!     $+3.516*1.0D7/((lambda*1.0D9)**4)		!refractive index of Silica
!	m(2)=1.32334+3479.0/((lambda*1.0D9)**2)	
!     $-5.111*1.0D7/((lambda*1.0D9)**4) 		!refractive index of water
	PE=1/(Cw*Lambda)
	call getnk(PE,m(3))			!refractive index of Au
c .........................................................................
c When using getnk, copy nk(gold).txt or nk(silver).txt to nk.txt.  
c .........................................................................
	m(4)=1					!refractive index of surrounding medium (air)
c .........................................................................
c the formula of refractive indices of silica and water were both given in 
c B. Khlebtsov et. al, Nanotechnology 17, 5167(2006)
c .........................................................................
	call GetQeffs(Lambda, nlayer, m, r0, Qext, Qscat, Qabs)	
	
c	输出文件的波长以nm为单位，各个量以制表符分割
	write(20,100) Lambda*1.0d9,char(9),Qext,char(9),Qscat,char(9),Qabs
100	format(ES12.5,a,ES12.5,a,ES12.5,a,ES12.5)
	end do

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine GetQeffs(Lambda, nlayer, m, r0, Qext, Qscat, Qabs)
c ...........................................................................................
c  subroutine GetQeffs caculates the factors of extinction, scattering and absorbation of 
c  multi-layer nanoparticle. 
c  This subroutine is based on Mie theory, the multilayer expansion coefficents an and bn
c  are calculated by J. Sinzig and M. Quinten, Appl. Phys. A 58, 157(1994)
c
c  Parameters:
c
c  Integer:
c
c    Lambda -- (input) the wavelength in the vacuum of the laser beam
c
c    nlayer -- (input) the number of layers of the spherical nanoparticle
c
c  complex*16:
c
c    m -- (input) dimension of nlayer+1, the refractive index of every layer and surrounding 
c         medium, m(0) is the core, m(nlayer+1) is the surrounding medium
c
c  double precision
c
c    r0 -- (input) dimension of nlayer, the radius of every layer, r0(0) is the core
c
c    Qext, Qscat, Qabs -- (output) the calculated factors of extinction, scattering and absorbation 
c    of the whole nanoparticle. 
c
c    
c ...................................................................................
	implicit none
      integer nmax               !实际截断的阶数
	integer nk, nlayer              
      double precision x(nlayer),temp1, Lambda, r0(nlayer),pi
	complex*16 m(nlayer+1),cc
      complex*16, allocatable :: an(:),bn(:)
	double precision Qext,Qscat,Qabs
      cc=(0.0D0,1.0D0)
	pi=4*ATAN(1.0D0)

      x=2.0D0*pi/Lambda*r0
	nmax=60
	allocate(an(1:nmax), bn(1:nmax))
	an=0
	bn=0

	call aceh(nlayer,nmax,x,m,an,bn)  !计算单球的散射场展开系数an/bn

	Qext=0
	Qscat=0
	Qabs=0

      do nk=1,nmax,1
         Qext=Qext+float(2*nk+1)*(an(nk)+bn(nk))
	   temp1=float(2*nk+1)*(an(nk)*conjg(an(nk))+bn(nk)*conjg(bn(nk)))
	   Qscat=Qscat+temp1
	end do
	deallocate(an, bn)

	Qext=-1*Qext*2.0/((m(nlayer+1)*x(nlayer))**2)
	Qscat=Qscat*2.0/((m(nlayer+1)*x(nlayer))**2)
	Qabs=Qext-Qscat
      return
	end subroutine

