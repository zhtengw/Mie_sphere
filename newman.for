      subroutine newman(n,np,bessy,x)
c .......................................................................
c .    subroutine to calculate the spherical bessel function of the	.
c .    second kind (or Newmann function) from order -1			.
c .    to order n; result stored in array bessy.  bess_y  is the common	.
c .    notation for the Newmann function.				.
c .......................................................................

	implicit double precision (a-h,o-z)
      integer n,np
      double precision x,xx,pi
      double precision bessy(-1:np)

      if (n.lt.3) then 
            n=3
      endif
      pi=1.0
      pi=4.*datan(pi)
      xx=x-2.*pi*int(x/2./pi)
c ..................................................................
c .Must use this xx because it's a value between 0 and 2pi, and it .
c .works as the argument for the sin and cos; if we just stuck in  .
c .x instead, it won't work.					   .
c ..................................................................
      bessy(-1)=dsin(xx)/x
      bessy(0)=-dcos(xx)/x
      bessy(1)=-dcos(xx)/x/x-dsin(xx)/x
      do 10 i=2,n
            bessy(i)=(2.*(i-1.)+1.)/x*bessy(i-1)-bessy(i-2)
 10   continue

      return
      end

c	**************************************************

      subroutine c_newman(n,np,bessy,x)
c .......................................................................
c . subroutine to calculate the complex spherical bessel function 	.
c . from order 0 to order n, the result is stored in array bessy.  	.
c . This is the Newman function, aka, bessel function of order 2.	.
c .......................................................................

	implicit double precision(a-h,o-z)

      integer n,np
      complex *16 x,xx,bessy(-1:np)
      double precision pi

      if (n.lt.3) then 
            n=3
      endif
      pi=1.0
      pi=4.*datan(pi)
      xx=x-2.*pi*int(x/2./pi)
      bessy(-1)=cdsin(xx)/x
      bessy(0)=-cdcos(xx)/x
      bessy(1)=-cdcos(xx)/x/x-cdsin(xx)/x
      do 10 i=2,n
            bessy(i)=(2.*(i-1.)+1.)/x*bessy(i-1)-bessy(i-2)
 10   continue
      return
      end

c	**************************************************