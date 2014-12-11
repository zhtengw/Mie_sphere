      subroutine c_bessel(n,np,bessj,x)
c .......................................................................
c .subroutine to calculate complex spherical bessel functions of the 	.
c .first kind from order -1 to order n.  The result is stored in array	.
c .bessj.  This is just like the bessel function, except we now have 	.
c .complex arguments, i.e. x is a complex variable	       		.
c .......................................................................
      implicit double precision (a-h,o-z)
      parameter(iacc=500,bigno=1d100,bigni=1d-100)
	double precision pi,xx
      integer n,np
      complex*16 x,xxxx,sum 
      complex*16 bessj(-1:np)
      complex*16 bjp,bj,factor,bjm
      
      if (n.lt.3) then 
            n=3
      endif

      xx=cdabs(x)
      pi=1.0
      pi=4.*datan(pi)
      xxxx=x
      if (xx.gt.32760) then
             xxxx=x-2*pi*int(x/2./pi)
      endif

      if (xx.lt.n) then
c .............................
c . Use downward recurrence   .
c .............................
         if (xx.eq.0) then
            bessj(0)=1.
            do 10 i=1,n
                  bessj(i)=0.
 10         continue
         else
            m=2*((n+int(sqrt(float(iacc*n))))/2)
            do 15 j=0,n
                  bessj(j)=0.
 15         continue
            bjp=dcmplx(0.d0,0.)
            bj=dcmplx(1.d0,1.)
            do 20 j=m,1,-1
                  bjm=(2.*j+1.)/x*bj-bjp
                  bjp=bj
                  bj=bjm
                  xxx=cdabs(bj)
                  if(xxx.gt.bigno) then
c .......................................
c . Renormalize to prevent overflow	.
c .......................................
                        bj=bj*bigni
                        bjp=bjp*bigni
                        do 17 i=0,n
                              bessj(i)=bessj(i)*bigni
 17                     continue
                  endif
                  if (j.le.n) bessj(j)=bjp
 20         continue
            bessj(0)=bj
	    sum=0.
	    do 25 j=0,n
		sum=sum+(2.*j+1.)*bessj(j)*bessj(j)
 25	    continue
            do 30 j=0,n
                  bessj(j)=bessj(j)/cdsqrt(sum)
 30         continue
            bessj(-1)=cdcos(xxxx)/x
        endif

      else
c ...........................
c . Use upward recurrence   .
c ...........................
            bessj(-1)=cdcos(xxxx)/x
            bessj(0)=cdsin(xxxx)/x
            do 40 i=0,n-1
                 bessj(i+1)=bessj(i)*(2.*i+1)/x-bessj(i-1)
 40         continue
      endif
      return
      end                                         