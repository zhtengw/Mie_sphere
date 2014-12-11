	subroutine bessel(n,np,bessj,x)
c .......................................................................
c .subroutine to calculate spherical bessel functions of the first      .
c .kind from order -1 to order n, the result is stored in the array 	  .
c .bessj.  n is related to the size parameter (nbg), and np is the	  .
c .physical dimension of the bessj array (isize).  x is a real argument .
c .......................................................................

      implicit double precision (a-h,o-z)
      parameter(iacc=500,bigno=1d100,bigni=1d-100)
      integer n,np
      double precision bessj(-1:np),x,xx,pi
      double precision bjp,bj,factor,ax

      pi=1.0
      pi=4.*datan(pi)
      xx=x-2*pi*int(x/2./pi)
      if (n.lt.3) then 
            n=3
      endif

c .......................................................................
c . If x is less than n, we must use the downward recurrence relation, 	.
c . otherwise our results for bessj will be unstable.  In fact, we will .
c . use the downward alot more than the upward because our x will 	.
c . often be less than n, i.e. the size parameter is less than the order.
c .......................................................................
      if (x.lt.n) then
      ax=dabs(x)
      axx=x-2*pi*int(ax/2./pi)
         if (ax.eq.0) then
            bessj(0)=1.
            do 10 i=1, n
                  bessj(i)=0.
 10         continue
         else
            m=2*((n+int(sqrt(float(iacc*n))))/2)
c ........................................................
c . m is the highest order to calculate; the higher m is .
c . the more accuracy we'll get, but we'll keep n terms	 .
c ........................................................
            do 15 j=0,n
                  bessj(j)=0.
 15         continue
            bjp=0.
            bj=1.
            do 20 j=m,1,-1
                  bjm=(2.*j+1.)/ax*bj-bjp
                  bjp=bj
                  bj=bjm
                  if(dabs(bj).gt.bigno) then
c ..........................................................
c . If bj is huge, we must renormalize to prevent overflow .
c ..........................................................
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
                  bessj(j)=bessj(j)/dsqrt(sum)
                  if (x.lt.0.and.mod(n,2).eq.1) bessj(j)=-bessj(j) 
 30        continue
            bessj(-1)=dcos(xx)/x
         endif
      
      else
c .......................................................................
c . Use upward recurrence for when x is greater than n			.
c .......................................................................
            bessj(-1)=dcos(xx)/x
            bessj(0)=dsin(xx)/x
            do 40 i=0,n-1
                 bessj(i+1)=bessj(i)*(2.*i+1.)/x-bessj(i-1)
 40         continue
      endif
      return
      end                                         

c	****************************************************

   
   
   
   
   
   
   


c	****************************************************