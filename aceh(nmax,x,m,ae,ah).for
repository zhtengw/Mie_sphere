      subroutine aceh(nlayer,nmax,x,m,an,bn)
	integer nlayer,n,nmax,nk
      double precision x(nlayer)
      complex*16 m(nlayer+1),xi,xo,cc          ! mΪ������
	complex*16 jnx(-1:nmax),ynx(-1:nmax),jnz(-1:nmax),ynz(-1:nmax)
	complex*16 hnx(-1:nmax),hnz(-1:nmax)
	complex*16 tns(1:nmax),dns(1:nmax),tns1,dns1
      complex*16 djnx,djnz,dhnx,dhnz,dynx,dynz
      complex*16 :: an(1:nmax),bn(1:nmax)
	!complex*16 :: ce(1:nmax),ch(1:nmax)
	cc=(0.0,1.0)
	tns=0
	dns=0

	do nk=1, nlayer-1

	   xi=x(nk)*m(nk)
	   xo=x(nk)*m(nk+1)
	   n=nmax !bessel�������ù淶����һ�����������ñ���
	   call c_bessel(n,nmax,jnx,xi)
	   n=nmax
	   call c_newman(n,nmax,ynx,xi)
	   do n=-1,nmax
	      hnx(n)=jnx(n)+cc*ynx(n)
	   end do
	   n=nmax
	   call c_bessel(n,nmax,jnz,xo) 
	   n=nmax
	   call c_newman(n,nmax,ynz,xo)
	   do n=-1,nmax,1
	      hnz(n)=jnz(n)+cc*ynz(n)
	   end do
	   do n=1,nmax,1
	      djnx=xi*jnx(n-1)-n*jnx(n)               ! d��ʾ��(xi*jnxi)��
	      dynx=xi*ynx(n-1)-n*ynx(n)               ! d��ʾ��(xi*ynxi)��
            dhnx=xi*hnx(n-1)-n*hnx(n)               ! d��ʾ��(xi*hnxi)��
	      djnz=xo*jnz(n-1)-n*jnz(n)               ! d��ʾ��(xo*jnxo)��
	      dynz=xo*ynz(n-1)-n*ynz(n)               ! d��ʾ��(xo*ynxo)��
            dhnz=xo*hnz(n-1)-n*hnz(n)               ! d��ʾ��(xo*hnxo)��
	      tns1=tns(n)
	      dns1=dns(n)
	tns(n)=-(m(nk+1)*m(nk+1)*jnz(n)*(djnx+tns1*dynx)
     $-m(nk)*m(nk)*djnz*(jnx(n)+tns1*ynx(n)))
	tns(n)=tns(n)/(m(nk+1)*m(nk+1)*ynz(n)*(djnx+tns1*dynx)
     $-m(nk)*m(nk)*dynz*(jnx(n)+tns1*ynx(n)))
	dns(n)=-(jnz(n)*(djnx+dns1*dynx)-djnz*(jnx(n)+dns1*ynx(n)))
	dns(n)=dns(n)/(ynz(n)*(djnx+dns1*dynx)-dynz*(jnx(n)+dns1*ynx(n)))
         end do

	end do

      an(:)=0
	bn(:)=0
	xi=x(nlayer)*m(nlayer)
	xo=x(nlayer)*m(nlayer+1)
	n=nmax !bessel�������ù淶����һ�����������ñ���
	call c_bessel(n,nmax,jnx,xi)
	n=nmax
	call c_newman(n,nmax,ynx,xi)
	do n=-1,nmax
	   hnx(n)=jnx(n)+cc*ynx(n)
	end do
	n=nmax
	call c_bessel(n,nmax,jnz,xo) 
	n=nmax
	call c_newman(n,nmax,ynz,xo)
	do n=-1,nmax,1
	   hnz(n)=jnz(n)+cc*ynz(n)
	end do
	do n=1,nmax,1
	   djnx=xi*jnx(n-1)-n*jnx(n)               ! d��ʾ��(xi*jnxi)��
	   dynx=xi*ynx(n-1)-n*ynx(n)               ! d��ʾ��(xi*ynxi)��
         dhnx=xi*hnx(n-1)-n*hnx(n)               ! d��ʾ��(xi*hnxi)��
	   djnz=xo*jnz(n-1)-n*jnz(n)               ! d��ʾ��(xo*jnxo)��
	   dynz=xo*ynz(n-1)-n*ynz(n)               ! d��ʾ��(xo*ynxo)��
         dhnz=xo*hnz(n-1)-n*hnz(n)               ! d��ʾ��(xo*hnxo)��
	   tns1=tns(n)
	   dns1=dns(n)
	an(n)=-(m(nlayer+1)*m(nlayer+1)*jnz(n)*(djnx+tns1*dynx)
     $-m(nlayer)*m(nlayer)*djnz*(jnx(n)+tns1*ynx(n)))
	an(n)=an(n)/(m(nlayer+1)*m(nlayer+1)*hnz(n)*(djnx+tns1*dynx)
     $-m(nlayer)*m(nlayer)*dhnz*(jnx(n)+tns1*ynx(n)))
	bn(n)=-(jnz(n)*(djnx+dns1*dynx)-djnz*(jnx(n)+dns1*ynx(n)))
	bn(n)=bn(n)/(hnz(n)*(djnx+dns1*dynx)-dhnz*(jnx(n)+dns1*ynx(n)))
      end do

      return
      end