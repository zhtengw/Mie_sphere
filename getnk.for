      subroutine getnk(PE,nk)
	implicit none
	integer i
	double precision :: nn(37),kk(37),PEI(37)
	double precision PE,nninterp,kkinterp
	complex*16 nk
      
	open(50,file='nk.txt',status='old')
	do i=1,37
	read(50,*) PEI(i),nn(i),kk(i)
	end do
	close(50)

	CALL PCHIP(37,PEI,nn,1,PE,nninterp)
      CALL PCHIP(37,PEI,kk,1,PE,kkinterp)
	nk=cmplx(nninterp,kkinterp)
	return
	end

c***************************************************************
	SUBROUTINE PCHIP(N, X, F, NE, XE, FE)
	IMPLICIT NONE
	INTEGER N,NE,IERR
	DOUBLE PRECISION X(N),F(N),D(N)
	DOUBLE PRECISION XE(NE),FE(NE)
      logical SKIP
      SKIP=.TRUE.
	CALL DPCHIM (N, X, F, D, 1, IERR)
	CALL DPCHFE (N, X, F, D, 1, SKIP, NE, XE, FE, IERR)
	RETURN
	END


	