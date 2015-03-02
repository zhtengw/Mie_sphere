	subroutine getnk(PE,nk)
	implicit none
	integer :: i,lines
	double precision, allocatable :: nn(:),kk(:),PEI(:)
	double precision PE,nninterp,kkinterp
	complex*16 nk
	
	call getfilelines('nk.txt', lines)
	allocate(nn(lines),kk(lines),PEI(lines))
        
	open(50,file='nk.txt',status='old')
	do i=1,lines
	read(50,*) PEI(i),nn(i),kk(i)
	end do
	close(50)
c	DPCHIM(N, X, F, D, 1, IERR) 要求数组X必须是从小到大排列的，所以数据nk.txt必须是从低频到高频排序
	CALL PCHIP(lines,PEI,nn,1,PE,nninterp)
	CALL PCHIP(lines,PEI,kk,1,PE,kkinterp)
	nk=cmplx(nninterp,kkinterp)
	
	deallocate(nn,kk,PEI)

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


	
