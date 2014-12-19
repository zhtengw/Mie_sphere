	subroutine getfilelines(filename,lines)
	implicit none
	integer :: lines
	character*6 :: filename
	character(len=79) :: buffer
	
	lines=0
	open(10, file=filename, status='old')

	do while(.true.)
	read(10,*,end=100) buffer
	lines=lines + 1
c	
c	It's different for gfortran compiler, 
c	if use gfortran newer than 4.9,
c	it should be:
c	read(10,*,end=100) buffer
c	lines=lines + 1

	end do
100	close(10)
	return
	end
