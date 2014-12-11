	subroutine getfilelines(filename,lines)
	implicit none
	integer :: lines
	character*6 :: filename
	character(len=79) :: buffer
	
	lines=0
	open(10, file=filename, status='old')
	do while(.true.)
	lines=lines + 1
	read(10,*,end=100) buffer
	end do
100	close(10)
	return
	end
