module usefull_function
	implicit none
	interface delta
		module procedure delta_real
		module procedure delta_int
		module procedure delta_com
	end interface
	interface sort!Bubble Sort
		module procedure sort1
		module procedure sort2
	end interface
!**********************************************************
!		A is a allocatable array and B is array too,when use
!	A=B before allocating storage space for A,A will get
! 	no data.This happen fortran90 but do not in fortran77.
	interface assignments
		module procedure assignment_int_dim1
		module procedure assignment_real8_dim2
		module procedure assignment_real4_dim2
		module procedure assignment_com_dim2
		module procedure assignment_com_dim1				
	end interface
	interface maxvalue
		module procedure maxvalue_real
	end interface
contains
!**********************************************************
!		A is a allocatable array and B is array too,when use
!	A=B before allocating storage space for A,A will get
! 	no data.This happen fortran90 but do not in fortran77.
!
	subroutine assignment_int_dim1(a,b)
		integer,allocatable,intent(out)::a(:)
		integer,intent(in)::b(:)
		integer::m
		m=size(b)
		allocate(a(m))
		a=b
		return
	end subroutine
	subroutine assignment_real8_dim2(a,b)
		real*8,allocatable,intent(out)::a(:,:)
		real*8,intent(in)::b(:,:)
		integer::m,n
		m=size(b,1)
		n=size(b,2)
		allocate(a(m,n))
		a=b
		return
	end subroutine
	subroutine assignment_real4_dim2(a,b)
		real*4,allocatable,intent(out)::a(:,:)
		real*4,intent(in)::b(:,:)
		integer::m,n
		m=size(b,1)
		n=size(b,2)
		allocate(a(m,n))
		a=b
		return
	end subroutine
	
	subroutine assignment_com_dim2(a,b)
		complex*16,allocatable,intent(out)::a(:,:)
		complex*16,intent(in)::b(:,:)
		integer::m,n
		m=size(b,1)
		n=size(b,2)
		allocate(a(m,n))
		a=b
		return
	end subroutine
	subroutine assignment_com_dim1(a,b)
		complex*16,allocatable,intent(out)::a(:)
		complex*16,intent(in)::b(:)
		integer::m
		m=size(b)
		allocate(a(m))
		a=b
		return
	end subroutine
!*****************delta *******************************
! delta function
	real*8 function delta_real(a,b)
		real*8,intent(in)::a,b
		if(a.eq.b) then
			delta_real=1d0
		else
			delta_real=0d0
		end if
		return
	end function
	complex*16 function delta_com(a,b)
		complex*16,intent(in)::a,b
		if(a.eq.b) then
			delta_com=dcmplx(1,0)
		else
			delta_com=dcmplx(0,0)
		end if
		return
	end function
	integer function delta_int(a,b)
		integer,intent(in)::a,b
		if(a.eq.b) then
			delta_int=1
		else
			delta_int=0
		end if
		return
	end function
!*****************  Bubble Sort  ****************
!inde output the order of the output
!sort the data form small to big
!if realpart=.true.,then sort base on the real part of the data
! else the imag part
	subroutine sort1(a,inde,realpart)
		complex*16,intent(inout) :: a(:)
		integer,allocatable,intent(out):: inde(:)
		logical,intent(in)::realpart
		complex*16 :: temp
		integer :: i,j,n,tempi
		n=size(a)
		allocate(inde(n))
		do i=1,n
			inde(i)=i
		end do
		if(realpart) then
			do i=1,n-1
				 do j=i+1,n
					  if (dreal(a(i)) .gt. dreal(a(j))) then
						   temp = a(i)
						   tempi=inde(i)
						   a(i) = a(j)
						   inde(i)=inde(j)
						   a(j) = temp
						   inde(j)=tempi
					  endif
				 enddo
			enddo
		else
			do i=1,n-1
				 do j=i+1,n
					  if (aimag(a(i)) .gt. aimag(a(j))) then
						   temp = a(i)
						   tempi=inde(i)
						   a(i) = a(j)
						   inde(i)=inde(j)
						   a(j) = temp
						   inde(j)=tempi
					  endif
				 enddo
			enddo
		end if
		return
	end subroutine
	
	subroutine sort2(a,realpart)
		complex*16,intent(inout) :: a(:)
		logical,intent(in)::realpart
		complex*16 :: temp
		integer :: i,j,n
		n=size(a)
		if(realpart) then
			do i=1,n-1
				 do j=i+1,n
					  if (dreal(a(i)) .gt. dreal(a(j))) then
						   temp = a(i)
						   a(i) = a(j)
						   a(j) = temp
					  endif
				 enddo
			enddo
		else
			do i=1,n-1
				 do j=i+1,n
					  if (aimag(a(i)) .gt. aimag(a(j))) then
						   temp = a(i)
						   a(i) = a(j)
						   a(j) = temp
					  endif
				 enddo
			enddo
		end if
		return
	end subroutine
! find the max value of a
! the line element of a is the max
	subroutine maxvalue_real(a,line,maxel)
		real*8,intent(in)::a(:)
		real*8,intent(inout)::maxel
		integer,intent(inout)::line
		real*8 :: temp
		integer :: i,n
		n=size(a)
		temp=a(1)
		line=1
		do i=2,n
			if(a(i).gt.temp) then
				temp=a(i)
				line=i
			end if
		end do
		maxel=temp
		return
	end subroutine

end module
