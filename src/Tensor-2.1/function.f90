module usefull_function
	implicit none
	integer,private,save::randomseed
	logical,private,save::seed_flag=.false.
	integer,private,save::max_len_of_char=100
	interface delta
		module procedure delta_real
		module procedure delta_int
		module procedure delta_com
	end interface
	interface sort!Bubble Sort
		module procedure sort1
		module procedure sort2
		module procedure sort11
		module procedure sort22
		module procedure sort3
		module procedure sort4
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
	interface allocateCheck! if size(A)<lenA then allocate A,else do nothing
		module procedure allocateCheck_int
		module procedure allocateCheck_real
		module procedure allocateCheck_com
	end interface
	
	interface operator(+)
		module procedure charAdd
		module procedure charAddint
		module procedure intAddchar
		module procedure charAddreal
		module procedure realAddchar
	end interface
	
	interface assignment(=)
		module procedure charset
	end interface	
!icopy(N,SX,INCX,SY,INCY) : BLAS scopy
contains
	subroutine set_max_len_of_cha(maxlen)
		integer,intent(in)::maxlen
		max_len_of_char=maxlen
		return
	end subroutine
	
	
	subroutine allocateCheck_int(A,lenA)! if size(A)<lenA then allocate A,else do nothing
		integer,allocatable,intent(inout)::A(:)
		integer::lenA
		if(allocated(A)) then
			if(size(A).lt.lenA) then
				deallocate(A)
				allocate(A(lenA))
			end if
		else
			allocate(A(lenA))
		end if
		return
	end subroutine
	subroutine allocateCheck_real(A,lenA)! if size(A)<lenA then allocate A,else do nothing
		real*8,allocatable,intent(inout)::A(:)
		integer::lenA
		if(allocated(A)) then
			if(size(A).lt.lenA) then
				deallocate(A)
				allocate(A(lenA))
			end if
		else
			allocate(A(lenA))
		end if
		return
	end subroutine
	subroutine allocateCheck_com(A,lenA)! if size(A)<lenA then allocate A,else do nothing
		complex*16,allocatable,intent(inout)::A(:)
		integer::lenA
		if(allocated(A)) then
			if(size(A).lt.lenA) then
				deallocate(A)
				allocate(A(lenA))
			end if
		else
			allocate(A(lenA))
		end if
		return
	end subroutine
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
	subroutine sort11(a,inde,realpart,increase)
		complex*16,intent(inout) :: a(:)
		integer,allocatable,intent(out):: inde(:)
		logical,intent(in)::realpart,increase
		complex*16 :: temp
		integer :: i,j,n,tempi
		n=size(a)
		allocate(inde(n))
		do i=1,n
			inde(i)=i
		end do
		if(increase) then
			call sort1(a,inde,realpart)
			return
		end if
		if(realpart) then
			do i=1,n-1
				 do j=i+1,n
					  if (dreal(a(i)) .lt. dreal(a(j))) then
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
					  if (aimag(a(i)) .lt. aimag(a(j))) then
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
	subroutine sort22(a,realpart,increase)
		complex*16,intent(inout) :: a(:)
		logical,intent(in)::realpart,increase
		complex*16 :: temp
		integer :: i,j,n
		n=size(a)
		if(increase) then
			call sort2(a,realpart)
			return
		end if
		if(realpart) then
			do i=1,n-1
				 do j=i+1,n
					  if (dreal(a(i)) .le. dreal(a(j))) then
						   temp = a(i)
						   a(i) = a(j)
						   a(j) = temp
					  endif
				 enddo
			enddo
		else
			do i=1,n-1
				 do j=i+1,n
					  if (aimag(a(i)) .le. aimag(a(j))) then
						   temp = a(i)
						   a(i) = a(j)
						   a(j) = temp
					  endif
				 enddo
			enddo
		end if
		return
	end subroutine
! use for type real*8	
	subroutine sort3(a,inde,increase)
		real*8,intent(inout) :: a(:)
		integer,allocatable,intent(out):: inde(:)
		logical,intent(in)::increase
		real*8 :: temp
		integer :: i,j,n,tempi
		n=size(a)
		allocate(inde(n))
		do i=1,n
			inde(i)=i
		end do
		if(increase) then
			do i=1,n-1
				 do j=i+1,n
					  if (a(i) .gt. a(j)) then
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
					  if (a(i) .lt. a(j)) then
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
	
	subroutine sort4(a,increase)
		real*8,intent(inout) :: a(:)
		logical,intent(in)::increase
		real*8 :: temp
		integer :: i,j,n
		n=size(a)
		if(increase) then
			do i=1,n-1
				 do j=i+1,n
					  if (a(i) .gt. a(j)) then
							temp = a(i)
							a(i) = a(j)
							a(j) = temp
					  endif
				 enddo
			enddo
		else
			do i=1,n-1
				 do j=i+1,n
					  if (a(i) .lt. a(j)) then
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
	
	SUBROUTINE icopy(N,SX,INCX,SY,INCY)
		INTEGER::incx,incy,n
		INTEGER::sx(*),sy(*)
		INTEGER i,ix,iy,m,mp1
		INTRINSIC mod
		IF (n.LE.0) RETURN
		IF (incx.EQ.1 .AND. incy.EQ.1) THEN
      !
      !        code for both increments equal to 1
      !
      !
      !        clean-up loop
      !
			m = mod(n,7)
			IF (m.NE.0) THEN
				DO i = 1,m
					sy(i) = sx(i)
				END DO
				IF (n.LT.7) RETURN
			END IF   
			mp1 = m + 1
			DO i = mp1,n,7
				sy(i) = sx(i)
				sy(i+1) = sx(i+1)
				sy(i+2) = sx(i+2)
				sy(i+3) = sx(i+3)
				sy(i+4) = sx(i+4)
				sy(i+5) = sx(i+5)
				sy(i+6) = sx(i+6)
			END DO
		ELSE      
      !
      !        code for unequal increments or equal increments
      !          not equal to 1
      !
			ix = 1
			iy = 1
			IF (incx.LT.0) ix = (-n+1)*incx + 1
			IF (incy.LT.0) iy = (-n+1)*incy + 1
			DO i = 1,n
				sy(iy) = sx(ix)
				ix = ix + incx
				iy = iy + incy
			END DO
		END IF
		RETURN
  	END SUBROUTINE
		
!	character(len=(max(len(w1)+len(w2),max_len_of_char))) function charAdd(w1,w2)
	character(len=max_len_of_char) function charAdd(w1,w2)
		character(len=*),intent(in)::w1,w2
		charAdd=(trim(w1))//(trim(w2))
		!A='A'+'B'  					~~~~ using time:1d-2
		!A=(trim('A'))//(trim('B')) ~~~~ using time:1d-4
		return
	end function
	subroutine charset(w,inte)
		character(len=*),intent(inout)::w
		integer,intent(in)::inte
		write(w,'(I4)')inte
		return
	end subroutine
!	character(len=(max(len(w1),max_len_of_char))) function charAddint(w1,inte)
	character(len=max_len_of_char) function charAddint(w1,inte)
		character(len=*),intent(in)::w1
		integer,intent(in)::inte
		character*500::temp
		write(temp,'(I4)')inte
		charAddint=(trim(w1))//(trim(adjustl(temp)))
		return
	end function
!	character(len=(max(len(w1),max_len_of_char))) function intAddchar(inte,w1)
	character(len=max_len_of_char) function intAddchar(inte,w1)
		character(len=*),intent(in)::w1
		integer,intent(in)::inte
		character*500::temp
		write(temp,'(I4)')inte
		intAddchar=(trim(adjustl(temp)))//(trim(w1))
		return
	end function
	
	character(len=max_len_of_char) function charAddreal(w1,B)
		character(len=*),intent(in)::w1
		real*8,intent(in)::B
		character*500::temp
		write(temp,'(3F20.16)')B
		charAddreal=(trim(w1))//(trim(adjustl(temp)))
		return
	end function
	
	character(len=max_len_of_char) function realAddchar(B,w1)
		character(len=*),intent(in)::w1
		real*8,intent(in)::B
		character*500::temp
		write(temp,'(3F10.6)')B
		realAddchar=(trim(adjustl(temp)))//(trim(w1))
		return
	end function
	
!  Purpose: Generate seed for each process
!
!====================================================	
	subroutine seed_gen(myseed)
		implicit none
		integer::seed0,counter,myseed
		call system_clock(count=counter)
		seed0=-counter
		myseed= int(-5970*rand(seed0)) 
	end subroutine 
!=========================================================================
!  Purpose:   
!  "Minimal" random number generator of Park and Miller with
!  Bays-Durham shuffle and added safeguards. Returns a uniform
!  random deviate between 0.0 and 1.0 (exclusive of the end points).
!  Call with idum a negative integer to initialize; thereafter do
!  not alter idum between successive deviates in a sequence. RNMX
!  should approximate the largest floating point value that is 
!  less than 1.
!
!  Input: idum, the seed
!  Output:  random number between 0.0 and 1.0
!============================================================================
	real*8 function rand(idum)
		implicit none
		  
		integer::idum     
		
		integer ia, im, iq, ir, ntab, ndiv
		real am, eps, rnmx
		
		parameter (ia=16807, im=2147483647, am=1./im)
		parameter (iq=127773, ir=2836, ntab=32, ndiv=1+(im-1)/ntab)
		parameter (eps=1.2e-7, rnmx=1.-eps)  

		integer j, k, iv(ntab), iy
		save iv, iy

		data iv /ntab*0/, iy /0/

		if (idum.le.0 .or. iy.eq.0) then
		   idum = max(-idum,1)
		   do j = ntab+8, 1, -1
		      k = idum/iq
		      idum = ia*(idum-k*iq)-ir*k
		      if (idum .lt. 0) idum = idum + im
		      if (j .le. ntab) iv(j) = idum
		   enddo
		   iy = iv(1)
		endif

		k = idum/iq
		idum = ia*(idum-k*iq) - ir*k
		if (idum.lt.0) idum = idum + im
		j = 1 + iy/ndiv
		iy = iv(j)
		iv(j) = idum
		rand = min(am*iy,rnmx)
		return
	end function
	
	real*8 function randomnumber()
		if(seed_flag) then
			randomnumber=rand(randomseed)
		else
			call seed_gen(randomseed)
			seed_flag=.true.
			randomnumber=rand(randomseed)
		end if
		return
	end function
	
	subroutine set_seed(idum)
		integer,intent(in)::idum
		if(idum.eq.0)then
			seed_flag=.false.
			return
		end if
		randomseed=idum
		seed_flag=.true.
		return
	end subroutine
	integer function out_randomseed()
		out_randomseed=randomseed
	end function
	integer function out_and_set_seed()
		call seed_gen(randomseed)
		seed_flag=.true.
		out_and_set_seed=randomseed
	end function
end module
