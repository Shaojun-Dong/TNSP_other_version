module fix_type
	use Tensor_type
	use usefull_function
	implicit none
contains

	type(Tensor) function linear(X,B)
		type(Tensor),intent(in)::X,B
		linear= (.inv.((.H.X)*X))* ((.H.X)*B)
	end function
	type(Tensor) function  runfix2(A)
		type(Tensor),intent(in)::A
		type(Tensor)::X,B
		integer::i,j
		call X%allocate([A%dim(1),2],A%getType())
		call B%allocate([A%dim(1)],A%getType())
		do i=1,A%dim(1)
			call X%setValue([i,1],A%i([i,1]))
			call X%setValue([i,2],1)
			call B%setValue([i],A%i([i,2]))
		end do
		runfix2=linear(X,B)
		return
	end function
	type(Tensor) function  runfix3(A)
		type(Tensor),intent(in)::A
		type(Tensor)::X,B
		integer::i,j
		call X%allocate([A%dim(1),3],A%getType())
		call B%allocate([A%dim(1)],A%getType())
		do i=1,A%dim(1)
			call X%setValue([i,1],A%i([i,1])*A%i([i,1]))
			call X%setValue([i,2],A%i([i,1]))
			call X%setValue([i,3],1)
			call B%setValue([i],A%i([i,2]))
		end do
		runfix3=linear(X,B)
		return
	end function
	subroutine plotdata3(X)
		type(Tensor),intent(in)::X
		real*8::a,b,c,L,r
		integer::i
		a=X%di(1)
		b=X%di(2)
		c=X%di(3)
		open(unit=123,file='fixData',status='replace')
		L=0
		do i=1,60
			r=a*L*L+b*L+c
			write(123,*)L,r
			L=L+0.005
		end do
		close(123)
		return
	end subroutine
	subroutine plotdata2(X)
		type(Tensor),intent(in)::X
		real*8::a,b,c,d,L,r
		integer::i
		a=X%di(1)
		b=X%di(2)
		open(unit=123,file='fixData',status='replace')
		L=0
		do i=1,60
			r=a*L+b
			write(123,*)L,r
			L=L+0.005
		end do
		close(123)
		
		return
	end subroutine
end module

program aaa
	use Tensor_type
	use usefull_function
	use fix_type
	implicit none
	type(Tensor)::A,X,Zero,A0
	integer::i,Nrow,flag
	real*8::E0,L,Eh
	
	
	call set_error_pointer
	
	
	
	write(*,*)"input row"
	read(*,*)Nrow
	call A%allocate([Nrow,2],'real*8')
	open(unit=234,file='input.dat',status='old')
	call A%readData(234)
	close(234)
	call A%print()
	do i=1,A%dim(1)
		call A%setValue([i,1],1d0/A%di([i,1]))
	end do
	write(*,*)"line:1. x^2:2"
	read(*,*)flag
	if(flag.eq.1)then
		X=runfix2(A)
		call plotdata2(X)
	else
		X=runfix3(A)
		call plotdata3(X)
	end if
	call X%print()
end 












