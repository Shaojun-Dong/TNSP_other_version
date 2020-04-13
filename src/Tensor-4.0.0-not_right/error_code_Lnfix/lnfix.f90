module fix_type
	use Tensor_type
	use Tools
	implicit none
contains

	type(Tensor) function linear(X,B)
		type(Tensor),intent(in)::X,B
		linear= (.inv.((.H.X)*X))* ((.H.X)*B)
	end function
	type(Tensor) function  runfixn(A,n)
		type(Tensor),intent(in)::A
		integer,intent(in)::n
		type(Tensor)::X,B
		integer::i,j
		real*8::xn,x0
		call X%allocate([A%dim(1),n],A%getType())
		call B%allocate([A%dim(1)],A%getType())
		do i=1,A%dim(1)
			x0=A%di([i,1])
			xn=1d0
			do j=1,n
				call X%setValue([i,j],xn)
				xn=xn*x0
			end do
			call B%setValue([i],A%i([i,2]))
		end do
		runfixn=linear(X,B)
		return
	end function
	subroutine plotdatan(A,minx,maxx,step)
		type(Tensor),intent(in)::A
		real*8,intent(in)::minx,maxx
		integer,intent(in)::step
		real*8::L,delta,r
		integer::i
		open(unit=123,file='fixData',status='replace')
		delta=(maxx-minx)/(step-1)
		L=minx
		do i=1,step+step/10
			r=point(A,L)
			write(123,*)L,r
			L=L+delta
		end do
		close(123)
		return
	end subroutine
	
	real*8 function point(A,X_)
		type(Tensor),intent(in)::A
		real*8,intent(in)::X_
		type(Tensor)::X
		real*8::x0
		integer::i,n
		n=A%getTotalData()
		call X%allocate([n],'real*8')
		x0=1d0
		do i=1,n
			call X%setValue([i],x0)
			x0=x0*X_
		end do
		point=A*X
		return
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
	subroutine plotdata3(X,minx,maxx,step)
		type(Tensor),intent(in)::X
		real*8,intent(in)::minx,maxx
		integer,intent(in)::step
		real*8::a,b,c,L,r,delta
		integer::i
		a=X%di(1)
		b=X%di(2)
		c=X%di(3)
		open(unit=123,file='fixData',status='replace')
		delta=(maxx-minx)/(step-1)
		L=minx
		do i=1,step
			r=a*L*L+b*L+c
			write(123,*)L,r
			L=L+delta
		end do
		close(123)
		return
	end subroutine
	subroutine plotdata2(X,minx,maxx,step)
		type(Tensor),intent(in)::X
		real*8::a,b,c,d,L,r,delta
		real*8,intent(in)::minx,maxx
		integer,intent(in)::step
		integer::i
		a=X%di(1)
		b=X%di(2)
		open(unit=123,file='fixData',status='replace')
		delta=(maxx-minx)/(step-1)
		L=minx
		do i=1,step
			r=a*L+b
			write(123,*)L,r
			L=L+delta
		end do
		close(123)
		
		return
	end subroutine
end module

program aaa
	use Tensor_type
	use Tools
	use fix_type
	implicit none
	type(Tensor)::A,X,Zero,A0
	integer::i,Nrow,flag,step,n
	real*8::E0,L,Eh,minx,maxx,x0
	call set_error_pointer
	write(*,*)"input row"
	read(*,*)Nrow
	call A%allocate([Nrow,2],'real*8')
	open(unit=234,file='input.dat',status='old')
	call A%readData(234)
	close(234)
	do i=1,A%dim(1)
		call A%setValue([i,1],1d0/A%di([i,1]))
	end do
	call A%print()
	write(*,*)"Fix point x=0"
	x0=0
	write(*,*)"fix n?"
	read(*,*)n
	n=n+1
	if(n.gt.Nrow)then
		call writemess('ERROR n')
		stop
	end if
	A0=A%subTensor(2,1)
	minx=min(x0,A0%dmin())
	maxx=max(x0,A0%dmax())
	step=300
	X=runfixn(A,n)
	call plotdatan(X,minx,maxx,step)
	call X%print()
	E0=point(X,x0)
	write(*,*)E0
	
end 












