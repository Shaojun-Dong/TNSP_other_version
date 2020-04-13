!****************************************************
!****************************************************
!*****************     note      ********************
!****************************************************
!****************************************************
!
!		1. This is the code of Gradient Method. One should
!	write the code of the function,say func(), to be
!	maximum or minimum when calling Gradient_Method.
!	The input of the function should be a one dimension
!	Type(Tensor), and output a real*8 value.
!
!	
!		2. The data in GMparameter are:
!----------------------------------------------------------------------------------------------
!diff_delta			1d-8		differentiation,(f(x+delta)-f(x))/delta,diff_delta is the delta
!max_step			5d0		go along the dirction with max_step,new_x=x+max_step
!delta_step			0.9		if f(new_x)<f(x),max_step=max_step*delta_step,delta_step should smaller than 1 
!stop_error			1d-8		if abs(f(new_x)-f(x))<stop_error,output
!max_running		5000		max running
!print_flag			1			if 1,print the progress
!running_type      2       !=1,along the Gradient,if can not find the min(max), step=step*delta_step
!                          !=2,along the Gradient, go for num=num_run_per_step, then step=step*delta_step
!num_run_per_step  100
!----------------------------------------------------------------------------------------------
!
!		4. Send email to tyrants@qq.com to report any bugs.
!
!****************************************************
!****************************************************

module GradientMethod
	use Tensor_complex
	implicit none
	real*8,private::zero_number=1d-10
	real*8,private::max_step_in_linear_search=1.
	real*8,private::InfinityNumber=1d20
!**********************************************
!*****   GradientMethod parameter   ***********
	private
	public CGMethod,Example_searchFunction,Example_GradientFunction
	type CGMethod
		logical,private::print_flag=.false.
		real*8,private::diff_delta=1d-5
		real*8,private::max_step=0.5
		real*8,private::delta_step=0.9
		real*8,private::stop_error=1d-8
		integer,private::max_running=2000
		integer,private::running_type=1!=1,along the Gradient,if can not find the min(max), step=step*delta_step
!                                    !=2,along the Gradient, go for num=num_run_per_step, then step=step*delta_step
		integer,private::num_run_per_step=20
		character*2,private::diretion='GM'! Do not finished this part
		                              !GM:Gradient
		                              !CG:ConjugateGradient
	contains
		generic,public::run =>Gradient_Method1,Gradient_Method2,Gradient_Method_para,Gradient_Method2_para
		generic,public::initial =>intGradientMethod1,intGradientMethod2,intGradientMethod3
		generic,public::step =>step1,step2,step3,step1_para,step2_para,step3_para
		procedure,public::set_step,set_max_running
		procedure::Gradient_Method1,Gradient_Method2,Gradient_Method_para,Gradient_Method2_para
		procedure::intGradientMethod1,intGradientMethod2,intGradientMethod3
		procedure::step1,step2,step3,step1_para,step2_para,step3_para
	end type CGMethod

!**********************************************


!**********************************************
!		Example of searchFunction
!**********************************************
!
!	real*8 function searchFunction(point)
!		Type(Tensor),intent(in)::point
!
!.......  write the function here .......
!		
!		return
!	end function
!---------------------------------------------
!or
!	real*8 function searchFunction(point,para)
!		Type(Tensor),intent(in)::point,para
!
!.......  write the function here .......
!		
!		return
!	end function
!**********************************************
!**********************************************


!example program
!	real*8 maxpoint
!	Type(Tensor)::point
!	maxpoint=Gradient_Method(point,.true.,3,Example_searchFunction) 
!or
!	!	maxpoint=Gradient_Method(point,para,.true.,3,Example_searchFunction) 
!---------------------------------------------------------

!**********************************************
contains
	type(Tensor) function Example_searchFunction(Ten,para)!this is the function to be minimize or maximize
		type(Tensor),intent(in)::Ten
		type(Tensor),intent(in)::para
		real*8::x,y,z
		x=Ten.i.1
		y=Ten.i.2
		z=Ten.i.3
		Example_searchFunction=para%i(1)*(x-1)**2-(y-2)**2+1+para%i(2)
		return
	end function
	type(Tensor) function Example_GradientFunction(Ten,para)!this is the function to be minimize or maximize
		type(Tensor),intent(in)::Ten
		type(Tensor),intent(in)::para
		real*8::x,y,z
		x=Ten.i.1
		y=Ten.i.2
		z=Ten.i.3
		call Example_GradientFunction%allocate((/3/),'real*8')
		call Example_GradientFunction%setValue(1,para%i(1)*2*(x-1))
		call Example_GradientFunction%setValue(2,-2*(y-2))
		call Example_GradientFunction%setValue(3,0)
		return
	end function
	
	subroutine intGradientMethod1(CG)
		class(CGMethod),intent(inout)::CG
		CHARACTER*100::notused
		integer::p_flag
		open(unit=10,file="GMparameter",STATUS='old')
		read(10,*)notused,CG%diff_delta!differentiation,(f(x+delta)-f(x))/delta,diff_delta is the delta
		read(10,*)notused,CG%max_step!go along the dirction with max_step,new_x=x+max_step
		read(10,*)notused,CG%delta_step!if f(new_x)<f(x),max_step=max_step*delta_step
		read(10,*)notused,CG%stop_error!if abs(f(new_x)-f(x))<stop_error,output
		read(10,*)notused,CG%max_running!max running
		read(10,*)notused,p_flag!if 1,print the progress
		read(10,*)notused,CG%running_type
		read(10,*)notused,CG%num_run_per_step
		if(p_flag.eq.1)then
			CG%print_flag=.true.
		else
			CG%print_flag=.false.
		end if
		close(unit=10)
		return
	end subroutine
	subroutine intGradientMethod2(CG,GMparameter_address)
		class(CGMethod),intent(inout)::CG
		CHARACTER(len=*),intent(in)::GMparameter_address
		CHARACTER*100::notused
		integer::p_flag
		open(unit=10,file=GMparameter_address,STATUS='old')
		read(10,*)notused,CG%diff_delta!differentiation,(f(x+delta)-f(x))/delta,diff_delta is the delta
		read(10,*)notused,CG%max_step!go along the dirction with max_step,new_x=x+max_step
		read(10,*)notused,CG%delta_step!if f(new_x)<f(x),max_step=max_step*delta_step
		read(10,*)notused,CG%stop_error!if abs(f(new_x)-f(x))<stop_error,output
		read(10,*)notused,CG%max_running!max running
		read(10,*)notused,p_flag!if 1,print the progress
		if(p_flag.eq.1)then
			CG%print_flag=.true.
		end if
		close(unit=10)
		return
	end subroutine
	
	subroutine intGradientMethod3(CG,diff_delta_,max_step_,delta_step_,stop_error_,max_running_,p_flag_)
		class(CGMethod),intent(inout)::CG
		logical,intent(in)::p_flag_
		real*8,intent(in)::diff_delta_,max_step_,delta_step_,stop_error_
		integer,intent(in)::max_running_
		CG%diff_delta=diff_delta_!differentiation,(f(x+delta)-f(x))/delta,diff_delta is the delta
		CG%max_step=max_step_!go along the dirction with max_step,new_x=x+max_step
		CG%delta_step=delta_step_!if f(new_x)<f(x),max_step=max_step*delta_step
		CG%stop_error=stop_error_!if abs(f(new_x)-f(x))<stop_error,output
		CG%max_running=max_running_!max running
		CG%print_flag=p_flag_!if 1,print the progress
		return
	end subroutine
	subroutine set_step(CG,max_step_)
		class(CGMethod),intent(inout)::CG
		class(*),intent(in)::max_step_
		select type(max_step_)
			type is (integer)
				CG%max_step=max_step_
			type is (real(kind=4))
				CG%max_step=max_step_
			type is (real(kind=8))
				CG%max_step=max_step_
			type is (complex(kind=4))
				CG%max_step=max_step_
			type is (complex(kind=8))
				CG%max_step=max_step_
		end select
		return
	end subroutine
	
	subroutine set_max_running(CG,max_running_)
		class(CGMethod),intent(inout)::CG
		integer,intent(in)::max_running_
		CG%max_running=max_running_
		return
	end subroutine
		
		
	
	

!**********************************************************
!**********************************************************
!**********************************************************
!**********No not modify the function below****************
!**********************************************************
!**********************************************************
!**********************************************************
	
	
	
	
!*************************************************************
!  CG method: f(P) is the function to be min
!             g(P) is the 	gradient of f(P)
!	1. Find the direction, d1=-g(P1)
!  2. search x, that min f(P1+x*d1)
!    1) Linear search, input random x0
!       and output P2, and x1
!  3. Find the CG direction d_i=-g(P_i)+\lambda*d_{i-1},where 
!                 g(P_i)^T*g(P_i)
!   \lambda=  -----------------------
!              g(P_{i-1})^T*g(P_{i-1})
!   or 
!               g(P_i)^T * ( g(P_i)-g(P_{i-1}) )
!   \lambda=  --------------------------------
!                 g(P_{i-1})^T*g(P_{i-1})
!
!  4. search x, that min f(Pi+x*di)
!    1) Linear search, input x1(not random)
!       and output P_{i+1}, and x_i
!  5. go to 3. untill g(P_i)=0 or max_step
!  6. output P. and f
!=============================================================
!  Linear search: f(x)=f(P+x*d) is a function of x
!                 g(x)= vec{g(P+x*d)} \cdot \vec{d}
!                 g(P+x*d) is a vector but  g(x) not
!
! 1. input x0, it is a random number at the first time
! 2. input         f(P)    --->  f(0)
!                  g(P)    --->  g(0)
! 3. calculate   f(P+x0*d) --->  f(x0)
! 4. calculate g(P+x0*d)*d --->  g(x0)
! 5. suppose f(x)=a*x^3+b*x^2+c*x+d ,use x=0,f(0),g(0),x0,f(x0),g(x0)
!   solve a,b,c,d
! 6. When we have a,b,c,d min f(x)=a*x^3+b*x^2+c*x+d--->output x1
!   1). f'(x)=3*a*x^2+2*b*x+c=0 ==> x1
!   2). if Delta= (2*b)^2-4*(3*a)*c <0, no root of f'(x)=0,
!      use quadratic interpolation
!   3).Quadratic interpolation: 
!      suppose g(x)=kx+b, use x=0,g(0),x0,g(x0) solve k and b
!      x1=-b/k
! 7. calculate   f(P+x1*d) --->  f(x1)
! 8. calculate g(P+x1*d)*d --->  g(x1)
! 9. Now there are 3 piont
!     xa =  0     xb =  x0     xc =  x1
!     fa = f(0)   fb = f(x0)   fc = f(x1)
!     ga = g(0)   gb = g(x0)   gb = g(x1)
!   1).if ga*gb>0, use Quadratic interpolation, output x2
!    (1).suppose g(x)=kx+b, use xa,ga,xb,gb solve k and b
!      x2=-b/k
!    (2). calculate   f(P+x2*d) --->  f(x2)
!    (3). calculate g(P+x2*d)*d --->  g(x2)
!    (4). if fa*fc<0 keep fa ,drop fb
!            xb=xc,fb=fc,gb=gc
!            xc=x2,fc=f2,gc=g2
!        else keep fb ,drop fa
!           xa=xc,fa=fc,ga=gc
!            xc=x2,fc=f2,gc=g2
!          if all point on the sameside,drop xb   
!    (5).  go to 9.
!     
!   2). if ga*gb<0, use cubic interpolation
!     (1). suppose g(x)=a*x^2+b*x+c and f(x)=\int_0^x g(t)dt=a/3*x^3+b/2*x^2+c*x
!     (2). xa,ga,xb,gb,xc,gc to solve a,b,c in g(x)=a*x^2+b*x+c
!     (3). g(x2)=0 find two x2, choose the one that min f(x) as output x2
!     (4). calculate   f(P+x2*d) --->  f(x2)
!     (5). calculate g(P+x2*d)*d --->  g(x2)
!     (6). if fa*fc<0 keep fa ,drop fb
!            xb=xc,fb=fc,gb=gc
!            xc=x2,fc=f2,gc=g2
!        else keep fb ,drop fa
!           xa=xc,fa=fc,ga=gc
!            xc=x2,fc=f2,gc=g2
!       if all point on the sameside,drop xb   
!     (7). go to 9.
! 10. repete 9 about LinearSearch_max_running time, output x
! 11. new P is P+x*d, x will be input x0 in next Linear search circle
!*************************************************************

!*************************************************************
!*************************************************************
!                 linear search
!*************************************************************
!*************************************************************
	!f(x)=ax^3+bx^2+cx+d
	!g(x)=3ax^2+2bx+c
	!input f(0),g(0),f(x1),g(x1),x1==>output a,b,c,d
	subroutine solvefx(a,b,c,d,f0,g0,f1,g1,x1)!use in LinearSearch_third_point
		real*8,intent(inout)::a,b,c,d
		real*8,intent(in)::f0,g0,f1,g1,x1
		real*8::x1_3,temp
		x1_3=x1*x1*x1
		c=g0
		d=f0
		a= ((g1+c)*x1-2d0*(f1-d))/x1_3
		b= (g1-c)/(2d0*x1)-1.5d0*x1*a
		return
	end subroutine
	
	!f(x)=ax^3+bx^2+cx+d
	!input a, b, c, d ,find minf(x) or max fx
	!outx should be larger than 0
	subroutine zerofx(outx,a,b,c,d,Delta,maxflag)!use in LinearSearch_third_point
		real*8,intent(in)::a,b,c,d
		real*8,intent(inout)::outx,Delta
		logical,intent(in)::maxflag
		real*8::x1,x2,f1,f2
		!Delta=4d0*b*b-12*a*c
		if(Delta.lt.zero_number) then
			write(*,*)"ERROR, should call zerolinearfx"
			stop
		end if
		if(dabs(a).le.zero_number)then!f(x)=bx^2+cx+d
			if(dabs(b).le.zero_number) then !f(x)=cx+d
				write(*,*)"ERROR in zerofx"
				stop
			end if
			outx=-c/(2d0*b)
			return
		end if
		Delta=dsqrt(Delta)
		x1=(-2d0*b+Delta)/(6d0*a)
		x2=(-2d0*b-Delta)/(6d0*a)
		f1=a*x1*x1*x1+b*x1*x1+c*x1+d
		f2=a*x2*x2*x2+b*x2*x2+c*x2+d
		if(f1.le.f2)then
			if(maxflag) then
				outx=x2
				if(outx.le.0) outx=x1
			else
				outx=x1
				if(outx.le.0) outx=x2
			end if
		else
			if(maxflag) then
				outx=x1
				if(outx.le.0) outx=x2
			else
				outx=x2
				if(outx.le.0) outx=x1
			end if
		end if
		return
	end subroutine
	!input two point
	!x0=0,f0,g0,x1,f1,g1
	!find the next point
	!suppose f(x)=ax^3+bx^2+cx+d,then g(x)=3ax^2+2bx+c
	!input x0=0,f0,g0,x1,f1,g1 can determin a,b,c,d
	!Then find xmin,g(xmin)=0
	!If there is no root in the equation g(x)=0
	!Then regard (0,g0) and (x1,g1) as two point on a line g(x)=kx+b
	!Then find xmin,g(xmin)=0
	subroutine LinearSearch_third_point(outx,f0,g0,f1,g1,x1,maxflag)
		real*8,intent(inout)::outx
		real*8,intent(in)::f0,g0,f1,g1,x1
		logical,intent(in)::maxflag
		real*8::Delta,a,b,c,d
		call solvefx(a,b,c,d,f0,g0,f1,g1,x1)
		Delta=4d0*b*b-12*a*c
		if(Delta.lt.zero_number) then
			call zerolinearfx(outx,0d0,g0,x1,g1)
			return
		end if
		call zerofx(outx,a,b,c,d,Delta,maxflag)
		return
	end subroutine
		

	!g(x)=kx+b
	!input x1,g1,x2,g2
	!g(x0)=0,output x0
	subroutine zerolinearfx(outx,x1,g1,x2,g2,stopflag)
		real*8,intent(inout)::outx
		real*8,intent(in)::x1,g1,x2,g2
		logical,optional,intent(inout)::stopflag
		real*8::k,b
		k=(g1-g2)/(x1-x2)
		if(dabs(k).le.zero_number)then
			if(present(stopflag)) then
				stopflag=.true.
				return
			else
				write(*,*)"ERROR in zerolinearfx,k=0"
				write(*,*)x1,g1
				write(*,*)x2,g2
				stop
			end if
		end if
		b=g1- (x1*k)
		outx=-b/k
		return
	end subroutine
	
	!on output,fc,gc are useless
	!on input, outx and xc can be the same variable,but other cannot
!
!if (ga*gb) >0
! 
! use b and c ==> d .because c may be obtain from a,b, if use a and b to get d,  c and d may be the same,go wrong
! a,b,c keep two point
!   if(ga*gc<0) keep a and c
!   if(gb*gc<0) keep b and c
!   ga * gb is larger than 0,so other case is ga ,gb ,gc are the same sign
!     in this case ,drop b, keep a ,c
!  
	
!If ga,gb,gc are the same sign, use Quadratic interpolation,otherwhile use cubic interpolation
	subroutine LinearSearch_fouth_point(outx,stopflag,xa,fa,ga,xb,fb,gb,xc_,fc,gc,maxflag)
		real*8,intent(inout)::outx,xa,fa,ga,xb,fb,gb
		logical,intent(inout)::stopflag
		real*8,intent(in)::xc_,fc,gc
		logical,intent(in)::maxflag
		real*8::xc
		real*8::a,b,c,Delta
		logical::Quadraticflag
		!Quadraticflag,Quadratic interpolation
		xc=xc_
		if(ga*gb.ge.0d0) then
			call zerolinearfx(outx,xa,ga,xc,gc,stopflag)
			if(stopflag.or.(outx.le.0d0))then
			 	call choose_fx(outx,xa,xb,xc,fa,fb,fc,maxflag)
			 	stopflag=.true.
			 	return
			end if
			fb=fc
			gb=gc
			xb=xc
			return
		end if
		call solvefx2(a,b,c,xa,ga,xb,gb,xc,gc)
		if(dabs(a).le.zero_number)then!g(x)=bx+c
			if(dabs(b).le.zero_number) then
				call choose_fx(outx,xa,xb,xc,fa,fb,fc,maxflag)
			else
				outx=-c/b
			end if
			stopflag=.true.
			if(outx.le.0d0) call choose_fx(outx,xa,xb,xc,fa,fb,fc,maxflag)
			return
		end if
		Delta=b*b-4.*a*c!g(x)=ax^2+bx+c
		if( (Delta.ge.0) .and.(Delta.lt.InfinityNumber)) then
			call zerolinearfx2(outx,a,b,c,Delta,maxflag)
			if(outx.le.0d0) then
				stopflag=.true.
				call choose_fx(outx,xa,xb,xc,fa,fb,fc,maxflag)
			end if
			if(ga*gc.gt.0d0) then
				ga=gc
				fa=fc
				xa=xc
			end if
			if(gb*gc.gt.0d0) then
				gb=gc
				fb=fc
				xb=xc
			end if
		else! if the g is too small or x are too close,Delta=NAM
			stopflag=.true.
			call choose_fx(outx,xa,xb,xc,fa,fb,fc,maxflag)
		end if
		return
	end subroutine
!outx should be larger than 0	
	subroutine choose_fx(outx,xa,xb,xc,fa,fb,fc,maxflag)
		real*8,intent(in)::xa,xb,xc,fa,fb,fc
		logical,intent(in)::maxflag
		real*8,intent(out)::outx
		real*8::f
		if(maxflag) then
			outx=xa
			f=fa
			if(f.le.fb) then	
				if(xb.ge.0) then
					 outx=xb
					 f=xb
				end if
			end if
			if(f.le.fc) then	
				if(xc.ge.0) then
				 outx=xc
				 f=xc
				end if
			end if
		else
			outx=xa
			f=fa
			if(f.ge.fb) then	
				if(xb.ge.0) then
				 outx=xb
				 f=xb
				end if
			end if
			if(f.ge.fc) then	
				if(xc.ge.0) then
				 outx=xc
				 f=xc
				end if
			end if
		end if
		if(outx.le.0d0)then
			outx=max_step_in_linear_search
		end if
		return
	end subroutine
		
	!g(x)=ax^2+bx+c
	!output a,b,c
	subroutine solvefx2(a,b,c,x1,g1,x2,g2,x3,g3)
		real*8,intent(inout)::a,b,c
		real*8,intent(in)::g1,x1,g2,x2,g3,x3
		real*8::temp1,x1_2,x2_2,x3_2
		x1_2=x1*x1
		x2_2=x2*x2
		x3_2=x3*x3
		temp1= (x2_2-x1_2)/(x3_2-x1_2)
		b=(temp1*(g3-g1)-g2+g1)/(temp1*(x3-x1)-x2+x1)
		temp1= (x2-x1)/(x3-x1)
		a=(temp1*(g3-g1)-g2+g1)/(temp1*(x3_2-x1_2)-x2_2+x1_2)
		c=g3-x3_2*a-x3*b
		return
	end subroutine
	
	!g(x)=ax^2+bx+c ==>f(x)=int_0^x g(t)dt
	!output minf(x) and outx should be larger than 0
	subroutine zerolinearfx2(outx,a,b,c,Delta,maxflag)
		real*8,intent(in)::a,b,c
		real*8,intent(out)::outx,Delta
		logical,intent(in)::maxflag
		real*8::x1,x2,f1,f2
		Delta=dsqrt(Delta)
		if(dabs(a).le.zero_number)then!g(x)=bx+c
			write(*,*)"ERROR in zerolinearfx2"
			stop
		end if
		x1=(-b+Delta)/(2d0*a)
		x2=(-b-Delta)/(2d0*a)
		f1=(a*x1*x1*x1/3d0)+(b*x1*x1*0.5d0)+c*x1
		f2=(a*x2*x2*x2/3d0)+(b*x2*x2*0.5d0)+c*x2
		if(f1.le.f2)then
			if(maxflag) then
				outx=x2
				if(outx.le.0) outx=x1
			else
				outx=x1
				if(outx.le.0) outx=x2
			end if
		else
			if(maxflag) then
				outx=x1
				if(outx.le.0) outx=x2
			else
				outx=x2
				if(outx.le.0) outx=x1
			end if
		end if
		return
	end subroutine
	
	
!differentiation of searchFunction
	type(Tensor) function diff(CG,inputP,ith,GMsearchFunction)
		class(CGMethod),intent(in)::CG
		type(Tensor),intent(in)::inputP
		integer,intent(in)::ith
		type(Tensor),external::GMsearchFunction
		type(Tensor)::delta_P,P2,P1
		integer::lenP
		lenP=inputP%getTotalData()
		call delta_P%allocate((/lenP/),inputP%getType())
		call delta_P%zero()
		call delta_P%setValue(ith,CG%diff_delta)
		P2=inputP+delta_P
		P1=inputP-delta_P
		diff=GMsearchFunction(P2)-GMsearchFunction(P1)
		diff=diff/(CG%diff_delta+CG%diff_delta)
		if(isnan(diff)) then
			write(*,*)"NAN error,diff"
			stop
		end if
		return
	end function
!the search function is GMsearchFunction(inoutP,point_parameter)
!Search inoutP
	type(Tensor) function diff_para(CG,inputP,point_parameter,ith,GMsearchFunction) result(diff)
		class(CGMethod),intent(in)::CG
		type(Tensor),intent(in)::inputP
		type(Tensor),intent(in)::point_parameter
		integer,intent(in)::ith
		type(Tensor),external::GMsearchFunction
		type(Tensor)::delta_P,P2,P1
		integer::lenP
		lenP=inputP%getTotalData()
		call delta_P%allocate((/lenP/),inputP%getType())
		call delta_P%zero()
		call delta_P%setValue(ith,CG%diff_delta)
		P2=inputP+delta_P
		P1=inputP-delta_P
		diff=GMsearchFunction(P2,point_parameter)-GMsearchFunction(P1,point_parameter)
		diff=diff/(CG%diff_delta+CG%diff_delta)
		if(isnan(diff)) then
			write(*,*)"NAN error,diff"
			stop
		end if
		return
	end function
! gradient vector of searchFunction
	type(Tensor) function gradient(CG,P,GMsearchFunction)
		class(CGMethod),intent(in)::CG
		type(Tensor),intent(in)::P
		type(Tensor),external::GMsearchFunction
		integer::lenP,i
		lenP=P%getTotalData()
		call gradient%empty()
		call gradient%allocate((/lenP/),P%getType())
		do i=1,lenP
			call gradient%setValue(i,diff(CG,P,i,GMsearchFunction))
		end do
		return
	end function
!the search function is GMsearchFunction(inoutP,point_parameter)
!Search inoutP
	type(Tensor) function gradient_para(CG,P,point_parameter,GMsearchFunction)result(gradient)
		class(CGMethod),intent(in)::CG
		type(Tensor),intent(in)::P,point_parameter
		type(Tensor),external::GMsearchFunction
		integer::lenP,i
		lenP=P%getTotalData()
		call gradient%empty()
		call gradient%allocate((/lenP/),P%getType())
		do i=1,lenP
			call gradient%setValue(i,diff_para(CG,P,point_parameter,i,GMsearchFunction))
		end do
		return
	end function
! Gradient Method,search direction  
	type(Tensor) function Sdirection(CG,P,max_flag,GMsearchFunction)
		class(CGMethod),intent(in)::CG
		type(Tensor)::S
		type(Tensor),intent(in)::P
		logical,intent(in)::max_flag
		type(Tensor),external::GMsearchFunction
		real*8::Snorm
		S=gradient(CG,P,GMsearchFunction)
		Snorm=S%dnorm()
		if(Snorm.eq.0)then
			Sdirection=S
			write(*,*)"gradient is 0"
			return
		end if
		if(max_flag)then
			Sdirection=S/Snorm
		else
			Sdirection=S/((-1.)*Snorm)
		end if
		return
	end function
!the search function is GMsearchFunction(inoutP,point_parameter)
!Search inoutP
	type(Tensor) function Sdirection_para(CG,P,point_parameter,max_flag,GMsearchFunction)result(Sdirection)
		class(CGMethod),intent(in)::CG
		type(Tensor)::S
		type(Tensor),intent(in)::point_parameter
		type(Tensor),intent(in)::P
		logical,intent(in)::max_flag
		type(Tensor),external::GMsearchFunction
		real*8::Snorm
		S=gradient_para(CG,P,point_parameter,GMsearchFunction)
		Snorm=S%dnorm()
		if(Snorm.eq.0)then
			Sdirection=S
			write(*,*)"gradient is 0"
			return
		end if
		if(max_flag)then
			Sdirection=S/Snorm
		else
			Sdirection=S/((-1.)*Snorm)
		end if
		return
	end function
	subroutine GMDirection(CG,Sdir,P,max_flag,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(in)::CG
		type(Tensor),intent(inout)::Sdir
		logical,intent(in)::max_flag
		type(Tensor),intent(in)::P
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		if(present(GradientFunction)) then
			Sdir=GradientFunction(P)
			Sdir=Sdir/Sdir%norm()
			if(.not.max_flag) Sdir= (-1.)*Sdir
		else
			Sdir=Sdirection(CG,P,max_flag,GMsearchFunction)
		end if
		return
	end subroutine
	
	subroutine GMDirection_para(CG,Sdir,P,point_parameter,max_flag,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(in)::CG
		type(Tensor),intent(inout)::Sdir
		logical,intent(in)::max_flag
		type(Tensor),intent(in)::point_parameter
		type(Tensor),intent(in)::P
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		if(present(GradientFunction)) then
			Sdir=GradientFunction(P,point_parameter)
			Sdir=Sdir/Sdir%norm()
			if(.not.max_flag) Sdir= (-1.)*Sdir
		else
			Sdir=Sdirection_para(CG,P,point_parameter,max_flag,GMsearchFunction)
		end if
		return
	end subroutine
	
	type(Tensor) function go_on_one_step(CG,P,max_flag,step,GMsearchFunction,GradientFunction)result(oneStep)
		class(CGMethod),intent(in)::CG
		type(Tensor),intent(inout)::P
		logical,intent(in)::max_flag
		real*8,intent(inout)::step
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::f0,f1,P0,Sdir
		integer::i
		select case(CG%running_type)
			case (1)
				P0=P
				call GMDirection(CG,Sdir,P,max_flag,GMsearchFunction,GradientFunction)
				P=P0+(step*Sdir)
				f0=GMsearchFunction(P0)
				f1=GMsearchFunction(P)
				if(max_flag)then
					do while(f0.gt.f1) 
						step=step*CG%delta_step
						P=P0+(step*Sdir)
						f1=GMsearchFunction(P)
					end do
				else
					do while(f0.lt.f1) 
						step=step*CG%delta_step
						P=P0+(step*Sdir)
						f1=GMsearchFunction(P)
					end do
				end if
				oneStep=f1
			case (2)
				do i=1,CG%num_run_per_step
					call GMDirection(CG,Sdir,P,max_flag,GMsearchFunction,GradientFunction)
					P=P+(step*Sdir)
				end do
				oneStep=GMsearchFunction(P)
				step=step*CG%delta_step
		end select
		return
	end function
	
	type(Tensor) function go_on_one_step_para(CG,P,point_parameter,max_flag,step,GMsearchFunction,GradientFunction)result(oneStep)
		class(CGMethod),intent(in)::CG
		type(Tensor),intent(inout)::P
		type(Tensor),intent(in)::point_parameter
		logical,intent(in)::max_flag
		real*8,intent(inout)::step
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::f0,f1,P0,Sdir
		integer::i
		select case(CG%running_type)
			case (1)
				P0=P
				call GMDirection_para(CG,Sdir,P,point_parameter,max_flag,GMsearchFunction,GradientFunction)
				P=P0+(step*Sdir)
				f0=GMsearchFunction(P0,point_parameter)
				f1=GMsearchFunction(P,point_parameter)
				if(max_flag)then
					do while(f0.gt.f1) 
						step=step*CG%delta_step
						P=P0+(step*Sdir)
						f1=GMsearchFunction(P,point_parameter)
					end do
				else
					do while(f0.lt.f1) 
						step=step*CG%delta_step
						P=P0+(step*Sdir)
						f1=GMsearchFunction(P,point_parameter)
					end do
				end if
				oneStep=f1
			case (2)
				do i=1,CG%num_run_per_step
					call GMDirection_para(CG,Sdir,P,point_parameter,max_flag,GMsearchFunction,GradientFunction)
					P=P+(step*Sdir)
				end do
				oneStep=GMsearchFunction(P,point_parameter)
				step=step*CG%delta_step
		end select
		return
	end function
	
	
	
! 
	type(Tensor) function Gradient_Method1(CG,inoutP,max_flag,lenofinP,classtype,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(inout)::CG
		type(Tensor),intent(inout)::inoutP
		logical,intent(in)::max_flag
		integer,intent(in)::lenofinP
		character(len=*),intent(in)::classtype
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::P,res
		type(Tensor)::newf,f0,tempf
		real*8::step
		integer::i,templen,lenofP
		if(.not.inoutP%getFlag()) then
				lenofP=lenofinP
				P=generate((/lenofP/),classtype)
		else
				write(*,*)"ERROR in Gradient_Method"
				write(*,*)"one shoould input a empty Tensor"
				stop
		end if
		step=CG%max_step
		f0=GMsearchFunction(P)
		newf=go_on_one_step(CG,P,max_flag,step,GMsearchFunction,GradientFunction)
		i=1
		do while(dabs(newf-f0).gt. CG%stop_error)
			f0=newf
			newf=go_on_one_step(CG,P,max_flag,step,GMsearchFunction,GradientFunction)
			i=i+1
			if(i.gt.CG%max_running) then
				tempf=(newf-f0)/(newf+f0)
				call tempf%print('wornning,delta f/f:')
				exit
			end if
			if(CG%print_flag)then
				call writemess(newf+','+((newf-f0)/(newf+f0)))
			end if
		end do
		inoutP=P
		Gradient_Method1=newf
		return
	end function
	type(Tensor) function Gradient_Method2(CG,inoutP,max_flag,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(inout)::CG
		type(Tensor),intent(inout)::inoutP
		logical,intent(in)::max_flag
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::P,res
		real*8::step
		type(Tensor)::newf,f0
		integer::i,templen
		if(.not.inoutP%getFlag()) then
				write(*,*)"ERROR in Gradient_Method"
				write(*,*)"Specify the length of input point,or input a initial point"
				stop
		else
			P=inoutP
		end if
		step=CG%max_step
		f0=GMsearchFunction(P)
		newf=go_on_one_step(CG,P,max_flag,step,GMsearchFunction,GradientFunction)
		i=1
		do while(dabs(newf-f0).gt. CG%stop_error)
			f0=newf
			newf=go_on_one_step(CG,P,max_flag,step,GMsearchFunction,GradientFunction)
			i=i+1
			if(i.gt.CG%max_running) then
				call writemess('wornning,delta f/f,'+((newf-f0)/(newf+f0)))
				exit
			end if
			if(CG%print_flag)then
				call writemess(newf+','+((newf-f0)/(newf+f0)))
			end if
		end do
		inoutP=P
		Gradient_Method2=newf
		return
	end function
!the search function is GMsearchFunction(inoutP,point_parameter)
!Search inoutP
	type(Tensor) function Gradient_Method_para(CG,inoutP,point_parameter,max_flag,lenofinP,classtype,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(inout)::CG
		type(Tensor),intent(inout)::inoutP
		character(len=*),intent(in)::classtype
		type(Tensor),intent(in)::point_parameter
		logical,intent(in)::max_flag
		integer,intent(in)::lenofinP
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::P,res
		real*8::step
		type(Tensor)::newf,f0
		integer::i,templen,lenofP
		if(.not.inoutP%getFlag()) then
				lenofP=lenofinP
				P=generate((/lenofP/),classtype)
		else
				write(*,*)"ERROR in Gradient_Method"
				write(*,*)"one shoould input a empty Tensor"
				stop
		end if
		step=CG%max_step
		f0=GMsearchFunction(P,point_parameter)
		newf=go_on_one_step_para(CG,P,point_parameter,max_flag,step,GMsearchFunction,GradientFunction)
		i=1
		do while(dabs(newf-f0).gt. CG%stop_error)
			f0=newf
			newf=go_on_one_step_para(CG,P,point_parameter,max_flag,step,GMsearchFunction,GradientFunction)
			i=i+1
			if(i.gt.CG%max_running) then
				call writemess('wornning,delta f/f,'+((newf-f0)/(newf+f0)))
				exit
			end if
			if(CG%print_flag)then
				call writemess(newf+','+((newf-f0)/(newf+f0)))
			end if
		end do
		inoutP=P
		Gradient_Method_para=newf
		return
	end function
	type(Tensor) function Gradient_Method2_para(CG,inoutP,point_parameter,max_flag,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(inout)::CG
		Type(Tensor),intent(inout)::inoutP
		Type(Tensor),intent(in)::point_parameter
		logical,intent(in)::max_flag
		Type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		Type(Tensor)::P,res
		real*8::step
		real*8::newf,f0,tempnum
		integer::i,templen
		if(.not.inoutP%getFlag()) then
				write(*,*)"ERROR in Gradient_Method"
				write(*,*)"Specify the length of input point,or input a initial point"
				stop
		else
			P=inoutP
		end if
		step=CG%max_step
		f0=GMsearchFunction(P,point_parameter)
		newf=go_on_one_step_para(CG,P,point_parameter,max_flag,step,GMsearchFunction,GradientFunction)
		i=1
		do while(dabs(newf-f0).gt. CG%stop_error)
			f0=newf
			newf=go_on_one_step_para(CG,P,point_parameter,max_flag,step,GMsearchFunction,GradientFunction)
			i=i+1
			if(i.gt.CG%max_running) then
				call writemess('wornning,delta f/f,'+((newf-f0)/(newf+f0)))
				exit
			end if
			if(CG%print_flag)then
				call writemess(newf+','+((newf-f0)/(newf+f0)))
			end if
		end do
		inoutP=P
		Gradient_Method2_para=newf
		return
	end function
	
	
	subroutine step1(CG,inoutP,max_flag,step_,num_step,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(inout)::CG
		type(Tensor),intent(inout)::inoutP
		logical,intent(in)::max_flag
		integer,intent(in)::num_step
		class(*),intent(in)::step_
		real*8::step
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::f0,f1,P0,Sdir
		integer::i
		select type(step_)
			type is (integer)
				step=step_
			type is (real(kind=4))
				step=step_
			type is (real(kind=8))
				step=step_
			type is (complex(kind=4))
				step=step_
			type is (complex(kind=8))
				step=step_
		end select
		do i=1,num_step
				call GMDirection(CG,Sdir,inoutP,max_flag,GMsearchFunction,GradientFunction)
				inoutP=inoutP+(step*Sdir)
		end do
		return		
	end subroutine
	subroutine step2(CG,inoutP,max_flag,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(inout)::CG
		type(Tensor),intent(inout)::inoutP
		logical,intent(in)::max_flag
		integer::num_step
		real*8::step
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::f0,f1,P0,Sdir
		integer::i
		step=CG%max_step
		num_step=CG%max_running
		do i=1,num_step
				call GMDirection(CG,Sdir,inoutP,max_flag,GMsearchFunction,GradientFunction)
				inoutP=inoutP+(step*Sdir)
		end do
		return		
	end subroutine
	subroutine step3(CG,inoutP,max_flag,num_step,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(inout)::CG
		type(Tensor),intent(inout)::inoutP
		logical,intent(in)::max_flag
		integer,intent(in)::num_step
		real*8::step
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::f0,f1,P0,Sdir
		integer::i
		step=CG%max_step
		do i=1,num_step
				call GMDirection(CG,Sdir,inoutP,max_flag,GMsearchFunction,GradientFunction)
				inoutP=inoutP+(step*Sdir)
		end do
		return		
	end subroutine
	subroutine step1_para(CG,inoutP,point_parameter,max_flag,step_,num_step,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(inout)::CG
		type(Tensor),intent(inout)::inoutP
		logical,intent(in)::max_flag
		Type(Tensor),intent(in)::point_parameter
		integer,intent(in)::num_step
		class(*),intent(in)::step_
		real*8::step
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::f0,f1,P0,Sdir
		integer::i
		select type(step_)
			type is (integer)
				step=step_
			type is (real(kind=4))
				step=step_
			type is (real(kind=8))
				step=step_
			type is (complex(kind=4))
				step=step_
			type is (complex(kind=8))
				step=step_
		end select
		do i=1,num_step
				call GMDirection_para(CG,Sdir,inoutP,point_parameter,max_flag,GMsearchFunction,GradientFunction)
				inoutP=inoutP+(step*Sdir)
		end do
		return		
	end subroutine
	subroutine step2_para(CG,inoutP,point_parameter,max_flag,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(inout)::CG
		type(Tensor),intent(inout)::inoutP
		logical,intent(in)::max_flag
		Type(Tensor),intent(in)::point_parameter
		integer::num_step
		real*8::step
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::f0,f1,P0,Sdir
		integer::i
		step=CG%max_step
		num_step=CG%max_running
		do i=1,num_step
				call GMDirection_para(CG,Sdir,inoutP,point_parameter,max_flag,GMsearchFunction,GradientFunction)
				inoutP=inoutP+(step*Sdir)
		end do
		return		
	end subroutine
	subroutine step3_para(CG,inoutP,point_parameter,max_flag,num_step,GMsearchFunction,GradientFunction)
		class(CGMethod),intent(inout)::CG
		type(Tensor),intent(inout)::inoutP
		logical,intent(in)::max_flag
		Type(Tensor),intent(in)::point_parameter
		integer,intent(in)::num_step
		real*8::step
		type(Tensor),external::GMsearchFunction
		Type(Tensor),optional,external::GradientFunction
		type(Tensor)::f0,f1,P0,Sdir
		integer::i
		step=CG%max_step
		do i=1,num_step
				call GMDirection_para(CG,Sdir,inoutP,point_parameter,max_flag,GMsearchFunction,GradientFunction)
				inoutP=inoutP+(step*Sdir)
		end do
		return		
	end subroutine
	
end module	




















