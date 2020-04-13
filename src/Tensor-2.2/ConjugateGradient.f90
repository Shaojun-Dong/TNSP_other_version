!
!                   _ooOoo_
!                  o8888888o
!                  88" . "88
!                  (| -_- |)
!                  O\  =  /O
!               ____/`---'\____
!             .'  \\|     |//  `.
!            /  \\|||  :  |||//  \
!           /  _||||| -:- |||||-  \
!           |   | \\\  -  /// |   |
!           | \_|  ''\---/''  |   |
!           \  .-\__  `-`  ___/-. /
!         ___`. .'  /--.--\  `. . __
!      ."" '<  `.___\_<|>_/___.'  >'"".
!     | | :  `- \`.;`\ _ /`;.`/ - ` : | |
!     \  \ `-.   \_ __\ /__ _/   .-` /  /
!======`-.____`-.___\_____/___.-`____.-'======
!                   `=---='
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!       Buddha blessed , no BUG 
!****************************************************
!****************************************************
!*****************     note      ********************
!****************************************************
!****************************************************
!
!		1. This is the code of Conjugate Gradient Method. 
!  One should write the code of the function,say func(),
!  or the subroutine that output function value and the
!  Gradient if type(DTensor), to be maximum or minimum
!  when calling Conjugate_Gradient.
!
!		2. When calling Conjugate_Gradient without initialation,
!	the program will read the parameter as default setting. 
!	Or one use intConjugateGradient to initialation the program.
!	
!		3. The data in CGparameter are:
!----------------------------------------------------------------------------------------------
!stop_error                      1d-13           !if abs(f(new_x)-f(x))<stop_error,output
!max_running                     300             !max running
!diff_delta                      1d-8            !differentiation,(f(x+delta)-f(x))/delta,diff_delta is the delta,If input Gradient function, it is useless
!LinearSearch_max_running        4               !it is not good  when too large,LinearSearch_max_running <7 is OK
!p_flag                          0               !if 1,print the progress
!max_step                        10.             !go along the dirction with max_step,x=max_step
!delta_step                      0.9             !if f(new_x)<f(x),max_step=max_step*delta_step,delta_step should smaller than 1 
!directionFlag                    3              !0:Gradient Method,1:g_{i+1}*g_i,2:g_{i+1}*(g_{i+1}-g_i)3: g_{i+1}*(g_{i+1}-g_i) or g_{i+1}*g_i
!num_same_step						  10             !See the note below
!min_step_when_num_same_step     5d-2				!See the note below
!run_type									1              !See the note below
!stop_type								2              !stop_type=1, stop_error will be |Gradient|^2,stop_type=2, stop_error will be abs(f(new_x)-f(x))
!num_run_type                     3            !It is useful when run_type=6,see the note of run_type=6
!step_run_type                   100 200       !It is useful when run_type=6,see the note of run_type=6
!which_type                       1 2 3        !It is useful when run_type=6,see the note of run_type=6
!----------------------------------------------------------------------------------------------
!   NOTE of PARAMETER
!run_type=1
!		use LinearSearch along the CG direction of GM direction, LinearSearch_max_running should be larger than 0
!		at the first time  max_step is the randomt tail point in LinearSearch
!run_type=2	
!		do not use LinearSearch. Go along the direction with length x :newP = P + x * direction
!	  when f(newP)>f(P) (find min), x=x*delta_step , do again
!	  at the first time x=max_step
!run_type=3
!		do not use LinearSearch. Go along the direction with length x :newP = P + x * direction
!	  do newP = P + x * direction for num_same_step times, and then x=x*delta_step. num_same_step
!	  should be larger than 0.
!run_type=4	  
!		the same as run_type=3,but  when x<min_step_when_num_same_step , set run_type=2,x=max_step
!	  search again
!run_type=5	  
!		the same as run_type=3,but  when x<min_step_when_num_same_step , set run_type=1,x=0d0
!	  search again	 
!
! directionFlag: if directionFlag=0 use the Gradient Method
!                if directionFlag=1 use the CG Method with search direction of :g_{i+1}*g_i/norm2(g_i)
!                if directionFlag=2 use the CG Method with search direction of :g_{i+1}*(g_{i+1}-g_i)/norm2(g_i)
!                if directionFlag=3 use the CG Method with search direction of :g_{i+1}*g_i/norm2(g_i) ,if its vlaue < zero_number
!                    change it to g_{i+1}*(g_{i+1}-g_i)/norm2(g_i)
!
!run_type=6	  
!		using the type of which_type,suppose num_run_type=N, then will read N-1 value of step_run_type and N value which_type.
!		suppse they are
!			num_run_type                     3            !It is useful when run_type=6
!			step_run_type                   100 200    !It is useful when run_type=6
!			which_type                       3 1 2        !It is useful when run_type=6		
!		then it means:
!		      do i=1 to 100, runtype=3
!		      and i=101 to 200 , runtype=1
!		      at lase ,i=201 to the end, runtype=2
!		if runtype=4 or 5 .the i will count from the sitiation when x<min_step_when_num_same_step
!		example
!			num_run_type                     2            !It is useful when run_type=6
!			step_run_type                   100    !It is useful when run_type=6
!			which_type                       4 1        !It is useful when run_type=6	
!		then it means:
!		      do runtype=4, until x<min_step_when_num_same_step, and change to run_type=2 doing i=1 to 100 
!		      at lase ,i=201 to the end, runtype=1		      
!
!		4. Send email to tyrants@qq.com to report any bugs.
!
!****************************************************
!****************************************************
module ConjugateGradientMethod
	use eigen_value
	implicit none
	real*8,parameter::zero_number=1d-16
	real*8,parameter::gradient_zero_number=5d-8!To small will error
	logical,private::print_flag=.false.
	real*8,private::diff_delta=1d-8!If input Gradient function, it is useless
	real*8,private::stop_error=1d-15
	integer,private::max_running=300
	integer,private::LinearSearch_max_running=3
	real*8,private::max_step=5d0
	real*8,private::delta_step=0.9
	integer,private::num_same_step=0
	real*8,private::min_step_when_num_same_step=1d-3
	integer,private::run_type=1
	integer,private::stop_type=1
	integer,private::num_run_type=0
	integer,allocatable,private::step_run_type(:)
	integer,allocatable,private::which_type(:)
	
	real*8,private::InfinityNumber=1d100
	!In LinearSearch_fouth_point, if case of two point are too closed 
	
!
!
	
	integer,private::directionFlag=1
	
	!0:  Gradient Method
	!1: g_{i+1}*g_i
	!2:g_{i+1}*(g_{i+1}-g_i)
	real*8,private::direction_corr=0d0
	interface Conjugate_Gradient
		module procedure CG_Method1
		module procedure CG_Method1_inputP
		module procedure CG_Method2
		module procedure CG_Method2_inputP
	end interface
!		real*8 a=Conjugate_Gradient(inoutP,lenofinP,GradSubroutine,maxflag,point_parameter)
!		real*8 a=Conjugate_Gradient(inoutP,GradSubroutine,maxflag,point_parameter)
!		real*8 a=Conjugate_Gradient(inoutP,lenofinP,searchFunction,maxflag,point_parameter)
!		real*8 a=Conjugate_Gradient(inoutP,searchFunction,maxflag,point_parameter)
!     where 
!
!		 inoutP : is the initial point,type(DTensor) for searching the max/min
!	               if inoutP is no data, one should set lenofinP as the
!		            length of inoutP, the program will create random data
!
!		 lenofinP : It is the length of inoutP,if this are no data in inoutP
!
!		 GradSubroutine: The subroutine that output function value and the Gradient
!
!                      call GradSubroutine(funcValue,grad_value,point,graflag,point_parameter)
!
!		                 and funcValue:real*8 the value of the searchFunction on point
!		                     grad_value:type(DTensor),the Gradient of the searchFunction on point
!		                     point: the variable to be fix to min or max searchFunction
!                         grafalg: .true. output funcValue and grad_value else output funcValue only
!                          point_parameter:optional, the parameter use in the searchFunction

!		 searchFunction: The function to be max or min
!                      real*8 f=searchFunction(point),or f=searchFunction(point,point_parameter
!                      the gradient is of searchFunction will be
!		                      (searchFunction(point+Delta_P)-searchFunction(point-Delta_P))/ (2*Delta_P)
!
!		 max_flag=.true.,find the max,if not ,min
!			
!		 point_parameter: optional, the parameter use in the searchFunction
!
	
	interface intConjugateGradient
		module procedure intConjugateGradient1
		module procedure intConjugateGradient2
		module procedure intConjugateGradient3
	end interface
!		intConjugateGradient(stop_error_,max_running_,diff_delta_,LinearSearch_max_running_,p_flag)
!  or intConjugateGradient(stop_error_,max_running_,,LinearSearch_max_running_,p_flag)
!	or	intConjugateGradient(CGparameter_address),the input file can copy from Tensor/MPI/CGparameter

!**********************************************
!		Example of searchFunction
!**********************************************
!
!	real*8 function searchFunction(point)
!		type(DTensor),intent(in)::point
!
!.......  write the function here .......
!		
!		return
!	end function
!---------------------------------------------
!or 
!	real*8 function searchFunction(point,para)!The code will not search on para
!		type(DTensor),intent(in)::point,para
!
!.......  write the function here .......
!		
!		return
!	end function
!**********************************************
!**********************************************
!**********************************************
!		Example of GradSubroutine
!**********************************************
!
!	subroutine Example_GradSubroutine(func,grad,point,flag)
!		real*8,intent(out)::func
!		type(DTensor),intent(out)::grad
!		type(DTensor),intent(in)::point
!     if(flag) then
!
!   .......  write the function output grad here .......
!
!     end if
!
!   .......  write the function output func here .......
!		
!		return
!	end subroutine
!---------------------------------------------
!or 
!	subroutine Example_GradSubroutine(func,grad,point,para)
!		real*8,intent(out)::func
!		type(DTensor),intent(out)::grad
!		type(DTensor),intent(in)::point,para
!     if(flag) then
!
!   .......  write the function output grad here .......
!
!     end if
!
!   .......  write the function output func here .......
!		
!		return
!	end subroutine
!**********************************************
!**********************************************

!example program
!	real*8 maxpoint
!	type(DTensor)::point
!  integer::lenP !lenP is the length of point,the code will random create a DTensor of length of lenP
!	maxpoint=Gradient_Method(point,lenP,.true.,Example_searchFunction) 
!or
!	maxpoint=Gradient_Method(point,lenP,.true.,Example_searchFunction,para)
!or
!	maxpoint=Gradient_Method(point,lenP,.true.,Example_GradSubroutine,para)
!or
!	maxpoint=Gradient_Method(point,lenP,.true.,Example_GradSubroutine,para)
!---------------------------------------------------------

!**********************************************
contains
	subroutine intConjugateGradient1(CGparameter_address)
		CHARACTER(len=*),intent(in)::CGparameter_address
		CHARACTER*100::notused
		integer::p_flag
		integer::i
		open(unit=10,file=CGparameter_address,STATUS='old')
		read(10,*)notused,stop_error!if abs(f(new_x)-f(x))<stop_error,output
		read(10,*)notused,max_running!max running
		read(10,*)notused,diff_delta!differentiation,(f(x+delta)-f(x))/delta,diff_delta is the delta
		                            !If input Gradient function, it is useless
		read(10,*)notused,LinearSearch_max_running!it is not good  when too large,LinearSearch_max_running <7 is OK
		read(10,*)notused,p_flag!if 1,print the progress
		if(p_flag.eq.1)then
			print_flag=.true.
		end if
		read(10,*)notused,max_step !go along the dirction with max_step,new_x=x+max_step
		read(10,*)notused,delta_step!if f(new_x)<f(x),max_step=max_step*delta_step,delta_step should smaller than 1 
		read(10,*)notused,directionFlag
		read(10,*)notused,num_same_step
		read(10,*)notused,min_step_when_num_same_step
		read(10,*)notused,run_type
		read(10,*)notused,stop_type
		if(run_type.eq.6)then
			read(10,*)notused,num_run_type
			if(num_run_type.le.1)then
				write(*,*)"ERROR in runtype=6"
				stop
			end if
			allocate(step_run_type(num_run_type))
			read(10,*)notused,(step_run_type(i),i=1,num_run_type-1)
			step_run_type(num_run_type)=-1
			allocate(which_type (num_run_type))
			read(10,*)notused,(which_type(i),i=1,num_run_type)
		end if
		close(unit=10)
		return
	end subroutine
	subroutine intConjugateGradient2(stop_error_,max_running_,diff_delta_,LinearSearch_max_running_,p_flag,&
														max_step_,delta_step_,directionFlag_)
		real*8,intent(in)::stop_error_,diff_delta_,max_step_,delta_step_
		integer,intent(in)::max_running_,LinearSearch_max_running_,directionFlag_
		logical,intent(in)::p_flag
		stop_error=stop_error_
		diff_delta=diff_delta_
		max_running=max_running_
		LinearSearch_max_running=LinearSearch_max_running_
		print_flag=p_flag
		max_step=max_step_
		delta_step=delta_step_
		directionFlag=directionFlag_
		return
	end subroutine
	subroutine intConjugateGradient3(stop_error_,max_running_,LinearSearch_max_running_,p_flag,&
							max_step_,delta_step_,directionFlag_)
		real*8,intent(in)::stop_error_,max_step_,delta_step_
		integer,intent(in)::max_running_,LinearSearch_max_running_,directionFlag_
		logical,intent(in)::p_flag
		stop_error=stop_error_
		max_running=max_running_
		LinearSearch_max_running=LinearSearch_max_running_
		print_flag=p_flag
		max_step=max_step_
		delta_step=delta_step_
		directionFlag=directionFlag_
		return
	end subroutine
	subroutine set_run_type(rtype,x)
		integer,intent(in)::rtype
		real*8,intent(inout)::x
		run_type=rtype
		if(run_type.eq.1) then
			x=0d0
		else
			x=max_step
		end if
		if(print_flag) write(*,*)"set run_type as" ,run_type
		return
	end subroutine
	subroutine set_stop_error(error)
		real*8,intent(in)::error
		stop_error=error
		return
	end subroutine
!*************************************************************
!  CG method: f(P) is the function to be min
!             g(P) is the 	gradient of f(P)
!	1. Find the direction, d1=-g(P1)
!  2. search x, that min f(P1+x*d1)
!    1) Linear search, input random x0
!       and output P2, and x1
!  3. Find the CG direction d_i=-g(P_i)+\lambda*d_{i-1},where 
!                 g(P_i)^T*g(P_{i-1})
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
			outx=max_step
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
	real*8 function ConjugateGradient_diff(inputP,ith,SearchFunction,point_parameter)result(diff)
		type(DTensor),intent(in)::inputP
		integer,intent(in)::ith
		real*8,external::SearchFunction
		type(DTensor),optional,intent(inout)::point_parameter
		type(Dtensor)::delta_P,P2,P1
		integer::lenP
		lenP=DgetTotaldata(inputP)
		call allocateDTensor(delta_P,(/lenP/))
		delta_P=0d0
		call Dmodify(delta_P,ith,diff_delta)
		P2=inputP+delta_P
		P1=inputP-delta_P
		if(present(point_parameter)) then
			diff=SearchFunction(P2,point_parameter)-SearchFunction(P1,point_parameter)
		else
			diff=SearchFunction(P2)-SearchFunction(P1)
		end if
		diff=diff/(diff_delta+diff_delta)
		if(isnan(diff)) then
			write(*,*)"NAN error,diff"
			!call DTMprint(P1)
			write(*,*)SearchFunction(P1)
			stop
		end if
		return
	end function
! gradient vector of searchFunction
	type(DTensor) function ConjugateGradient_gradientFunction(P,SearchFunction,point_parameter)result(gradientFunction)
		type(DTensor),intent(in)::P
		real*8,external::SearchFunction
		type(DTensor),optional,intent(inout)::point_parameter
		integer::lenP,i
		real*8::Pdatai
		lenP=DgetTotaldata(P)
		call allocateDTensor(gradientFunction,(/lenP/))
		do i=1,lenP
			Pdatai=ConjugateGradient_diff(P,i,SearchFunction,point_parameter)
			call Dmodify(gradientFunction,i,Pdatai)
		end do
		return
	end function	

!	GradSubroutine is a subroutine
!	call GradSubroutine(func,gradient,P,point_parameter)
!	output func,which is the function to be min or max
!	and gradient is the gradient of func, it is a vector
!  df(x)/dx= vec{gradient(p+x*dir) } * vec{dir}
!
!	input 
!	P  : the variable to be fix to min or max searchfunction
!	x	: the trial lenth of linear search,use the value of prior
!	    step, if x=0, x will be x=max_step
!	dir: The CG search direction
!	gra: the gradient of the search function on P
!
!	output
!	P  : The new P find from Linear search
!	x	: the new x, use for next step
!	f  : The value of searchfunction on new P
!	gra: the gradient of the search function on new P
!
!note on linear search, the variable to fix is x, in sf(P+x*dir)=searchfunction(P+x*dir)
!	1.gradient of  linear searchfunction g is
!	    g=d sf(P+x*dir) /d x= \vec{gf(P+x*dir)} * vec{dir}
!	    where \vec{gf(P+x*dir)} is the gradient of searchfunction
!	2. The first step, set x=0
	subroutine LinearSearch1(P,x,f,gra,dir,GradSubroutine,maxflag,point_parameter)
		type(DTensor),intent(inout)::P,gra
		type(DTensor),intent(in)::dir
		type(DTensor),optional,intent(inout)::point_parameter
		real*8,intent(inout)::x,f
		logical,intent(in)::maxflag
		real*8::x1,x2,x3,outx
		logical::stopflag
		external::GradSubroutine
		type(DTensor)::NewP
		real*8::f1,f2,f3,g1,g2,g3
		integer::i
		g1=gra.x.dir
		f1=f
		x1=0d0
		
		if(dabs(x).le.zero_number) then
			x2=max_step
		else 
			x2=x
		end if
		
		NewP=P+(x2*dir)
		if(present(point_parameter)) then
			call GradSubroutine(f2,gra,NewP,.true.,point_parameter)
		else
			call GradSubroutine(f2,gra,NewP,.true.)
		end if
		g2=gra.x.dir
		if(dabs(g2).le.gradient_zero_number)then
			P=NewP
			x=x2
			f=f2
			return
		end if
		
		call LinearSearch_third_point(x3,f1,g1,f2,g2,x2,maxflag)
		NewP=P+(x3*dir)
		if(present(point_parameter)) then
			call GradSubroutine(f3,gra,NewP,.true.,point_parameter)
		else
			call GradSubroutine(f3,gra,NewP,.true.)
		end if
		g3=gra.x.dir
		if(dabs(g3).le.gradient_zero_number)then
			P=NewP
			x=x3
			f=f3
			return
		end if
		stopflag=.false.
		do i=1,LinearSearch_max_running
			call LinearSearch_fouth_point(outx,stopflag,x1,f1,g1,x2,f2,g2,x3,f3,g3,maxflag)
			if(isnan(outx).or.(outx.lt.0d0))then
				write(*,*)"ERROR in LinearSearch"
				write(*,*)"x1,x2"
				write(*,*)x1
				write(*,*)x2
				stop
			end if
			NewP=P+(outx*dir)
			if(present(point_parameter)) then
				call GradSubroutine(f3,gra,NewP,.true.,point_parameter)
			else
				call GradSubroutine(f3,gra,NewP,.true.)
			end if
			g3=gra.x.dir
			x3=outx
			if((dabs(g3).le.gradient_zero_number).or. stopflag) exit
		end do
		P=NewP
		x=outx
		f=f3
		return
	end subroutine	
!  searchFunction is the function to be min or max
!  real*8 function searchFunction(P,point_parameter)
!  and the gradient of searchFunction will be 
!		 (searchFunction(P+Delta_P)-searchFunction(P-Delta_P))/ (2*Delta_P)

!	(P,x,f,gra,dir,GradSubroutine,maxflag,point_parameter)
	subroutine LinearSearch2(P,x,f,gra,dir,searchFunction,maxflag,point_parameter)
		type(DTensor),intent(inout)::P,gra
		type(DTensor),intent(in)::dir
		type(DTensor),optional,intent(inout)::point_parameter
		real*8,external::searchFunction
		real*8,intent(inout)::x,f
		logical,intent(in)::maxflag
		real*8::x1,x2,x3,outx
		type(DTensor)::NewP
		logical::stopflag
		real*8::f1,f2,f3,g1,g2,g3
		integer::i
		g1=gra.x.dir
		f1=f
		x1=0d0
		
		if(dabs(x).le.zero_number) then
			 x2=max_step
		else 
			x2=x
		end if
		NewP=P+(x2*dir)
		if(present(point_parameter)) then
			f2=searchFunction(NewP,point_parameter)
		else
			f2=searchFunction(NewP)
		end if
		gra=ConjugateGradient_gradientFunction(NewP,SearchFunction,point_parameter)
		g2=gra.x.dir
		if(dabs(g2).le.gradient_zero_number)then
			P=NewP
			x=x2
			f=f2
			return
		end if
		call LinearSearch_third_point(x3,f1,g1,f2,g2,x2,maxflag)
		NewP=P+(x3*dir)
		if(present(point_parameter)) then
			f3=searchFunction(NewP,point_parameter)
		else
			f3=searchFunction(NewP)
		end if
		gra=ConjugateGradient_gradientFunction(NewP,SearchFunction,point_parameter)
		g3=gra.x.dir
		if(dabs(g3).le.gradient_zero_number)then
			P=NewP
			x=x3
			f=f3
			return
		end if
		stopflag=.false.
		do i=1,LinearSearch_max_running
			call LinearSearch_fouth_point(outx,stopflag,x1,f1,g1,x2,f2,g2,x3,f3,g3,maxflag)
			if(isnan(outx).or.(outx.lt.0d0))then
				write(*,*)"ERROR in LinearSearch"
				write(*,*)"x1,x2"
				write(*,*)x1
				write(*,*)x2
				stop
			end if
			NewP=P+(outx*dir)
			if(present(point_parameter)) then
				f3=searchFunction(NewP,point_parameter)
			else
				f3=searchFunction(NewP)
			end if
			gra=ConjugateGradient_gradientFunction(NewP,SearchFunction,point_parameter)
			g3=gra.x.dir
			x3=outx
			if( (dabs(g3).le.gradient_zero_number) .or. stopflag) exit
		end do
		P=NewP
		x=outx
		f=f3
		return
	end subroutine	
!*************************************************************
!    Linearstep
!
!   
!	input 
!	P  : the variable to be fix to min or max searchfunction
!	x	: the trial lenth of linear search,use the value of prior
!	    step, if x=0, x will be max_step
!	dir: The CG search direction
!	gra: the gradient of the search function on P
!
!	output
!	P  : The new P find from Linear search
!	x	: the new x, use for next step
!	f  : The value of searchfunction on new P
!	gra: the gradient of the search function on new P
	subroutine Linearstep1(P,x,f,gra,dir,GradSubroutine,maxflag,point_parameter)
		type(DTensor),intent(inout)::P,gra
		type(DTensor),intent(in)::dir
		type(DTensor),optional,intent(inout)::point_parameter
		external::GradSubroutine
		real*8,intent(inout)::x,f
		logical,intent(in)::maxflag
		real*8::oldf
		type(DTensor)::OldP
		if(dabs(x).le.zero_number) x=max_step
		oldf=f
		OldP=P
		P=OldP+(x*dir)
		if(run_type.ne.2) then
			if(present(point_parameter)) then
				call GradSubroutine(f,gra,P,.true.,point_parameter)
			else
				call GradSubroutine(f,gra,P,.true.)
			end if
			return
		else
			if(present(point_parameter)) then
				call GradSubroutine(f,gra,P,.false.,point_parameter)
			else
				call GradSubroutine(f,gra,P,.false.)
			end if
		end if
		
		if(maxflag)then
			do while(oldf.ge.f) 
				x=x*delta_step
				P=OldP+(x*dir)
				if(present(point_parameter)) then
					call GradSubroutine(f,gra,P,.false.,point_parameter)
				else
					call GradSubroutine(f,gra,P,.false.)
				end if
			end do
		else
			do while(oldf.le.f) 
				x=x*delta_step
				P=OldP+(x*dir)
				if(present(point_parameter)) then
					call GradSubroutine(f,gra,P,.false.,point_parameter)
				else
					call GradSubroutine(f,gra,P,.false.)
				end if
			end do
		end if 
		if(present(point_parameter)) then
			call GradSubroutine(f,gra,P,.true.,point_parameter)
		else
			call GradSubroutine(f,gra,P,.true.)
		end if
		return
	end subroutine
	subroutine Linearstep2(P,x,f,gra,dir,searchFunction,maxflag,point_parameter)
		type(DTensor),intent(inout)::P,gra
		type(DTensor),intent(in)::dir
		type(DTensor),optional,intent(inout)::point_parameter
		real*8,external::searchFunction
		real*8,intent(inout)::x,f
		logical,intent(in)::maxflag
		real*8::oldf
		type(DTensor)::OldP
		if(dabs(x).le.zero_number) x=max_step
		oldf=f
		OldP=P
		P=OldP+(x*dir)
		if(present(point_parameter)) then
			f=searchFunction(P,point_parameter)
		else
			f=searchFunction(P)
		end if
		if(run_type.ne.2)then!do not modify x
			gra=ConjugateGradient_gradientFunction(P,SearchFunction,point_parameter)
			return
		end if
		if(maxflag)then
			do while(oldf.gt.f) 
				x=x*delta_step
				P=OldP+(x*dir)
				if(present(point_parameter)) then
					f=searchFunction(P,point_parameter)
				else
					f=searchFunction(P)
				end if
			end do
		else
			do while(oldf.lt.f) 
				x=x*delta_step
				P=OldP+(x*dir)
				if(present(point_parameter)) then
					f=searchFunction(P,point_parameter)
				else
					f=searchFunction(P)
				end if
			end do
		end if 
		gra=ConjugateGradient_gradientFunction(P,SearchFunction,point_parameter)
		return
	end subroutine
!*************************************************************
!*************************************************************
!                 CG direction
!*************************************************************
!*************************************************************

!
!
!	P is the value to search in
!	input P
!  input gradient(Gra), and direction(dir) of prior step
!
!  output the gradient of P, store in Gra
!	direction of P, store in dir
!	The search function value, store in f
!	the gradient use in linear search, store in g
!
!	at the first step, no dir , priorgra. they are empty STensor
!
!	input
!	P  : the variable to be fix to min or max searchfunction
!	dir: The CG search direction of prior step
!	priorgra: the gradient of the search function of prior step
!  gra: the gradient of the search function of P (This step)
!
!	output
!	dir: The CG search direction on P
!  priorgra: It will be rewrite be the the gradient of P of this step
	subroutine Searchdirection(dir,Gra,priorgra,maxflag)
		type(DTensor),intent(in)::Gra
		type(DTensor),intent(inout)::dir,priorgra
		logical,intent(in)::maxflag
		real*8::dirnorm
		if(DgetFlag(dir)) then
			direction_corr=correctPara(gra,priorGra)
			if(maxflag)then
				dir=gra+(direction_corr*dir)
			else
				dir=(direction_corr*dir)-gra
			end if
		else
			if(maxflag)then
				dir=gra
			else
				dir=(-1d0)*gra
			end if
		end if
		dirnorm=dnorm(dir)
		if(dirnorm.le.1d-16)then
			priorgra=Gra
			return
		end if
		dir=dir/dirnorm
		priorgra=Gra
		return
	end subroutine
	
	real*8 function correctPara(newgra,gra)
		type(DTensor),intent(in)::newgra,gra
		real*8::nor
		if(directionFlag.eq.0)then
			correctPara=0d0
			return
		end if
		if(directionFlag.eq.1)then
			correctPara=newgra.x.gra
			correctPara=correctPara/dnorm2(gra)
			if(isnan(correctPara))then
				call DTMprint(gra)
				write(*,*)dnorm2(gra),newgra.x.gra
				stop
			end if
			return
		end if
		if(directionFlag.eq.2)then
			correctPara=newgra.x.(newgra-gra)
			correctPara=correctPara/dnorm2(gra)
			return
		end if
		nor=1d0/dnorm2(gra)
		correctPara=newgra.x.gra
		correctPara=correctPara*nor
		if(correctPara.le.0.5) then
			correctPara=(dnorm2(newgra)*nor)-correctPara
		end if
		return
	end function
	
!*************************************************************
!*************************************************************
!                step 
!*************************************************************
!*************************************************************		

!	P is the value to search in
!  input gradient(Gra), and direction(dir) of prior step
!	x is the parameter to fix at linear search, the value of prior step
!  
!
!	Then the code do :
!		1. find the CG direction of the input P.
!		2. Do the linear search, find the min(or max) point along this direction
!
!	input
!	P  : the variable to be fix to min or max searchfunction
!	dir: The CG search direction of prior step
!	priorgra: the gradient of the search function of prior step 
!  gra: the gradient of the search function of this step
!	x	: the trial lenth of linear search,use the value of prior
!	    step, if x=0, x will be random value
!
!	on output
!	P  : the new point
!	f  : the min search function of new P
!	f0 : the search function of old P
!	dir: the CG direction of new P
!	Gra: the gradient of new P
!	x  : the parameter to fix at linear search, which is use in the next step as the first input

!	at the first step, no dir , Gra. they are empty STensor, and set x=0d0
	subroutine go_one_step1(P,f,f0,Gra,priorgra,dir,x,GradSubroutine,maxflag,point_parameter)
		type(DTensor),intent(inout)::P,dir,Gra,priorgra
		real*8,intent(inout)::f,x,f0
		external::GradSubroutine
		logical,intent(in)::maxflag
		type(DTensor),optional,intent(inout)::point_parameter
		real*8::g
		if(.not.DgetFlag(Gra))then
			if(present(point_parameter)) then
				call GradSubroutine(f0,gra,P,.true.,point_parameter)
			else
				call GradSubroutine(f0,gra,P,.true.)
			end if
		end if
		call  Searchdirection(dir,Gra,priorgra,maxflag)
		f=f0
		if(run_type.eq.1)then
			call LinearSearch1(P,x,f,Gra,dir,GradSubroutine,maxflag,point_parameter)
		else
			call Linearstep1(P,x,f,gra,dir,GradSubroutine,maxflag,point_parameter)
			if((dabs(x).le.zero_number).and.(directionFlag.ne.0)) then!The CG direction can not find any point usefull, use  Gradient Method
				x=max_step
				if(maxflag)then
					dir=priorgra
				else
					dir=(-1d0)*priorgra
				end if
				direction_corr=0d0
				dir=dir/dnorm(dir)
				call Linearstep1(P,x,f,gra,dir,GradSubroutine,maxflag,point_parameter)
			end if
		end if
		return
	end subroutine
	subroutine go_one_step2(P,f,f0,Gra,priorgra,dir,x,searchFunction,maxflag,point_parameter)
		type(DTensor),intent(inout)::P,dir,Gra,priorgra
		real*8,intent(inout)::f,x,f0
		real*8,external::searchFunction
		logical,intent(in)::maxflag
		type(DTensor),optional,intent(inout)::point_parameter
		real*8::g
		if(.not.DgetFlag(Gra))then
			if(present(point_parameter)) then
				f0=searchFunction(P,point_parameter)
			else
				f0=searchFunction(P)
			end if	
			gra=ConjugateGradient_gradientFunction(P,SearchFunction,point_parameter)
		end if
		call  Searchdirection(dir,Gra,priorgra,maxflag)
		f=f0
		if(run_type.eq.1)then
			call LinearSearch2(P,x,f,Gra,dir,searchFunction,maxflag,point_parameter)
		else
			call Linearstep2(P,x,f,gra,dir,searchFunction,maxflag,point_parameter)
			if((dabs(x).le.zero_number).and.(directionFlag.ne.0)) then!The CG direction can not find any point usefull, use  Gradient Method
				x=max_step
				if(maxflag)then
					dir=priorgra
				else
					dir=(-1d0)*priorgra
				end if
				direction_corr=0d0
				dir=dir/dnorm(dir)
				call Linearstep2(P,x,f,gra,dir,searchFunction,maxflag,point_parameter)
			end if
		end if
		return
	end subroutine
		
!*************************************************************
!*************************************************************
!                 Conjugate_Gradient
!*************************************************************
!*************************************************************		
	real*8 function CG_Method1(inoutP,lenofinP,GradSubroutine,maxflag,point_parameter)
		type(DTensor),intent(inout)::inoutP
		logical,intent(in)::maxflag
		integer,intent(in)::lenofinP
		external::GradSubroutine
		type(DTensor),optional,intent(inout)::point_parameter
		type(DTensor)::dir,Gra,priorgra
		real*8::newf,f0,x,stop_value
		integer::i,lenofP,runtyp6counter,typecounter
		logical::runtyp6flag
		if(.not.DgetFlag(inoutP)) then
				lenofP=lenofinP
				inoutP=Dgenerate((/lenofP/),(/-1d0,1d0/))
		else
				write(*,*)"ERROR in CG_Method"
				write(*,*)"one shoould input a empty Tensor"
				stop
		end if
		i=0
		x=0d0
		call emptyDTensor(dir)
		call emptyDTensor(Gra)
		call go_one_step1(inoutP,newf,f0,dir,Gra,priorgra,x,GradSubroutine,maxflag,point_parameter)
		if(stop_type.eq.1)then
			stop_value=dnorm2(Gra)
		else if(stop_type.eq.2)then
			stop_value=dabs((newf-f0)/(newf+f0))
		else
			write(*,*)"ERROR stop type"
		end if
		if(print_flag)then
			if(stop_type.eq.1)then
				write(*,*)"newf,|gradient|^2,len of step(x),direction_corr"
			else if(stop_type.eq.2)then
				write(*,*)"newf,dabs((newf-f0)/(newf+f0)),len of step(x),direction_corr"
			end if
		end if
		runtyp6flag=run_type.eq.6
		runtyp6counter=0
		typecounter=1
		if(runtyp6flag) call set_run_type(which_type(typecounter),x) 
		
		do while(stop_value.gt.stop_error)
			f0=newf
			i=i+1
			call go_one_step1(inoutP,newf,f0,dir,Gra,priorgra,x,GradSubroutine,maxflag,point_parameter)
			if(stop_type.eq.1)then
				stop_value=dnorm2(Gra)
			else if(stop_type.eq.2)then
				stop_value=dabs((newf-f0)/(newf+f0))
			end if
			if( (run_type.eq.3) .or. (run_type.eq.4 ) .or. (run_type.eq.5 )) then
				if(mod(i,num_same_step).eq.0) then
					x=x*delta_step
				end if
			end if
			if( (run_type.eq.4 ) .or. (run_type.eq.5 )) then
				if(x.lt.min_step_when_num_same_step)then
					max_step=min_step_when_num_same_step
					if(run_type.eq.4) then
						call set_run_type(2,x)
					else
						call set_run_type(1,x)
					end if
				end if
			end if
			if( (run_type.eq.3 ) ) then
				if(x.lt.min_step_when_num_same_step)then
					if(stop_type.eq.1) then
						write(*,*)"wornning,|gradient|^2",stop_value
					else if(stop_type.eq.2) then
						write(*,*)"wornning,dabs((newf-f0)/(newf+f0))",stop_value
					end if
					exit
				end if
			end if
			if(runtyp6flag)then
				if((run_type.eq.1) .or. (run_type.eq.2 ) .or. (run_type.eq.3 ))then
					runtyp6counter=runtyp6counter+1
					if(runtyp6counter.eq.step_run_type(typecounter))then
						typecounter=typecounter+1
						call set_run_type(which_type(typecounter),x) 
					end if
				end if
			end if
			if(print_flag)then
				write(*,*)newf,stop_value,x,direction_corr
			end if
			if(i.gt.max_running) then
				if(stop_type.eq.1)then
					write(*,*)"wornning,|gradient|^2",stop_value
				else if(stop_type.eq.2)then
					write(*,*)"wornning,dabs((newf-f0)/(newf+f0))",stop_value
				end if
				exit
			end if
		end do
		CG_Method1=newf
		if(print_flag)then
			write(*,*)"num_of_step_to_finish_search",i
		end if
		return
	end function
	real*8 function CG_Method1_inputP(inoutP,GradSubroutine,maxflag,point_parameter)
		type(DTensor),intent(inout)::inoutP
		logical,intent(in)::maxflag
		external::GradSubroutine
		type(DTensor),optional,intent(inout)::point_parameter
		type(DTensor)::dir,Gra,priorgra
		real*8::newf,f0,x,stop_value
		integer::i,lenofP,runtyp6counter,typecounter
		logical::runtyp6flag
		if(.not.DgetFlag(inoutP)) then
				write(*,*)"ERROR in CG_Method"
				write(*,*)"Specify the length of input point,or input a initial point"
				stop
		else
			lenofP=DgettotalData(inoutP)
		end if
		i=0
		x=0d0
		call emptyDTensor(dir)
		call emptyDTensor(Gra)
		call go_one_step1(inoutP,newf,f0,dir,Gra,priorgra,x,GradSubroutine,maxflag,point_parameter)
		if(stop_type.eq.1)then
			stop_value=dnorm2(Gra)
		else if(stop_type.eq.2)then
			stop_value=dabs((newf-f0)/(newf+f0))
		else
			write(*,*)"ERROR stop type"
		end if
		if(print_flag)then
			if(stop_type.eq.1)then
				write(*,*)"newf,|gradient|^2,len of step(x),direction_corr"
			else if(stop_type.eq.2)then
				write(*,*)"newf,dabs((newf-f0)/(newf+f0)),len of step(x),direction_corr"
			end if
		end if
		runtyp6flag=run_type.eq.6
		runtyp6counter=0
		typecounter=1
		if(runtyp6flag) call set_run_type(which_type(typecounter),x) 
		
		do while(stop_value.gt.stop_error)
			f0=newf
			i=i+1
			call go_one_step1(inoutP,newf,f0,dir,Gra,priorgra,x,GradSubroutine,maxflag,point_parameter)
			if(stop_type.eq.1)then
				stop_value=dnorm2(Gra)
			else if(stop_type.eq.2)then
				stop_value=dabs((newf-f0)/(newf+f0))
			end if
			if( (run_type.eq.3) .or. (run_type.eq.4 ) .or. (run_type.eq.5 )) then
				if(mod(i,num_same_step).eq.0) then
					x=x*delta_step
				end if
			end if
			if( (run_type.eq.4 ) .or. (run_type.eq.5 )) then
				if(x.lt.min_step_when_num_same_step)then
					max_step=min_step_when_num_same_step
					if(run_type.eq.4) then
						call set_run_type(2,x)
					else
						call set_run_type(1,x)
					end if
				end if
			end if
			if( (run_type.eq.3 ) ) then
				if(x.lt.min_step_when_num_same_step)then
					if(stop_type.eq.1) then
						write(*,*)"wornning,|gradient|^2",stop_value
					else if(stop_type.eq.2) then
						write(*,*)"wornning,dabs((newf-f0)/(newf+f0))",stop_value
					end if
					exit
				end if
			end if
			if(runtyp6flag)then
				if((run_type.eq.1) .or. (run_type.eq.2 ) .or. (run_type.eq.3 ))then
					runtyp6counter=runtyp6counter+1
					if(runtyp6counter.eq.step_run_type(typecounter))then
						typecounter=typecounter+1
						call set_run_type(which_type(typecounter),x) 
					end if
				end if
			end if
			if(print_flag)then
				write(*,*)newf,stop_value,x,direction_corr
			end if
			if(i.gt.max_running) then
				if(stop_type.eq.1)then
					write(*,*)"wornning,|gradient|^2",stop_value
				else if(stop_type.eq.2)then
					write(*,*)"wornning,dabs((newf-f0)/(newf+f0))",stop_value
				end if
				exit
			end if
		end do
		CG_Method1_inputP=newf
		if(print_flag)then
			write(*,*)"num_of_step_to_finish_search",i
		end if
		return
	end function


	real*8 function CG_Method2(inoutP,lenofinP,searchFunction,maxflag,point_parameter)
		type(DTensor),intent(inout)::inoutP
		logical,intent(in)::maxflag
		integer,intent(in)::lenofinP
		real*8,external::searchFunction
		type(DTensor),optional,intent(inout)::point_parameter
		type(DTensor)::dir,Gra,priorgra
		real*8::newf,f0,x,stop_value
		integer::i,lenofP,runtyp6counter,typecounter
		logical::runtyp6flag
		if(.not.DgetFlag(inoutP)) then
				lenofP=lenofinP
				inoutP=Dgenerate((/lenofP/),(/-1d0,1d0/))
		else
				write(*,*)"ERROR in CG_Method"
				write(*,*)"one shoould input a empty Tensor"
				stop
		end if
		i=0
		x=0d0
		call emptyDTensor(dir)
		call emptyDTensor(Gra)
		call go_one_step2(inoutP,newf,f0,dir,Gra,priorgra,x,searchFunction,maxflag,point_parameter)
		if(stop_type.eq.1)then
			stop_value=dnorm2(Gra)
		else if(stop_type.eq.2)then
			stop_value=dabs((newf-f0)/(newf+f0))
		else
			write(*,*)"ERROR stop type"
		end if
		if(print_flag)then
			if(stop_type.eq.1)then
				write(*,*)"newf,|gradient|^2,len of step(x),direction_corr"
			else if(stop_type.eq.2)then
				write(*,*)"newf,dabs((newf-f0)/(newf+f0)),len of step(x),direction_corr"
			end if
		end if
		runtyp6flag=run_type.eq.6
		runtyp6counter=0
		typecounter=1
		if(runtyp6flag) call set_run_type(which_type(typecounter),x) 
		
		do while(stop_value.gt.stop_error)
			i=i+1
			f0=newf
			call go_one_step2(inoutP,newf,f0,dir,Gra,priorgra,x,searchFunction,maxflag,point_parameter)
			if(stop_type.eq.1)then
				stop_value=dnorm2(Gra)
			else if(stop_type.eq.2)then
				stop_value=dabs((newf-f0)/(newf+f0))
			end if
			if( (run_type.eq.3) .or. (run_type.eq.4 ) .or. (run_type.eq.5 )) then
				if(mod(i,num_same_step).eq.0) then
					x=x*delta_step
				end if
			end if
			if( (run_type.eq.4 ) .or. (run_type.eq.5 )) then
				if(x.lt.min_step_when_num_same_step)then
					max_step=min_step_when_num_same_step
					if(run_type.eq.4) then
						call set_run_type(2,x)
					else
						call set_run_type(1,x)
					end if
				end if
			end if
			if( (run_type.eq.3 ) ) then
				if(x.lt.min_step_when_num_same_step)then
					if(stop_type.eq.1) then
						write(*,*)"wornning,|gradient|^2",stop_value
					else if(stop_type.eq.2) then
						write(*,*)"wornning,dabs((newf-f0)/(newf+f0))",stop_value
					end if
					exit
				end if
			end if
			if(runtyp6flag)then
				if((run_type.eq.1) .or. (run_type.eq.2 ) .or. (run_type.eq.3 ))then
					runtyp6counter=runtyp6counter+1
					if(runtyp6counter.eq.step_run_type(typecounter))then
						typecounter=typecounter+1
						call set_run_type(which_type(typecounter),x) 
					end if
				end if
			end if
			if(print_flag)then
				write(*,*)newf,stop_value,x,direction_corr
			end if
			if(i.gt.max_running) then
				if(stop_type.eq.1)then
					write(*,*)"wornning,|gradient|^2",stop_value
				else if(stop_type.eq.2)then
					write(*,*)"wornning,dabs((newf-f0)/(newf+f0))",stop_value
				end if
				exit
			end if
		end do
		CG_Method2=newf
		if(print_flag)then
			write(*,*)"num_of_step_to_finish_search",i
		end if
		return
	end function
	real*8 function CG_Method2_inputP(inoutP,searchFunction,maxflag,point_parameter)
		type(DTensor),intent(inout)::inoutP
		logical,intent(in)::maxflag
		real*8,external::searchFunction
		type(DTensor),optional,intent(inout)::point_parameter
		type(DTensor)::dir,Gra,priorgra
		real*8::newf,f0,x,stop_value
		integer::i,lenofP,runtyp6counter,typecounter
		logical::runtyp6flag
		if(.not.DgetFlag(inoutP)) then
				write(*,*)"ERROR in CG_Method"
				write(*,*)"Specify the length of input point,or input a initial point"
				stop
		else
			lenofP=DgettotalData(inoutP)
		end if
		i=0
		x=0d0
		call emptyDTensor(dir)
		call emptyDTensor(Gra)
		call go_one_step2(inoutP,newf,f0,dir,Gra,priorgra,x,searchFunction,maxflag,point_parameter)
		if(stop_type.eq.1)then
			stop_value=dnorm2(Gra)
		else if(stop_type.eq.2)then
			stop_value=dabs((newf-f0)/(newf+f0))
		else
			write(*,*)"ERROR stop type"
		end if
		if(print_flag)then
			if(stop_type.eq.1)then
				write(*,*)"newf,|gradient|^2,len of step(x),direction_corr"
			else if(stop_type.eq.2)then
				write(*,*)"newf,dabs((newf-f0)/(newf+f0)),len of step(x),direction_corr"
			end if
		end if
		runtyp6flag=run_type.eq.6
		runtyp6counter=0
		typecounter=1
		if(runtyp6flag) call set_run_type(which_type(typecounter),x) 
		do while(stop_value.gt.stop_error)
			i=i+1
			f0=newf
			call go_one_step2(inoutP,newf,f0,dir,Gra,priorgra,x,searchFunction,maxflag,point_parameter)
			if(stop_type.eq.1)then
				stop_value=dnorm2(Gra)
			else if(stop_type.eq.2)then
				stop_value=dabs((newf-f0)/(newf+f0))
			end if
			if( (run_type.eq.3) .or. (run_type.eq.4 ) .or. (run_type.eq.5 )) then
				if(mod(i,num_same_step).eq.0) then
					x=x*delta_step
				end if
			end if
			if( (run_type.eq.4 ) .or. (run_type.eq.5 )) then
				if(x.lt.min_step_when_num_same_step)then
					max_step=min_step_when_num_same_step
					if(run_type.eq.4) then
						call set_run_type(2,x)
					else
						call set_run_type(1,x)
					end if
				end if
			end if
			if( (run_type.eq.3 ) ) then
				if(x.lt.min_step_when_num_same_step)then
					if(stop_type.eq.1) then
						write(*,*)"wornning,|gradient|^2",stop_value
					else if(stop_type.eq.2) then
						write(*,*)"wornning,dabs((newf-f0)/(newf+f0))",stop_value
					end if
					exit
				end if
			end if
			if(runtyp6flag)then
				if((run_type.eq.1) .or. (run_type.eq.2 ) .or. (run_type.eq.3 ))then
					runtyp6counter=runtyp6counter+1
					if(runtyp6counter.eq.step_run_type(typecounter))then
						typecounter=typecounter+1
						call set_run_type(which_type(typecounter),x) 
					end if
				end if
			end if
			if(print_flag)then
					write(*,*)newf,stop_value,x,direction_corr
			end if
			if(i.gt.max_running) then
				if(stop_type.eq.1) then
					write(*,*)"wornning,|gradient|^2",stop_value
				else if(stop_type.eq.2) then
					write(*,*)"wornning,dabs((newf-f0)/(newf+f0))",stop_value
				end if
				exit
			end if
		end do
		CG_Method2_inputP=newf
		if(print_flag)then
			write(*,*)"num_of_step_to_finish_search",i
		end if
		return
	end function














end module
