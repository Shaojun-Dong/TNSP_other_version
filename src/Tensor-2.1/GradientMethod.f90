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
!	type(DTensor), and output a real*8 value.
!
!		2. When calling Gradient_Method without initialation,
!	the program will read the parameter in ./GMparameter.
!	So one should create this file. Or 
!			1.call intGradientMethod(GMparameter_address)
!			2.call intGradientMethod(diff_delta_,max_step_,
!					delta_step_,stop_error_,max_running_,p_flag_)
!	to initialation the program.
!	
!		3. The data in GMparameter are:
!----------------------------------------------------------------------------------------------
!diff_delta			1d-8		differentiation,(f(x+delta)-f(x))/delta,diff_delta is the delta
!max_step			5d0		go along the dirction with max_step,new_x=x+max_step
!delta_step			0.9		if f(new_x)<f(x),max_step=max_step*delta_step,delta_step should smaller than 1 
!stop_error			1d-8		if abs(f(new_x)-f(x))<stop_error,output
!max_running		5000		max running
!print_flag			1			if 1,print the progress
!----------------------------------------------------------------------------------------------
!
!		4. Send email to tyrants@qq.com to report any bugs.
!
!****************************************************
!****************************************************

include "/home/sjdong/Tensor/MPI/Tensor.f90"
!include "/home/sjdong/Tensor/MPI/constant.f90"
module GradientMethod
!	use	constant
	use Tensor_complex
	implicit none
!**********************************************
!*****   GradientMethod parameter   ***********
	logical,private::flag=.false.,print_flag=.false.
	real*8,private::diff_delta,max_step,delta_step,stop_error
	integer,private::max_running
	private::diff,gradient,Sdirection,Pstep
	interface Gradient_Method
		module procedure Gradient_Method1
		module procedure Gradient_Method2
		module procedure Gradient_Method_para
		module procedure Gradient_Method2_para
	end interface
!		real*8 a=Gradient_Method(inoutP,max_flag,lenofinP,func) 
!		real*8 a=Gradient_Method(inoutP,para,max_flag,lenofinP,func) 
!			inoutP is the initial point,type(DTensor) for searching the max/min
!		if inoutP is no data, one should set lenofinP as the
!		length of inoutP, the program will create random data

!			lenofinP is optional. It is the length of inoutP,if this are data in inoutP
!
!			max_flag=.true.,find the max,if not ,min
!			
!			The search function is searchFunction(Ten),or searchFunction(Ten,para),Search inoutP only
!
!
	interface intGradientMethod
		module procedure intGradientMethod1
		module procedure intGradientMethod2
		module procedure intGradientMethod3
	end interface
!		intGradientMethod(diff_delta_,max_step_,delta_step_,stop_error_,max_running_,p_flag_)
!	or	intGradientMethod(GMparameter_address),the input file can copy from Tensor/MPI/GMparameter
!	or intGradientMethod() read the file in ./GMparameter
!**********************************************
!**********************************************


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
!	real*8 function searchFunction(point,para)
!		type(DTensor),intent(in)::point,para
!
!.......  write the function here .......
!		
!		return
!	end function
!**********************************************
!**********************************************


!example program
!	real*8 maxpoint
!	type(DTensor)::point
!	maxpoint=Gradient_Method(point,.true.,3,Example_searchFunction) 
!or
!	!	maxpoint=Gradient_Method(point,para,.true.,3,Example_searchFunction) 
!---------------------------------------------------------

!**********************************************
contains
	real*8 function Example_searchFunction(Ten)!this is the function to be minimize or maximize
		type(DTensor),intent(in)::Ten
		real*8::x,y,z
		x=Ten.i.1
		y=Ten.i.2
		z=Ten.i.3
		Example_searchFunction=-(x-1)**2-(y-2)**2+1
		return
	end function
	
	subroutine intGradientMethod1()
		CHARACTER*100::notused
		integer::p_flag
		open(unit=10,file="GMparameter",STATUS='old')
		read(10,*)notused,diff_delta!differentiation,(f(x+delta)-f(x))/delta,diff_delta is the delta
		read(10,*)notused,max_step!go along the dirction with max_step,new_x=x+max_step
		read(10,*)notused,delta_step!if f(new_x)<f(x),max_step=max_step*delta_step
		read(10,*)notused,stop_error!if abs(f(new_x)-f(x))<stop_error,output
		read(10,*)notused,max_running!max running
		read(10,*)notused,p_flag!if 1,print the progress
		if(p_flag.eq.1)then
			print_flag=.true.
		end if
		flag=.true.!if flag，means one have run this intGradientMethod
		close(unit=10)
		return
	end subroutine
	subroutine intGradientMethod2(GMparameter_address)
		CHARACTER(len=*),intent(in)::GMparameter_address
		CHARACTER*100::notused
		integer::p_flag
		open(unit=10,file=GMparameter_address,STATUS='old')
		read(10,*)notused,diff_delta!differentiation,(f(x+delta)-f(x))/delta,diff_delta is the delta
		read(10,*)notused,max_step!go along the dirction with max_step,new_x=x+max_step
		read(10,*)notused,delta_step!if f(new_x)<f(x),max_step=max_step*delta_step
		read(10,*)notused,stop_error!if abs(f(new_x)-f(x))<stop_error,output
		read(10,*)notused,max_running!max running
		read(10,*)notused,p_flag!if 1,print the progress
		if(p_flag.eq.1)then
			print_flag=.true.
		end if
		flag=.true.!if flag，means one have run this intGradientMethod
		close(unit=10)
		return
	end subroutine
	
	subroutine intGradientMethod3(diff_delta_,max_step_,delta_step_,stop_error_,max_running_,p_flag_)
		logical,intent(in)::p_flag_
		real*8,intent(in)::diff_delta_,max_step_,delta_step_,stop_error_
		integer,intent(in)::max_running_
		diff_delta=diff_delta_!differentiation,(f(x+delta)-f(x))/delta,diff_delta is the delta
		max_step=max_step_!go along the dirction with max_step,new_x=x+max_step
		delta_step=delta_step_!if f(new_x)<f(x),max_step=max_step*delta_step
		stop_error=stop_error_!if abs(f(new_x)-f(x))<stop_error,output
		max_running=max_running_!max running
		print_flag=p_flag_!if 1,print the progress
		flag=.true.!if flag，means one have run this intGradientMethod
		return
	end subroutine
	
	
	

!**********************************************************
!**********************************************************
!**********************************************************
!**********No not modify the function below****************
!**********************************************************
!**********************************************************
!**********************************************************
	
!differentiation of searchFunction
	real*8 function diff(inputP,ith,GMsearchFunction)
		type(DTensor),intent(in)::inputP
		integer,intent(in)::ith
		real*8,external::GMsearchFunction
		type(Dtensor)::delta_P,P2,P1
		real*8,allocatable::Pdata(:)
		integer::lenP
		lenP=DgetTotaldata(inputP)
		allocate(Pdata(lenP))
		Pdata=0d0
		Pdata(ith)=diff_delta
		delta_P=Pdata
		P2=inputP+delta_P
		P1=inputP-delta_P
		diff=GMsearchFunction(P2)-GMsearchFunction(P1)
		diff=diff/(diff_delta+diff_delta)
		if(isnan(diff)) then
			write(*,*)"NAN error,diff"
			stop
		end if
		return
	end function
!the search function is GMsearchFunction(inoutP,point_parameter)
!Search inoutP
	real*8 function diff_para(inputP,point_parameter,ith,GMsearchFunction) result(diff)
		type(DTensor),intent(in)::inputP
		type(DTensor),intent(in)::point_parameter
		integer,intent(in)::ith
		real*8,external::GMsearchFunction
		type(Dtensor)::delta_P,P2,P1
		real*8,allocatable::Pdata(:)
		integer::lenP
		lenP=DgetTotaldata(inputP)
		allocate(Pdata(lenP))
		Pdata=0d0
		Pdata(ith)=diff_delta
		delta_P=Pdata
		P2=inputP+delta_P
		P1=inputP-delta_P
		diff=GMsearchFunction(P2,point_parameter)-GMsearchFunction(P1,point_parameter)
		diff=diff/(diff_delta+diff_delta)
		if(isnan(diff)) then
			write(*,*)"NAN error,diff"
			stop
		end if
		return
	end function
! gradient vector of searchFunction
	type(DTensor) function gradient(P,GMsearchFunction)
		type(DTensor),intent(in)::P
		real*8,allocatable::Pdata(:)
		real*8,external::GMsearchFunction
		integer::lenP,i
		lenP=DgetTotaldata(P)
		allocate(Pdata(lenP))
		do i=1,lenP
			Pdata(i)=diff(P,i,GMsearchFunction)
		end do
		gradient=Pdata
		return
	end function
!the search function is GMsearchFunction(inoutP,point_parameter)
!Search inoutP
	type(DTensor) function gradient_para(P,point_parameter,GMsearchFunction)result(gradient)
		type(DTensor),intent(in)::P,point_parameter
		real*8,allocatable::Pdata(:)
		real*8,external::GMsearchFunction
		integer::lenP,i
		lenP=DgetTotaldata(P)
		allocate(Pdata(lenP))
		do i=1,lenP
			Pdata(i)=diff_para(P,point_parameter,i,GMsearchFunction)
		end do
		gradient=Pdata
		return
	end function
! Gradient Method,search direction  
	type(DTensor) function Sdirection(P,max_flag,GMsearchFunction)
		type(DTensor)::S
		type(DTensor),intent(in)::P
		logical,intent(in)::max_flag
		real*8,external::GMsearchFunction
		real*8::Snorm
		S=gradient(P,GMsearchFunction)
		Snorm=dnorm(S)
		if(Snorm.eq.0)then
			Sdirection=S
			write(*,*)"gradient is 0"
			call DTMprint(S,1)
			call DTMprint(P,1)
			return
		end if
		if(max_flag)then
			Sdirection=S/Snorm
		else
			Sdirection=S/Snorm
			Sdirection=(-1d0)*Sdirection
		end if
		return
	end function
!the search function is GMsearchFunction(inoutP,point_parameter)
!Search inoutP
	type(DTensor) function Sdirection_para(P,point_parameter,max_flag,GMsearchFunction)result(Sdirection)
		type(DTensor)::S
		type(DTensor),intent(in)::point_parameter
		type(DTensor),intent(in)::P
		logical,intent(in)::max_flag
		real*8,external::GMsearchFunction
		real*8::Snorm
		S=gradient_para(P,point_parameter,GMsearchFunction)
		Snorm=dnorm(S)
		if(Snorm.eq.0)then
			Sdirection=S
			write(*,*)"gradient is 0"
			call DTMprint(S,1)
			call DTMprint(P,1)
			return
		end if
		if(max_flag)then
			Sdirection=S/Snorm
		else
			Sdirection=S/Snorm
			Sdirection=(-1d0)*Sdirection
		end if
		return
	end function
	real*8 function Pstep(P,max_flag,GMsearchFunction)
		type(DTensor),intent(inout)::P
		logical,intent(in)::max_flag
		real*8,external::GMsearchFunction
		real*8::f0,f1,step
		real*8::Tendata(3)
		type(DTensor)::Sdir,P1,P0
		P0=P
		step=max_step
		Sdir=Sdirection(P0,max_flag,GMsearchFunction)
		f0=GMsearchFunction(P0)
		P=P0+(step*Sdir)
		f1=GMsearchFunction(P)
		if(max_flag)then
			do while(f0.gt.f1) 
				step=step*delta_step
				P=P0+(step*Sdir)
				f1=GMsearchFunction(P)
			end do
		else
			do while(f0.lt.f1) 
				step=step*delta_step
				P=P0+(step*Sdir)
				f1=GMsearchFunction(P)
			end do
		end if
		Pstep=f1
		return
	end function
!the search function is GMsearchFunction(inoutP,point_parameter)
!Search inoutP
	real*8 function Pstep_para(P,point_parameter,max_flag,GMsearchFunction)result(Pstep)
		type(DTensor),intent(inout)::P
		type(DTensor),intent(in)::point_parameter
		logical,intent(in)::max_flag
		real*8,external::GMsearchFunction
		real*8::f0,f1,step
		real*8::Tendata(3)
		type(DTensor)::Sdir,P1,P0
		P0=P
		step=max_step
		Sdir=Sdirection_para(P0,point_parameter,max_flag,GMsearchFunction)
		f0=GMsearchFunction(P0,point_parameter)
		P=P0+(step*Sdir)
		f1=GMsearchFunction(P,point_parameter)
		if(max_flag)then
			do while(f0.gt.f1) 
				step=step*delta_step
				P=P0+(step*Sdir)
				f1=GMsearchFunction(P,point_parameter)
			end do
		else
			do while(f0.lt.f1) 
				step=step*delta_step
				P=P0+(step*Sdir)
				f1=GMsearchFunction(P,point_parameter)
			end do
		end if
		Pstep=f1
		return
	end function
! 
	real*8 function Gradient_Method1(inoutP,max_flag,lenofinP,GMsearchFunction)
		type(DTensor),intent(inout)::inoutP
		logical,intent(in)::max_flag
		integer,intent(in)::lenofinP
		real*8,external::GMsearchFunction
		type(DTensor)::P,res
		real*8::newf,f0,tempnum
		integer::i,templen,lenofP
		if(.not.flag) then
			call intGradientMethod()
		end if
		if(.not.DgetFlag(inoutP)) then
				lenofP=lenofinP
				P=Dgenerate((/lenofP/))
		else
				write(*,*)"ERROR in Gradient_Method"
				write(*,*)"one shoould input a empty Tensor"
				stop
		end if
		f0=GMsearchFunction(P)
		newf=Pstep(P,max_flag,GMsearchFunction)
		i=1
		do while(dabs(newf-f0).gt.stop_error)
			f0=newf
			newf=Pstep(P,max_flag,GMsearchFunction)
			i=i+1
			if(i.gt.max_running) then
				write(*,*)"wornning,delta f/f",(newf-f0)/(newf+f0)
				exit
			end if
			if(print_flag)then
				write(*,*)newf,(newf-f0)/(newf+f0)
			end if
		end do
		inoutP=P
		Gradient_Method1=newf
		return
	end function
	real*8 function Gradient_Method2(inoutP,max_flag,GMsearchFunction)
		type(DTensor),intent(inout)::inoutP
		logical,intent(in)::max_flag
		real*8,external::GMsearchFunction
		type(DTensor)::P,res
		real*8::newf,f0,tempnum
		integer::i,templen,lenofP
		if(.not.flag) then
			call intGradientMethod()
		end if
		if(.not.DgetFlag(inoutP)) then
				write(*,*)"ERROR in Gradient_Method"
				write(*,*)"Specify the length of input point,or input a initial point"
				stop
		else
			P=inoutP
			lenofP=DgettotalData(inoutP)
		end if
		f0=GMsearchFunction(P)
		newf=Pstep(P,max_flag,GMsearchFunction)
		i=1
		do while(dabs(newf-f0).gt.stop_error)
			f0=newf
			newf=Pstep(P,max_flag,GMsearchFunction)
			i=i+1
			if(i.gt.max_running) then
				write(*,*)"wornning,delta f/f",(newf-f0)/(newf+f0)
				exit
			end if
			if(print_flag)then
				write(*,*)newf,(newf-f0)/(newf+f0)
			end if
		end do
		inoutP=P
		Gradient_Method2=newf
		return
	end function
!the search function is GMsearchFunction(inoutP,point_parameter)
!Search inoutP
	real*8 function Gradient_Method_para(inoutP,point_parameter,max_flag,lenofinP,GMsearchFunction)
		type(DTensor),intent(inout)::inoutP
		type(DTensor),intent(in)::point_parameter
		logical,intent(in)::max_flag
		integer,intent(in)::lenofinP
		real*8,external::GMsearchFunction
		type(DTensor)::P,res
		real*8::newf,f0,tempnum
		integer::i,templen,lenofP
		if(.not.flag) then
			call intGradientMethod()
		end if
		if(.not.DgetFlag(inoutP)) then
				lenofP=lenofinP
				P=Dgenerate((/lenofP/))
		else
				write(*,*)"ERROR in Gradient_Method"
				write(*,*)"one shoould input a empty Tensor"
				stop
		end if
		f0=GMsearchFunction(P,point_parameter)
		newf=Pstep_para(P,point_parameter,max_flag,GMsearchFunction)
		i=1
		do while(dabs(newf-f0).gt.stop_error)
			f0=newf
			newf=Pstep_para(P,point_parameter,max_flag,GMsearchFunction)
			i=i+1
			if(i.gt.max_running) then
				write(*,*)"wornning,delta f/f",(newf-f0)/(newf+f0)
				exit
			end if
			if(print_flag)then
				write(*,*)newf,(newf-f0)/(newf+f0)
			end if
		end do
		inoutP=P
		Gradient_Method_para=newf
		return
	end function
	real*8 function Gradient_Method2_para(inoutP,point_parameter,max_flag,GMsearchFunction)
		type(DTensor),intent(inout)::inoutP
		type(DTensor),intent(in)::point_parameter
		logical,intent(in)::max_flag
		real*8,external::GMsearchFunction
		type(DTensor)::P,res
		real*8::newf,f0,tempnum
		integer::i,templen,lenofP
		if(.not.flag) then
			call intGradientMethod()
		end if
		if(.not.DgetFlag(inoutP)) then
				write(*,*)"ERROR in Gradient_Method"
				write(*,*)"Specify the length of input point,or input a initial point"
				stop
		else
			P=inoutP
			lenofP=DgettotalData(inoutP)
		end if
		f0=GMsearchFunction(P,point_parameter)
		newf=Pstep_para(P,point_parameter,max_flag,GMsearchFunction)
		i=1
		do while(dabs(newf-f0).gt.stop_error)
			f0=newf
			newf=Pstep_para(P,point_parameter,max_flag,GMsearchFunction)
			i=i+1
			if(i.gt.max_running) then
				write(*,*)"wornning,delta f/f",(newf-f0)/(newf+f0)
				exit
			end if
			if(print_flag)then
				write(*,*)newf,(newf-f0)/(newf+f0)
			end if
		end do
		inoutP=P
		Gradient_Method2_para=newf
		return
	end function
end module	




















