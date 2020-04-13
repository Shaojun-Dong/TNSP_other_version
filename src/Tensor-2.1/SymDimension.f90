!************************************************************
!************* START OF SymDimension *******************
!************************************************************
!*****************     worning    *********************************
!				1.In fortran,the data store in memory as:A_{1,1}->A_{2,1}->A_{3,1}
!			->...A{M,1}->A_{1,2}->A_{2,2}->..->A_{M,2}->..->A_{M-1,N}->A_{M,N}.
!			But in matlab or C,it is as:A{0,0}->A_{0,1}->A_{0,2}->..A{0,N}
!			->A{1,0}->..->A{M,N-1}->A_{M,N}.The code of  SymDimConstract and  SymDimDecompose
!			is base on the rule of storing data.
!				2.A is a allocatable array and B is array too,when use A=B before
!			allocating storage space for A,A will get no data.This happen in 
!			fortran90 but do not in fortran77.The interface assignments is to  
!			solve this program.
!				3.The code is fortran90 version.
!******************   note   ********************
!				DimData :store the SymDimension data
!				boundary: store the boundary of every SymDimension,for
!			 example,the SymDimension D1=[2,2,3,4,5],after constract the
!			 first and second index ,it will be D2=[4,3,4,5],then
!			 the data of D1 :
!						boundarysize	=6
!						Dimsize			=5
!						boundary		=[0,1,2,3,4,5]
!						DimData			=[2,2,3,4,5]
!					which means the SymDimension is [2,2,3,4,5]
!			 the data of D2 :
!						boundarysize	=5
!						Dimsize			=4
!						boundary		=[0,2,3,4,5]
!						DimData			=[2,2,3,4,5]
!				which means the SymDimension is [4,3,4,5] (or[(2*2),3,2,4] )
!				the DimData will not change,but will the boundary,the boundary(i)+1
!			 and boundary(i+1) is the i dimenison,for the case D2,boundary(1)+1
!			 =0+1=1,boundary(1+1)=2,then the DimData(1) and DimData(2) is the first
!			 dimenison,such dimenison is 2*2=4
!				Dimsize is not the length of DimData. But the rank of The Tensor.
!			 	such of storing data is convenient for constraction and decomposion.
!				Send email to tyrants@qq.com to report any bugs.
!			This is the code for MPI.It is the same as the non-MPI case exacpt the
!		MPI code at the end of this file.
!***************************************************************************	
!		2014.11.13
!			Add a new type SymDimensionName into type(SymDimension)  
!						DimData			=[2,2,3,4,5]
!						SymDimName		=[1,2,3,4,5]  
!                              [a,a,a,a,a] beause this is the same Tensor,indexname are the same
!    SymDimName are
!		[1,a],[2,a],[3,a],[4,a],[5,a]
!					which means the SymDimension is [2,2,3,4,5],the first index call "a.1" second call "a.2",
!		When permute if to 
!												[2,3,2,4,5]
!     The value of indexID will be
!												[1,3,2,4,5]
!		 SymDimName are
!		[1,a],[3,a],[2,a],[4,a],[5,a]
!		The name will change according to the permutation, but the value will not change. So is nameID.
!		nameID is the same in a SymDimension. When there are two dimesion. 
!		It is easy to get the order of the index whose name is '3'.
!		If there are two Tensor [2,2] \times [2,3,2,4]
!		the first Tensor :
!						[1,a],[2,a]
!		the second Tensor :
!						[1,b],[2,b],[3,b],[4,b]
!		the contract the i2 of the first and i3 of the second,result Tensor is 
!			[1,a],[1,b],[2,b],[4,b]
!
!			The best way to achieve the function of Name is to redefine DimData to store the name.
!		That	is
!	type SymDimensionData
!		CHARACTER(len=len_of_Name)::TensorName
!		CHARACTER(len=len_of_Name)::SymDimName
!	end type SymDimensionData
!
!	type SymDimension
!		integer,private :: boundarysize=0
!		integer,private :: Dimsize=0
!		integer,allocatable,private :: boundary(:)
!		type(SymDimensionData),allocatable,private::dimData(:)
!		logical,private::nameflag=.false.
!	end type SymDimension
!			But it cost all lot to modify in the files.
!***************
!
! original SymDimension means that the SymDimension do not do any contract 
!

!index of the Tensor is 1,2,3...N
!index of the physis Tensor should be related to the quantum number(QN)
!If the index of the Tensor is(QN=0.5)
!index of Tensor		[1]
!index of SymTen		[0]
!
!1/2 spin:
!index of Tensor		 [1,2]  
!index of SymTen		[-0.5,0.5]
!
!1 spin:
!index of Tensor		 [1,2,3]
!index of SymTen		[-1,0,1]

!3/2 spin:
!index of Tensor		 [1,2,3,4]
!index of SymTen		[-3/2,-1/2,1/2,3/2]
!
!index of Tensor		 [1,2,...,N]
!index of SymTen		[-(N-1)*Qnum_,-(N-1)*Qnum_+Delta_Q_,..,(N-1)*Qnum_-Delta_Q_,(N-1)*Qnum_]
!the rule for partition,in U(1),the rule is \sum S_in=\sum S_out(prb 83,115125;pra 82,050301)
!==>\sum S_in - \sum S_out =0
!==>rule(1)*(N_1-1)*Qnum_+rule(1)*(N_2-1)*Qnum_+..+rule(n)*(N_n-1)*Qnum_=0
!N_i is the ith dimension
!
!   s1 -->--[Tensor]--->-- s3    : a [2*3*4] Tensor,the non-zero partitioned Tensor is 1*s1+1*s2-1*s3=0
!              |                rule:[1,1,-1] ,1 mean the QN gose in, and -1 mean gose out
!              ^                In svd ,use for determine the new dimenion of QTensor
!              | s2
!
!example of degeneracy:
!QN:1,0.5,1.5
!   1 -->--[Tensor]--->-- 1.5 
!              |        
!              ^ 
!              | 0.5  
!first index of SymTen[-1,0,1],degeneracy=[1,2,1]( [d1_1,d1_2,d1_3])
!second index of SymTen[-0.5,0.5],degeneracy=[1,1]( [d2_1,d2_2])
!third index of SymTen[-1.5,-0.5,0.5,1.5],degeneracy=[1,2,2,1]( [d3_1,d3_2,d3_3,d3_4])
!The non-zero partitioned Tensors of T_{s1,s2,s3} are
!T_{-1,-0.5,-1.5} a Tesnor of 1*1*1 ,or d1_1*d2_1*d3_1
!T_{-1, 0.5,-0.5} a Tesnor of 1*1*2 ,or d1_1*d2_2*d3_2
!T_{ 0,-0.5,-0.5} a Tesnor of 2*1*2 ,or d1_2*d2_1*d3_2
!T_{ 0, 0.5, 0.5} a Tesnor of 2*1*2 ,or d1_2*d2_2*d3_3
!T_{ 1,-0.5, 0.5} a Tesnor of 1*1*2 ,or d1_3*d2_1*d3_3
!T_{ 1, 0.5, 1.5} a Tesnor of 1*1*1 ,or d1_3*d2_2*d3_3
module SymDimension_typede
	use eigen_value
	implicit none
	CHARACTER*1,save,private::indexsymbol='.'
	!integer,parameter::len_of_Name=20
	integer,parameter::Delta_QN_=1
	real*8,parameter::zero_error_=1d-15
	type SymDimensionName
		CHARACTER(len=len_of_Name)::TensorName
		CHARACTER(len=len_of_Name)::dimenName
		real*4,allocatable::QN(:)!quantum number
		integer,allocatable::degeneracy(:)
		real*4,private::max_QN
		integer::rule! rule=1 or -1  in U(1) Sin=Sout==> Sin-Sout=0==>Sin*rule(Sin)+Sout*rule(Sout)=0
	end type SymDimensionName
	type SymDimension
		integer,private :: boundarysize=0
		integer,private :: Dimsize=0
		integer,allocatable,private :: boundary(:)
		integer,allocatable,private :: DimData(:)
		type(SymDimensionName),allocatable,private::DimName(:)
		logical,private::nameflag=.false.
	end type SymDimension
	
	interface SymNameinit
		module procedure SymNameinit1
		module procedure SymNameinit2
	end interface	
	interface setDimQN
		module procedure setDimQN1
		module procedure setDimQN2
		module procedure setDimQN3
	end interface
	interface reverseDimrule
		module procedure reverseDimrule1
		module procedure reverseDimrule2
	end interface	
	!find the index,whose name is w	,output the order of it in the SymDimension
	!If can not find , output 0
	interface SymNameorder
		module procedure SymNameorder2!(dimen,character)
		module procedure SymNameorder3!(dimen,name)
		module procedure SymNameorderArray!(dimen,character(:)),output a integer(:)
	end interface
	
	interface SymDimDegeneracy
		module procedure SymDimDegeneracy1
		module procedure SymDimDegeneracy2
		module procedure SymDimDegeneracy3
		module procedure SymDimDegeneracy4
	end interface
	interface symResetdeg
		module procedure symResetdeg1
		module procedure symResetdeg2
		module procedure symResetdeg3
	end interface
	interface SymQN
		module procedure SymQN1
		module procedure SymQN2
	end interface
	interface SymQNum
		module procedure SymQNum1
		module procedure SymQNum2
	end interface
	interface SymRule
		module procedure SymRule1
		module procedure SymRule2
	end interface
	interface SymU1Dim
		module procedure SymU1Dim1
		module procedure SymU1Dim2
		module procedure SymU1Dim3
	end interface
	
!Find all the SymDimensionNames, whose indexname is oldname
!change them to newname
!If Cannot find it , do nothing
!If input 'A' 'B', all name such as 'A.1','A.2','A.as' will change to 'B.1','B.2','B.as'
!if input 'A.1' 'B.1',the name 'A.1' will change to 'B.1'
	interface SymResetname
		module procedure SymResetname1
		module procedure SymResetname2
	end interface
	interface assignment(=)
		module procedure SymDimInitialization
		module procedure SymDimInitialization2
		!module procedure SymgetDim
		module procedure SymDimToint
		module procedure SymDimName!(dimen,w),set a name to the dimenison
		module procedure SymcopyName!name1=name2
		module procedure SymcopyNameArray!name1(:)=name2(:),willnot allocate name1
		module procedure SymcharName!name=character
		module procedure SymNamechar!character=name
		module procedure SymcharNameArray!name(:)=character(:)
		module procedure SymNamecharArray!character(:)=name(:)
	end interface
	interface operator(+)
		module procedure SymDimadd
		module procedure SymDimadd2
		module procedure SymDimadd3
	end interface
	interface operator(.sub.)
		module procedure SymgetSubDim2
		module procedure SymgetSubDim2_name
	end interface
	interface operator(.i.)
		module procedure SymDim_i
	!	module procedure SymNameorderAllocateArray
		!the same function with SymNameorderArray,but input output a allocatable vector
	end interface
	interface operator(.equ.)
		module procedure  Symequal_of_dim
		module procedure  Symequal_name1
		module procedure  Symequal_name2
		module procedure  Symequal_name3
	end interface
	
!**********************************************************

!
!**********************************************************	
	contains
	logical function Symif_original_dim(dimen)
		type(SymDimension),intent(in)::dimen
		if(size(dimen%DimData).eq.dimen%Dimsize) then	
			Symif_original_dim=.true.
		else
			Symif_original_dim=.false.
		end if
		return
	end function
	subroutine name_QN(dimen,QN,degeneracy,rule)
		type(SymDimensionName),intent(inout)::dimen
		integer,intent(in)::rule,degeneracy(:)
		real*4,intent(in)::QN
		integer::i,lendegen
		real*4::QN_i
		lendegen=size(degeneracy)
		if(allocated(dimen%QN))then
			deallocate(dimen%QN)
		end if
		if(allocated(dimen%degeneracy))then
			deallocate(dimen%degeneracy)
		end if
		allocate(dimen%degeneracy(lendegen))
		allocate(dimen%QN(lendegen))
		QN_i=-QN
		do i=1,size(degeneracy)
			dimen%degeneracy(i)=degeneracy(i)
			dimen%QN(i)=QN_i
			QN_i=QN_i+1
		end do
		dimen%rule=rule
		dimen%max_QN=QN
		return
	end subroutine
!*****************************************************
!  input A0 output 'A' , 0
!*******  function or subroutine for name   **************
	subroutine cleanSymDimensionName(Dimen)
		type(SymDimension),intent(inout)::	Dimen
		if(allocated(Dimen%DimName))then
			deallocate(Dimen%DimName)
		end if
		Dimen%nameflag=.false.
		return
	end subroutine
!  input asd.ad output 'asd' and 'ad', SymchartoName_log=.true.
!	!  input asd output 'asd' , '0'  and SymchartoName_log=.false.
	logical function SymchartoName_log(TensorName,DimenName,w_)
		character(len=*),intent(in)::w_
		character(len=len(trim(adjustl(w_))))::w
		character(len=len_of_Name),intent(out)::DimenName,TensorName
		integer::lenw,lenw2
		w=trim(adjustl(w_))
		lenw=index(w,indexsymbol)
		lenw2=len(trim(w))
		if(lenw.eq.0)then!input asd
			SymchartoName_log=.false.
			TensorName=w
			DimenName='0'
			return
		end if
		SymchartoName_log=.true.
		!input 'A2.1'(if indexsymbol='.') 
		TensorName=w(1:lenw-1)
		DimenName=w(lenw+1:lenw2)
		return
	end function
!*******  function or subroutine for name   **************
	type(SymDimensionName) function SymNameinit1(TensorName,DimenName)
		character(len=*),intent(in)::TensorName
		character(len=*),intent(in)::DimenName
		SymNameinit1%TensorName=TensorName
		SymNameinit1%DimenName=DimenName
		return
	end function
	type(SymDimensionName) function SymNameinit2(TensorName,DimenName)
		character(len=*),intent(in)::TensorName
		integer,intent(in)::DimenName
		SymNameinit2%TensorName=TensorName
		SymNameinit2%DimenName=DimenName
		return
	end function
!*******  function or subroutine for name   **************
	subroutine SymcopyName(Nameout,Namein)
		type(SymDimensionName),intent(inout)::Nameout
		type(SymDimensionName),intent(in)::Namein
		Nameout%TensorName=Namein%TensorName
		Nameout%DimenName=Namein%DimenName
		if(allocated(Nameout%QN))then
			deallocate(Nameout%QN)
			deallocate(Nameout%degeneracy)
		end if
		if(allocated(Namein%QN))then
			allocate(Nameout%QN(size(Namein%QN)))
			allocate(Nameout%degeneracy(size(Namein%degeneracy)))
			Nameout%QN=Namein%QN
			Nameout%degeneracy=Namein%degeneracy
			Nameout%rule=Namein%rule
			Nameout%max_QN=Namein%max_QN
		end if
		return
	end subroutine		
!*******  function or subroutine for name   **************
	subroutine SymcopyNameArray(Nameout,Namein)
		type(SymDimensionName),intent(inout)::Nameout(:)
		type(SymDimensionName),intent(in)::Namein(:)
		integer::Lenout,lenIN,i
		lenIN=size(Namein)
		if(size(Nameout).lt.lenIN)then
			write(*,*)"ERROR in assignment of two Name array "
			write(*,*)"Name1(:)=Name2(:),size(Name1)<size(Name2)"
			write(*,*)size(Nameout),lenIN
			stop
		end if
		do i=1,lenIN
			Nameout(i)=Namein(i)
		end do
		return
	end subroutine	
!Name=w
!w='A','A1','A1_1'
	subroutine SymcharName(outName,w)
		type(SymDimensionName),intent(inout)::outName
		character(len=*),intent(in)::w
		logical::flag
		flag=SymchartoName_log(outName%TensorName,outName%DimenName,w)
		return
	end subroutine
	subroutine SymNamechar(w,inName)
		type(SymDimensionName),intent(in)::inName
		character(len=*),intent(inout)::w
		w=inName%TensorName+indexsymbol+inName%DimenName
		return
	end subroutine
	
!	input array ,for example (/'A1.12','A1.2 ','A1.3 '/)
!	output a array of SymDimensionName
!	the array input should be the same length
!	'A1.12' length is 5, 'A1.2 ' is 5 too
!name(:)=w(:)
	 subroutine SymcharNameArray(R,w)
	 	type(SymDimensionName),intent(inout)::R(:)
		character(len=*),intent(in)::w(:)
		integer::lenw,i
		lenw=size(w)
		if(size(R).lt.lenw)then
			write(*,*)"ERROR in assignment of two Name array "
			write(*,*)"R(:)=w(:),size(R)<size(w)"
			write(*,*)size(R),lenw
			stop
		end if
		do i=1,lenw
			R(i)=w(i)
		end do
		return
	end subroutine
!character(:)=SymDimensionName(:)
	subroutine SymNamecharArray(w,inName)
	 	type(SymDimensionName),intent(in)::inName(:)
		character(len=*),intent(inout)::w(:)
		integer::lenw,i
		lenw=size(inName)
		if(size(w).lt.lenw)then
			write(*,*)"ERROR in assignment of two Name array "
			write(*,*)"R(:)=w(:),size(R)<size(w)"
			write(*,*)size(w),lenw
			stop
		end if
		do i=1,lenw
			w(i)=inName(i)
		end do
		return
	end subroutine
!outName(:)=inName(:),will allocate outName	
	subroutine SymcopyNameAllocate(outName,inName)
		type(SymDimensionName),allocatable,intent(inout)::outName(:)
		type(SymDimensionName),intent(in)::inName(:)
		integer::Lenout,lenIN,i
		lenIN=size(inName)
		if(allocated(outName))then
			if(size(outName).ne.lenIN)then
				deallocate(outName)
				allocate(outName(lenIN))
			end if
		else
			allocate(outName(lenIN))
		end if
		do i=1,lenIN
			outName(i)=inName(i)
		end do
	end subroutine
	subroutine SymprintQN(na)
		type(SymDimensionName),intent(in)::na
		write(*,*)"quantum number(QN):",na%max_QN
		write(*,*)"rule",na%rule
		write(*,*)"QN_i :"
		write(*,*)na%QN
		write(*,*)"degeneracy:"
		write(*,*)na%degeneracy
		return
	end subroutine
	subroutine SymprintName(na)
		type(SymDimensionName),intent(in)::na
		character*1000::w
		w=na
		write(*,*)trim(adjustl(w))
		return
	end subroutine	
	
	
	
	
!*******  function or subroutine for name   **************
!set a name to the diemsnion
!w may be something like 'A' or 'A1'
	subroutine SymDimName(dimen,w)
		type(SymDimension),intent(inout)::dimen
		character(len=*),intent(in)::w
		character(len=len_of_Name)::TensorName,DimenName
		integer::lenD,i
		logical::flag
		flag=SymchartoName_log(TensorName,DimenName,w)
		if(flag)then
			write(*,*)"ERROR in Set the name to the Tensor"
			write(*,*)"input:",w
			write(*,*)"could not contains the indexsymbol:",indexsymbol
			stop
		end if
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension"
			stop
		end if
		lenD=size(dimen%DimData)
		if(allocated(dimen%DimName)) then
			if(lenD.ne.size(dimen%DimName)) then
				deallocate(dimen%DimName)
				allocate(dimen%DimName(lenD))
			end if
		else
			allocate(dimen%DimName(lenD))
		end if
		if(dimen%nameflag)then!There are name in the SymDimension,rename
			do i=1,lenD
				dimen%DimName(i)%TensorName=TensorName
			end do
		else
			do i=1,lenD
				dimen%DimName(i)=SymNameinit(TensorName,i)
			end do
			dimen%nameflag=.true.
		end if
		return
	end subroutine
!subroutine or function for quantum number
	subroutine set_QN1(dimen,ith,QN,degeneracy,rule)
		type(SymDimension),intent(inout)::dimen
		integer,intent(in)::rule,ith,degeneracy(:)
		real*4,intent(in)::QN
		integer::i,lendegen
		real*4::QN_i
		if(.not.dimen%nameflag)then
			call SymDimName(dimen,'0')
		end if
		lendegen=size(degeneracy)
		if(lendegen.ne.SymDim_i(dimen,ith))then
			write(*,*)"ERROR in set_QN"
			stop
		end if
		if(allocated(dimen%DimName(ith)%QN))then
			deallocate(dimen%DimName(ith)%QN)
		end if
		if(allocated(dimen%DimName(ith)%degeneracy))then
			deallocate(dimen%DimName(ith)%degeneracy)
		end if
		allocate(dimen%DimName(ith)%degeneracy(lendegen))
		allocate(dimen%DimName(ith)%QN(lendegen))
		QN_i=-QN
		do i=1,size(degeneracy)
			dimen%DimName(ith)%degeneracy(i)=degeneracy(i)
			dimen%DimName(ith)%QN(i)=QN_i
			QN_i=QN_i+1
		end do
		dimen%DimName(ith)%rule=rule
		dimen%DimName(ith)%max_QN=QN
		return
	end subroutine
!QN and Deg
!0.5:       1  1
!1  :      1  2  1
!1.5:     1  3  3  1
!2  :    1  4  6  4  1		
!2.5:   1 5 10 10  5  1 
!3  :  1 6 15 20 15 6  1
!3.5: 1 7 21 35 35 21 7 1
!These QNs or N*Deg give a good result on MPS
!input QN,Num
!output a random symTensor with degeneracy deg(i)=Num*C_{N-1}^{i-1},i=0,N-1
!if noDeg, all the degeneracy will be 1
	subroutine set_QN2(dimen,ith,QN,lendegen,rule,Num)
		type(SymDimension),intent(inout)::dimen
		integer,intent(in)::rule,ith,lendegen
		real*4,intent(in)::QN
		integer,optional,intent(in)::Num
		integer::i,N
		real*4::QN_i
		if(.not.dimen%nameflag)then
			call SymDimName(dimen,'0')
		end if
		if(lendegen.ne.SymDim_i(dimen,ith))then
			write(*,*)"ERROR in set_QN"
			stop
		end if
		if(allocated(dimen%DimName(ith)%QN))then
			deallocate(dimen%DimName(ith)%QN)
		end if
		if(allocated(dimen%DimName(ith)%degeneracy))then
			deallocate(dimen%DimName(ith)%degeneracy)
		end if
		allocate(dimen%DimName(ith)%degeneracy(lendegen))
		allocate(dimen%DimName(ith)%QN(lendegen))
		QN_i=-QN
		do i=1,lendegen
			if(present(Num))then
				dimen%DimName(ith)%degeneracy(i)=Num*C_N_I(lendegen-1,i-1)
			else
				dimen%DimName(ith)%degeneracy(i)=1
			end if
			dimen%DimName(ith)%QN(i)=QN_i
			QN_i=QN_i+1
		end do
		dimen%DimName(ith)%rule=rule
		dimen%DimName(ith)%max_QN=QN
		return
	end subroutine
!C_N^i
	integer function C_N_I(N,i)
		integer,intent(in)::N,i
		integer::j
		if(i.eq.0)then
			C_N_I=1
			return
		end if
		C_N_I=1
		do j=1,i
			C_N_I=C_N_I*(N-j+1)/j
		end do
		return
	end function
!***************   cleanSymDimension   *****************
	subroutine cleanSymDimension(Dimen)
		type(SymDimension),intent(inout)::	Dimen
		Dimen%boundarysize=0
		Dimen%Dimsize=0
		if(allocated(Dimen%boundary)) then
			deallocate(Dimen%boundary)
		end if
		if(allocated(Dimen%DimData)) then
			deallocate(Dimen%DimData)
		end if
		call cleanSymDimensionName(Dimen)
		return
	end subroutine
!****************Symprint the information of SymDimension ****************************************	
	subroutine SymDprint(Dimen)
		type(SymDimension),intent(in) ::Dimen
		integer,allocatable :: dimenVec(:)
		CHARACTER*5000,allocatable::w(:)
		integer::i
		write(*,*) "***   START   ***"
		write(*,*) Dimen%DimData
		call SymgetDim(dimenVec,Dimen)
		write(*,*) "			--				"
		write(*,*)	dimenVec
		if(Dimen%nameflag)then
			allocate(w(size(Dimen%DimName)))
			w=Dimen%DimName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,size(Dimen%DimName))
			call SymRule(dimenVec,dimen)
			write(*,*)"quantum number are:"
			write(*,*)(dimenVec(i)*Dimen%DimName(i)%max_QN,i=1,size(Dimen%DimName))
			write(*,*)"Degeneracy:"
			do i=1,size(Dimen%DimName)
				write(*,*)Dimen%DimName(i)%QN
				write(*,*)Dimen%DimName(i)%Degeneracy
				write(*,*)"  ---===---  "
			end do
		end if
		write(*,*) "***   END   ***"
	end subroutine
				
	subroutine SymDprint0(Dimen)
		type(SymDimension),intent(in) ::Dimen
		integer,allocatable :: dimenVec(:)
		CHARACTER*5000,allocatable::w(:)
		integer::i
		write(*,*) "***   SymDimension Data    ***"
		write(*,*) Dimen%DimData
		call SymgetDim(dimenVec,Dimen)
		write(*,*)	dimenVec
		if(Dimen%nameflag)then
			allocate(w(size(Dimen%DimName)))
			w=Dimen%DimName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,size(Dimen%DimName))
			call SymRule(dimenVec,dimen)
			write(*,*)"quantum number are:"
			write(*,*)(dimenVec(i)*Dimen%DimName(i)%max_QN,i=1,size(Dimen%DimName))
			write(*,*)"Degeneracy:"
			do i=1,size(Dimen%DimName)
				write(*,*)Dimen%DimName(i)%QN
				write(*,*)Dimen%DimName(i)%Degeneracy
				write(*,*)"  ---===---  "
			end do
		end if
		write(*,*) "***   SymDimension END   ***"
	end subroutine
				
	subroutine SymDprint2(Dimen)
		type(SymDimension),intent(in) ::Dimen
		integer,allocatable :: dimenVec(:)
		CHARACTER*500,allocatable::w(:)
		integer::i
		write(*,*) "***   START   ***"
		write(*,*) "Dimsize"
		write(*,*) Dimen%Dimsize
		write(*,*) "boundarysize"
		write(*,*) Dimen%boundarysize
		write(*,*) "boundary"
		write(*,*) Dimen%boundary
		write(*,*) "DimData"
		write(*,*) Dimen%DimData
		call SymgetDim(dimenVec,Dimen)
		write(*,*) "			--				"
		write(*,*)	dimenVec
		if(Dimen%nameflag)then
			allocate(w(size(Dimen%DimName)))
			w=Dimen%DimName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,size(Dimen%DimName))
			call SymRule(dimenVec,dimen)
			write(*,*)"quantum number are:"
			write(*,*)(dimenVec(i)*Dimen%DimName(i)%max_QN,i=1,size(Dimen%DimName))
			write(*,*)"Degeneracy:"
			do i=1,size(Dimen%DimName)
				write(*,*)Dimen%DimName(i)%QN
				write(*,*)Dimen%DimName(i)%Degeneracy
				write(*,*)"  ---===---  "
			end do
		end if
		write(*,*) "***   END   ***"
		write(*,*) " "
	end subroutine	
	subroutine SymprintdimenQN(na,ith)
		type(SymDimension),intent(in)::na
		integer,intent(in)::ith
		write(*,*)"quantum number(QN):",na%DimName(ith)%max_QN
		write(*,*)"rule",na%DimName(ith)%rule
		write(*,*)"QN_i :"
		write(*,*)na%DimName(ith)%QN
		write(*,*)"degeneracy:"
		write(*,*)na%DimName(ith)%degeneracy
		return
	end subroutine	
	type(SymDimension) function SymDimensioniniti(boundarysize,Dimsize,boundary,DimData)
		integer,intent(in) :: boundarysize
		integer,intent(in) :: Dimsize
		integer,intent(in) :: boundary(:)
		integer,intent(in) :: DimData(:)
		allocate(SymDimensioniniti%boundary(size(boundary,1)))
		allocate(SymDimensioniniti%DimData(size(DimData,1)))
		SymDimensioniniti%boundary=boundary
		SymDimensioniniti%DimData=DimData
		SymDimensioniniti%boundarysize=boundarysize
		SymDimensioniniti%Dimsize=Dimsize
		SymDimensioniniti%Nameflag=.false.
		return
	end function
	subroutine SymDimInitialization(Dimen,DimData)
		type(SymDimension),intent(inout) ::Dimen
		integer,intent(in) :: DimData(:)
		integer::i
		call cleanSymDimension(Dimen)
		Dimen%Dimsize=size(DimData)
		Dimen%boundarysize=Dimen%Dimsize+1
		allocate(Dimen%boundary(Dimen%boundarysize))
		do i=1,Dimen%boundarysize
			Dimen%boundary(i)=i-1
		end do
		allocate(Dimen%DimData(Dimen%Dimsize))
		Dimen%DimData=DimData
		Dimen%nameflag=.false.
		return
	end subroutine
	
	
	
	subroutine SymDimInitialization2(Dimen,Dimen2)
		type(SymDimension),intent(inout) ::Dimen
		type(SymDimension),intent(in) ::Dimen2
		call cleanSymDimension(Dimen)
		Dimen%Dimsize=Dimen2%Dimsize
		Dimen%boundarysize=Dimen2%boundarysize
		allocate(Dimen%DimData(size(Dimen2%Dimdata)))
		allocate(Dimen%boundary(size(Dimen2%boundary)))
		Dimen%boundary=Dimen2%boundary
		Dimen%DimData=Dimen2%DimData
		Dimen%nameflag=Dimen2%nameflag
		if(Dimen2%nameflag)then
			call SymcopyNameAllocate(Dimen%DimName,Dimen2%DimName)
		end if
		return
	end subroutine
	integer	function SymDimSize(Dimen)
		type(SymDimension),intent(in) :: Dimen
		SymDimSize=Dimen%dimSize
		return
	end function
!*******  function or subroutine for name   **************		
!output the TensorName of the ith SymDimension
	CHARACTER(len=len_of_Name) function  SymoutNameTen(dimen,ith)
		type(SymDimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(.not.Dimen%nameflag)then
			write(*,*)"There is no name in the SymDimension"
			stop
		end if
		if(ith.gt.size(Dimen%DimName))then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,size(Dimen%DimName)
			stop
		end if
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"NameID is use in original SymDimension"
			stop
		end if
		 SymoutNameTen=Dimen%DimName(ith)%TensorName
		return
	end function
!*******  function or subroutine for name   **************		
!output the index of the ith SymDimension
	CHARACTER(len=len_of_Name+len_of_Name) function SymoutName(dimen,ith)
		type(SymDimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(.not.Dimen%nameflag)then
			write(*,*)"There is no name in the SymDimension"
			stop
		end if
		if(ith.gt.size(Dimen%DimName))then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,size(Dimen%DimName)
			stop
		end if
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"indexID is use in original SymDimension"
			stop
		end if
		SymoutName=Dimen%DimName(ith)
		return
	end function

!*******  function or subroutine for name   **************		
!find the index,whose name is w	,output the order of it in the SymDimension
!If can not find , output 0
	integer function SymNameorder2(dimen,w)
		type(SymDimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w
		integer::i
		type(SymDimensionName)::nam
		if(.not.Dimen%nameflag)then
			write(*,*)"There is no name in the SymDimension,SymNameorder2"
			stop
		end if
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"SymNameorder2 is use in original SymDimension"
			stop
		end if
		nam=w
		do i=1,size(dimen%DimName)
			if(dimen%DimName(i).equ.nam)then
				SymNameorder2=i
				return
			end if
		end do
		SymNameorder2=0
		return
	end function
	integer function SymNameorder3(dimen,SymDimName)
		type(SymDimension),intent(in) :: Dimen
		integer::i
		type(SymDimensionName),intent(in)::SymDimName
		if(.not.Dimen%nameflag)then
			write(*,*)"There is no name in the SymDimension,SymNameorder3"
			stop
		end if
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"SymNameorder2 is use in original SymDimension"
			stop
		end if
		do i=1,size(dimen%DimName)
			if(dimen%DimName(i).equ.SymDimName)then
				SymNameorder3=i
				return
			end if
		end do
		SymNameorder3=0
		return
	end function
	function SymNameorderAllocateArray(dimen,w)
		integer,allocatable::SymNameorderAllocateArray(:)
		type(SymDimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w(:)
		integer::i,lenn
		if(.not.Dimen%nameflag)then
			write(*,*)"There is no name in the SymDimension,SymNameorderAllocateArray"
			stop
		end if
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"SymNameorder2 is use in original SymDimension"
			stop
		end if
		lenn=size(w)
		if(lenn.gt.dimen%Dimsize)then
			write(*,*)"ERROR i SymNameorderAllocateArray"
			write(*,*)lenn,dimen%Dimsize
			stop
		end if
		allocate(SymNameorderAllocateArray(lenn))
		do i=1,lenn
			SymNameorderAllocateArray(i)=SymNameorder(dimen,w(i))
		end do
		return
	end function
	function SymNameorderArray(dimen,w)
		type(SymDimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w(:)
		integer::SymNameorderArray(size(w))
		integer::i,lenn
		if(.not.Dimen%nameflag)then
			write(*,*)"There is no name in the SymDimension,SymNameorderArray"
			stop
		end if
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"SymNameorder2 is use in original SymDimension"
			stop
		end if
		lenn=size(w)
		if(lenn.gt.dimen%Dimsize)then
			write(*,*)"ERROR i SymNameorderArray"
			write(*,*)lenn,dimen%Dimsize
			stop
		end if
		do i=1,lenn
			SymNameorderArray(i)=SymNameorder(dimen,w(i))
		end do
		return
	end function

	
!*******  function or subroutine for name   **************		
!Find all the SymDimensionNames, whose indexname is oldname
!change them to newname
!If Cannot find it , do nothing
	subroutine SymResetname1(dimen,oldname,newname)
		type(SymDimension),intent(inout) :: dimen
		CHARACTER(len=*),intent(in)::oldname,newname
		integer::i
		CHARACTER(len=len_of_Name)::oldTensorName,olddimenName
		CHARACTER(len=len_of_Name)::newTensorName,newdimenName
		logical::oldfullname,newfullname
		if(.not.Dimen%nameflag)then
			write(*,*)"There is no name in the SymDimension,resetindexname"
			stop
		end if
		oldfullname=SymchartoName_log(oldTensorName,olddimenName,oldname)
		!(TensorName,DimenName,w_)
		newfullname=SymchartoName_log(newTensorName,newdimenName,newname)
		if(oldfullname.neqv.newfullname)then
			write(*,*)"ERROR in resetIndexname1"
			stop
		end if
		if(newfullname)then
			do i=1,size(dimen%DimName)
				if(dimen%DimName(i)%TensorName.equ.oldTensorName)then
					if(dimen%DimName(i)%dimenName.equ.olddimenName)then
						dimen%DimName(i)%TensorName=newTensorName
						dimen%DimName(i)%dimenName=newdimenName
					end if
				end if
			end do
		else
			do i=1,size(dimen%DimName)
				if(dimen%DimName(i)%TensorName.equ.oldTensorName)then
						dimen%DimName(i)%TensorName=newTensorName
				end if
			end do
		end if
		return
	end subroutine
!*******  function or subroutine for name   **************		
!set the NameID and Indexname of the ith index
	subroutine SymResetname2(dimen,ith,newname)
		type(SymDimension),intent(inout) :: dimen
		CHARACTER(len=*),intent(in)::newname
		integer,intent(in)::ith
		CHARACTER(len=len_of_Name)::newTensorName,newdimenName
		logical::fullname
		integer::i,lenName
		if(ith.gt.size(Dimen%DimData))then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,size(Dimen%DimName)
			stop
		end if
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"resetIndexNameID is use in original SymDimension"
			stop
		end if
		fullname=SymchartoName_log(newTensorName,newdimenName,newname)
		if(.not.Dimen%nameflag)then!There is no name in SymDimension
			lenName=size(dimen%DimData)
			allocate(dimen%DimName(lenName))
			do i=1,lenName
				dimen%DimName(i)%TensorName=newTensorName
				dimen%DimName(i)%dimenName=i
			end do
			if(fullname)then
				dimen%DimName(ith)%dimenName=newdimenName
			end if
			Dimen%nameflag=.true.
			return
		end if
		
		if(fullname)then
			dimen%DimName(ith)%TensorName=newTensorName
			dimen%DimName(ith)%dimenName=newdimenName
		else
			dimen%DimName(ith)%TensorName=newTensorName
		end if
		return
	end subroutine

	integer function  SymouttotalData(Dimen)
		type(SymDimension),intent(in) :: Dimen
		 SymouttotalData=product(Dimen%DimData)
		return
	end function				
!	return the inde SymDimension	,outpout in a integer
	integer function SymDim_i(Dimen,inde)
		type(SymDimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		integer :: i,D1
		if(inde.gt.dimen%dimSize) then 
			write(*,*) "ERROR in SymDim_i"
			write(*,*)"stop"
			stop
		end if
		D1=1
		do i=Dimen%boundary(inde) + 1,Dimen%boundary(inde+1)
			D1=D1*Dimen%DimData(i)
		end do
		SymDim_i=D1
	end function 
!	return  all the  SymDimension	,outpout in a vector	
	subroutine SymgetDim(dimenVec,Dimen)
		integer,allocatable,intent(inout) :: dimenVec(:)
		type(SymDimension),intent(in) :: Dimen
		integer :: i
		if(allocated(dimenVec)) then
			if(Dimen%DimSize.ne.size(dimenVec)) then 
				deallocate(dimenVec)
				allocate(dimenVec(Dimen%DimSize))
			end if
		else
			allocate(dimenVec(Dimen%DimSize))
		end if
		do i=1,Dimen%Dimsize
			dimenVec(i)=SymDim_i(Dimen,i)
		end do
	end 	subroutine	
	subroutine SymDimToint(dimenVec,Dimen)
		integer,intent(inout) :: dimenVec(:)
		type(SymDimension),intent(in) :: Dimen
		integer :: i
		do i=1,Dimen%Dimsize
			dimenVec(i)=SymDim_i(Dimen,i)
		end do
	end 	subroutine	

!	return the inde  SymDimension	,outpout in a vector of the dimenison    
! If do the contract, onedimenison will have more than one value
! [2,3,4,5]	-->contract the 2,3 index -->[2,(3,4),5]!the SymDimension is 2,12,5
! then SymgetSubDim(Dimen,2,dimenVec)==>dimenVec=[3,4]
	subroutine SymgetSubDim(Dimen,inde,dimenVec)
		type(SymDimension),intent(in) :: Dimen
		integer,allocatable,intent(inout) :: dimenVec(:)
		integer,intent(in) :: inde
		integer::i,j,Dlen
		Dlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
		if(inde.gt.dimen%dimSize) then 
			write(*,*) "ERROR in SymgetSubDim"
			return
		end if
		if(allocated(dimenVec)) then
			if(Dlen.ne.size(dimenVec)) then 
				deallocate(dimenVec)
				allocate(dimenVec(Dlen))
			end if
		else
			allocate(dimenVec(Dlen))
		end if
		j=1
		do i=Dimen%boundary(inde) + 1,Dimen%boundary(inde+1)
			dimenVec(j)=Dimen%DimData(i)
			j=j+1
		end do
	end subroutine
!*******  function or subroutine for name   **************	
	subroutine getSubSymDimName(Dimen,inde,dimenName)
		type(SymDimension),intent(in) :: Dimen
		type(SymDimensionName),allocatable,intent(inout) :: dimenName(:)
		integer,intent(in) :: inde
		integer::i,j,Dlen
		Dlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
		if(inde.gt.dimen%dimSize) then 
			write(*,*) "ERROR in SymgetSubDim"
			return
		end if
		if(allocated(dimenName)) then
			if(Dlen.ne.size(dimenName)) then 
				deallocate(dimenName)
				allocate(dimenName(Dlen))
			end if
		else
			allocate(dimenName(Dlen))
		end if
		j=1
		do i=Dimen%boundary(inde) + 1,Dimen%boundary(inde+1)
			dimenName(j)=Dimen%DimName(i)
			j=j+1
		end do
	end subroutine
!input SymDimension [1,1,2,1,3,1,1,4,1]
!output dimenison [2,3,4]
! if input [1,1,1,1,1]
!  ouput [1] without name
	type(SymDimension) function  SymRNDim(dimen) 
		type(SymDimension),intent(inout) :: Dimen
		integer::i,lenD,lenNewD
		integer,allocatable::Dimindex(:)
		if(.not.Symif_original_dim(dimen)) then
			write(*,*)"ERROR IN  SymRNDim,in Diemnsion.f90"
			stop
		end if
		lenD=size(Dimen%Dimdata)
		allocate(Dimindex(lenD))
		lenNewD=0
		do i=1,lenD
			if(Dimen%Dimdata(i).ne.1)then
				lenNewD=lenNewD+1
				Dimindex(lenNewD)=i
			end if
		end do
		if(lenNewD.eq.0)then
			 SymRNDim=(/1/)
			return
		end if
			
		allocate( SymRNDim%DimData(lenNewD))
		do i=1,lenNewD
			 SymRNDim%DimData(i)=dimen%DimData(Dimindex(i))
		end do
		 SymRNDim%Dimsize=lenNewD
		 SymRNDim%boundarysize= SymRNDim%Dimsize+1
		allocate( SymRNDim%boundary( SymRNDim%boundarysize))
		do i=1, SymRNDim%boundarysize
			 SymRNDim%boundary(i)=i-1
		end do
		 SymRNDim%NameFlag=dimen%NameFlag
		if(dimen%NameFlag)then
			allocate( SymRNDim%DimName(lenNewD))
			do i=1,lenNewD
				 SymRNDim%DimName(i)=dimen%DimName(Dimindex(i))
			end do
		end if
		return
	end function
!*******  function or subroutine for name   **************	
!	return the inde  SymDimension	,outpout in a type(dimenison)
	type(SymDimension) function  SymgetSubDim2(Dimen,inde)
		type(SymDimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		integer::boundarysize,boundary(2),Dimsize
		integer,allocatable :: Dimdata(:)
		integer::i,j,Dlen
		call SymgetSubDim(Dimen,inde,Dimdata)
		boundary(1)=0
		boundary(2)=size(Dimdata)
		boundarysize=2
		SymgetSubDim2=SymDimensioniniti(boundarysize,1,boundary,DimData)
		if(Dimen%nameFlag)then
			call getSubSymDimName(Dimen,inde,SymgetSubDim2%DimName)
			SymgetSubDim2%NameFlag=.true.
		else
			SymgetSubDim2%NameFlag=.false.
		end if
	end function
	type(SymDimension) function  SymgetSubDim2_name(Dimen,w)
		type(SymDimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w
		integer::inde
		inde=SymNameorder(Dimen,w)
		SymgetSubDim2_name=SymgetSubDim2(Dimen,inde)
		return
	end function
!*******  function or subroutine for name   **************	
! The type of input are both SymDimension	
!If dimen have a name but Dimen2 do not or Dimen2 have but not dimen
!then the name for the one with no name will be ['0',0,i]
	type(SymDimension) function  SymDimadd(Dimen,Dimen2)
		type(SymDimension),intent(in) :: Dimen,Dimen2
		!integer::boundarysize,Dimsize
		!integer,allocatable :: Dimdata(:),boundary(:)
		integer::i,j,l1,l2
		SymDimadd%boundarysize=Dimen%boundarysize+Dimen2%boundarysize-1
		allocate(SymDimadd%boundary(SymDimadd%boundarysize))
		do i=1,Dimen%boundarysize
			SymDimadd%boundary(i)=Dimen%boundary(i)
		end do
		l1=Dimen%boundary(Dimen%boundarysize)
		do i=1,Dimen2%boundarysize-1
			SymDimadd%boundary(i+Dimen%boundarysize)=Dimen2%boundary(i+1)+l1
		end do
		l1=size(Dimen%DimData)
		l2=size(Dimen2%DimData)
		allocate(SymDimadd%Dimdata(l1+l2))
		SymDimadd%Dimdata(:l1)=Dimen%DimData
		SymDimadd%Dimdata(l1+1:)=Dimen2%DimData
		SymDimadd%Dimsize=Dimen%Dimsize+Dimen2%Dimsize
		if(Dimen%nameFlag.or.Dimen2%nameFlag)then
			allocate(SymDimadd%DimName(l1+l2))
			if(Dimen%nameFlag)then
				SymDimadd%DimName(1:l1)=Dimen%DimName
			else
				do i=1,l1
					SymDimadd%DimName(i)=SymNameinit('0',i)
				end do
			end if
			if(Dimen2%nameFlag)then
				SymDimadd%DimName(l1+1:)=Dimen2%DimName
			else
				do i=l1+1,l1+l2
					SymDimadd%DimName(i)=SymNameinit('0',i-l1)
				end do
			end if
			SymDimadd%nameflag=.true.
		else
			SymDimadd%nameflag=.false.
		end if
	!	SymDimadd=SymDimensioniniti(boundarysize,Dimen%Dimsize+Dimen2%Dimsize,boundary,DimData)
	end function
! The type of input are one SymDimension and the other vector	
	type(SymDimension) function  SymDimadd2(Dimen,Dimenvec)
		type(SymDimension),intent(in) :: Dimen
		integer,intent(in) :: Dimenvec(:)
		!integer::boundarysize,Dimsize
		!integer,allocatable :: Dimdata(:),boundary(:)
		integer::i,j,l1,l2
		SymDimadd2%boundarysize=Dimen%boundarysize+size(Dimenvec)
		allocate(SymDimadd2%boundary(SymDimadd2%boundarysize))
		do i=1,Dimen%boundarysize
			SymDimadd2%boundary(i)=Dimen%boundary(i)
		end do
		l1=Dimen%boundary(Dimen%boundarysize)
		do i=1,size(Dimenvec)
			SymDimadd2%boundary(i+Dimen%boundarysize)=i+l1
		end do
		l1=size(Dimen%DimData)
		l2=size(Dimenvec)
		allocate(SymDimadd2%Dimdata(l1+l2))
		SymDimadd2%Dimsize=Dimen%Dimsize+size(Dimenvec)
		SymDimadd2%Dimdata(:l1)=Dimen%DimData
		SymDimadd2%Dimdata(l1+1:)=Dimenvec
		if(Dimen%nameFlag)then
			allocate(SymDimadd2%DimName(l1+l2))
			SymDimadd2%DimName(1:l1)=Dimen%DimName
			do i=l1+1,l1+l2
				SymDimadd2%DimName(i)=SymNameinit('0',i-l1)
			end do
			SymDimadd2%nameflag=.true.
		else
			SymDimadd2%nameflag=.false.
		end if
		!SymDimadd2=SymDimensioniniti(boundarysize,Dimen%Dimsize+size(Dimenvec),boundary,DimData)
	end function
	type(SymDimension) function  SymDimadd3(Dimenvec,Dimen)
		type(SymDimension),intent(in) :: Dimen
		integer,intent(in) :: Dimenvec(:)
		type(SymDimension)::Dimenvec_
		Dimenvec_=Dimenvec
		SymDimadd3=SymDimadd(Dimenvec_,Dimen)
	end function
	
	integer function SymgetSubDimlen(Dimen,inde)
		type(SymDimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		SymgetSubDimlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
	end function
!		1,2,3,4,..inde,inde+1,...rank-->1,2,3,4,..(inde*inde+1*..inde+num),...rank     	
!		if inde+num>boundarysize,it will constract all the index that larger than inde
	type(SymDimension) function  SymDimConstract(dimen,inde,num)
		type(SymDimension),intent(in) :: Dimen
		integer,intent(in):: inde,num
		!integer :: boundarysize
		integer :: i!,Dimsize
		!integer,allocatable :: boundary(:)
		if(inde.gt.dimen%Dimsize) then 
			write(*,*) "ERROR in  SymDimConstract"
			write(*,*)inde,dimen%Dimsize
			stop
		end if
		if((num.le.0).or.(inde.eq.dimen%Dimsize)) then
			 SymDimConstract=dimen
			return
		end if
		if(inde+num.ge.dimen%boundarysize) then
			 SymDimConstract%boundarysize=inde+1
			 SymDimConstract%Dimsize=inde
			allocate( SymDimConstract%boundary( SymDimConstract%boundarysize))
			 SymDimConstract%boundary(:inde)=dimen%boundary(:inde)
			 SymDimConstract%boundary(inde+1)=dimen%boundary(dimen%boundarysize)
			 allocate(SymDimConstract%DimData(size(dimen%DimData)))
			 SymDimConstract%DimData=dimen%DimData
			if(dimen%Nameflag)then
				call SymcopyNameAllocate( SymDimConstract%DimName,dimen%DimName)
				 SymDimConstract%Nameflag=.true.
			else
				 SymDimConstract%Nameflag=.false.
			end if
			return
		end if
		 SymDimConstract%boundarysize=dimen%boundarysize-num
		 SymDimConstract%Dimsize=dimen%Dimsize-num
		allocate( SymDimConstract%boundary( SymDimConstract%boundarysize))
		 SymDimConstract%boundary(:inde)=dimen%boundary(:inde)
		do i=inde+1, SymDimConstract%boundarysize
			 SymDimConstract%boundary(i)=dimen%boundary(i+num)
		end do
		 allocate(SymDimConstract%DimData(size(dimen%DimData)))
		 SymDimConstract%DimData=dimen%DimData
		if(dimen%Nameflag)then
			call SymcopyNameAllocate( SymDimConstract%DimName,dimen%DimName)
			 SymDimConstract%Nameflag=.true.
			! SymDimConstract=SymDimensioniniti_name(boundarysize,Dimsize,boundary,dimen%DimData,dimen%DimName)
		else
			 SymDimConstract%Nameflag=.false.
			! SymDimConstract=SymDimensioniniti(boundarysize,Dimsize,boundary,dimen%DimData)
		end if
		return
	end function
     	

     	
!		(inde,inde+1,..midindde,midindde+1,...rank)-->(inde,inde+1,..midindde),(midindde+1,...rank)	     		
!		if (inde,inde+1 ..midindde),	midindde larger then next elmemt,do nothing
!		for example:	
!			D=[2,2,(3,4,5),2,(3,4)],	boundary=[0,1,2,5,6,8]
!		 SymDimDecompose(D,3,2)=[2,2,(3,4),5,2,(3,4)],boundary=[0,1,2,4,5,6,8]
!		 SymDimDecompose(D,3,4) will do nothing
	type(SymDimension) function  SymDimDecompose(dimen,inde,midindde)
		type(SymDimension),intent(in) :: Dimen
		integer,intent(in):: inde,midindde
	!	integer :: boundarysize
		integer :: i!,Dimsize
	!	integer,allocatable :: boundary(:)
		if(dimen%boundary(inde) +midindde .ge.dimen%boundary(inde+1)) then
			 SymDimDecompose=dimen
			return
		end if
		 SymDimDecompose%boundarysize=dimen%boundarysize+1
		 SymDimDecompose%Dimsize=dimen%Dimsize+1
		allocate( SymDimDecompose%boundary( SymDimDecompose%boundarysize))
		 SymDimDecompose%boundary(:inde)=dimen%boundary(:inde)
		 SymDimDecompose%boundary(inde+1)=dimen%boundary(inde) + midindde
		do i=inde+1, SymDimDecompose%boundarysize-1
			 SymDimDecompose%boundary(i+1)=dimen%boundary(i)
		end do
		allocate(SymDimDecompose%DimData(size(dimen%DimData)))
		 SymDimDecompose%DimData=dimen%DimData
		if(dimen%Nameflag)then
			call SymcopyNameAllocate( SymDimDecompose%DimName,dimen%DimName)
			 SymDimDecompose%Nameflag=.true.
		!	 SymDimDecompose=SymDimensioniniti_name(boundarysize,Dimsize,boundary,dimen%DimData,dimen%DimName)
		else
			 SymDimDecompose%Nameflag=.false.
			! SymDimDecompose=SymDimensioniniti(boundarysize,Dimsize,boundary,dimen%DimData)
		end if
		
		end function
		
		type(SymDimension) function  SymDimDecomposeAll(dimen)
		type(SymDimension),intent(in) :: Dimen
	!	integer :: boundarysize
		integer :: i!,Dimsize
!		integer,allocatable :: boundary(:)
		if(dimen%boundarysize.eq.size(dimen%DimData)+1) then
			 SymDimDecomposeAll=dimen
			return
     	end if
		 SymDimDecomposeAll%boundarysize=size(dimen%DimData)+1
		 SymDimDecomposeAll%Dimsize=size(dimen%DimData)
		allocate( SymDimDecomposeAll%boundary( SymDimDecomposeAll%boundarysize))
		do i=1, SymDimDecomposeAll%boundarysize
			 SymDimDecomposeAll%boundary(i)=i-1
		end do
		allocate(SymDimDecomposeAll%DimData(size(dimen%DimData)))
		 SymDimDecomposeAll%DimData=dimen%DimData
		if(dimen%Nameflag)then
			call SymcopyNameAllocate( SymDimDecomposeAll%DimName,dimen%DimName)
			 SymDimDecomposeAll%Nameflag=.true.
		!	 SymDimDecomposeAll=SymDimensioniniti_name(boundarysize,Dimsize,boundary,dimen%DimData,dimen%DimName)
		else
		!	 SymDimDecomposeAll=SymDimensioniniti(boundarysize,Dimsize,boundary,dimen%DimData)
			 SymDimDecomposeAll%Nameflag=.false.
		end if
	end function
     
     			
	type(SymDimension) function  SymDimpermute(dimen,v)
		type(SymDimension),intent(in) :: Dimen
		integer,intent(in):: v(Dimen%dimSize)
		type(SymDimensionName),allocatable::dimenName(:)
		integer::i,datalen,subDlen,k!,boundarysize,Dimsize
	!	integer,allocatable :: boundary(:),DimData(:)
		integer,allocatable :: dimenVec(:)
		 SymDimpermute%Dimsize=Dimen%Dimsize
		 SymDimpermute%boundarysize=Dimen%boundarysize
		datalen=size(Dimen%DimData)
		allocate( SymDimpermute%boundary( SymDimpermute%boundarysize))
		allocate( SymDimpermute%DimData(datalen))
		if(dimen%NameFlag)then
			allocate( SymDimpermute%DimName(datalen))
		end if
		 SymDimpermute%boundary(1)=0
		k=1
		do i=1,Dimen%dimSize
			subDlen=SymgetSubDimlen(dimen,v(i))
			 SymDimpermute%boundary(i+1)= SymDimpermute%boundary(i)+subDlen
			call SymgetSubDim(Dimen,v(i),dimenVec)
			 SymDimpermute%DimData(k:subDlen+k-1)=dimenVec
			if(dimen%NameFlag)then
				call getSubSymDimName(Dimen,v(i),dimenName)
				 SymDimpermute%DimName(k:subDlen+k-1)=dimenName
			end if
			k=k+subDlen
		end do
		 SymDimpermute%Nameflag=dimen%NameFlag
		! SymDimpermute=SymDimensioniniti(boundarysize,Dimsize,boundary,DimData)
	end function


!******************************************************************************
!if flag=1,inv[1,2,3,4,...,ith,ith+1,...]->[ith,1,2,3,..,ith-1,ith+1,ith+2,...]
!if flag=0,inv[1,2,3,4,...,ith,ith+1,...]->[2,3,..,ith-1,ith,1,ith+1,ith+2,...]
	subroutine  Sympermuteorder(outv,inv_,ith,flag)
		integer,allocatable,intent(inout)::outv(:)
		integer,intent(in)::inv_(:)
		integer,intent(in)::ith,flag
		integer,allocatable::inv(:)
		integer::lenv,i
		call assignments(inv,inv_)
		lenv=size(inv)
		if(allocated(outv)) then
			if(lenv.ne.size(outv)) then
				deallocate(outv)
				allocate(outv(lenv))
			end if
		else
			allocate(outv(lenv))
		end if
		if(flag.eq.1) then
			outv(1)=inv(ith)
			do i=1,ith-1
				outv(i+1)=inv(i)
			end do
			do i=ith+1,lenv
				outv(i)=inv(i)
			end do
		else if(flag.eq.0) then
			do i=2,ith
				outv(i-1)=inv(i)
			end do
			outv(ith)=inv(1)
			do i=ith+1,lenv
				outv(i)=inv(i)
			end do
		else
			write(*,*)"ERROR in  Sympermuteorder"
			stop
		end if
		return
	end subroutine
	logical function Symequal_of_dim(dim1,dim2)
		type(SymDimension),intent(in) :: dim1,dim2
		integer,allocatable :: dimenVec1(:)
		integer,allocatable :: dimenVec2(:)
		call SymgetDim(dimenVec1,dim1)
		call SymgetDim(dimenVec2,dim2)
		Symequal_of_dim=equal_of_array(dimenVec1,dimenVec2)
	end function
		logical function  Symequal_name2(dimname,dimname2)
		type(SymDimensionName),intent(in)::dimname,dimname2
		if(trim(adjustl(dimname%TensorName)).ne.trim(adjustl(dimname2%TensorName)))then
			 Symequal_name2=.false.
			return
		end if
		if(trim(adjustl(dimname%dimenName)).ne.trim(adjustl(dimname2%dimenName)))then
			 Symequal_name2=.false.
			return
		end if
		 Symequal_name2=.true.
		return
	end function
	logical function  Symequal_name1(dimname,w)
		type(SymDimensionName),intent(in)::dimname
		CHARACTER(len=*),intent(in)::w
		type(SymDimensionName)::tempName
		tempName=w
		 Symequal_name1= Symequal_name2(dimname,tempName)
		return
	end function
	logical function  Symequal_name3(w,dimname)
		type(SymDimensionName),intent(in)::dimname
		CHARACTER(len=*),intent(in)::w
		type(SymDimensionName)::tempName
		tempName=w
		 Symequal_name3= Symequal_name2(dimname,tempName)
		return
	end function	
!|S,rs>,S=-1,0,1,degeneracy: 1,2,1
!in the symmetry base
! |-1,1>,|0,1>,|0,2>,|1,1>
!transfrom to non-symmety index will be
! |1>  , |2>  ,|3>  , |4>
!input s, output imin and imax
! example,input s=0,output imin=2,imax=3
! example,input s=1,output imin=4,imax=4
	subroutine QN2nonSyminde_one(dimen,ith,QN,imin,imax)
		type(SymDimension),intent(in)::dimen
		real,intent(in)::QN
		integer,intent(in)::ith
		integer,intent(out)::imin,imax
		integer::i
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,QN2nonSyminde"
			stop
		end if
		imin=1
		do i=1,size(dimen%dimName(ith)%QN)
			if(abs(dimen%dimName(ith)%QN(i)-QN).gt.zero_error_)then
				imin=imin+dimen%dimName(ith)%degeneracy(i)
			else
				imax=dimen%dimName(ith)%degeneracy(i)
				exit
			end if
		end do
		imax=imax+imin-1
		return
	end subroutine
!|S,rs>,S=-1,0,1,degeneracy: 1,2,1
!in the symmetry base
! |-1,1>,|0,1>,|0,2>,|1,1>
!transfrom to non-symmety index will be
! |1>  , |2>  ,|3>  , |4>
!input s, output imin and imax,thay are array
! example,dimen=[2,2,3],QN=[0.5,0.5,1],degeneracy=[(1,1),(1,1),(1,2,1)]
! input [0.5,-0.5,0] output imin:[1,2,2] imax:[1,2,3]
	subroutine QN2nonSyminde(dimen,QN,imin,imax)
		type(SymDimension),intent(in)::dimen
		real,intent(in)::QN(:)
		integer,intent(out)::imin(:),imax(:)
		integer::i
		do i=1,dimen%Dimsize
			call QN2nonSyminde_one(dimen,i,QN(i),imin(i),imax(i))
		end do
		return
	end subroutine
!U(1) symmetry Tensor to Tensor
!calculate the dimension in Tensor, according to the SymDimension
	subroutine QN2nonSymmaxdim(dimen,outdim)
		type(SymDimension),intent(in)::dimen
		integer,intent(out)::outdim(:)
		integer::i
		do i=1,dimen%Dimsize
			outdim(i)=sum(dimen%dimName(i)%degeneracy)
		end do
		return
	end subroutine
	
!the dimensionis[2,2,3] the quantum number are [(-0.5,0.5),(-0.5,0.5),(-1,0,1)]
!output the index of the dimension,whoes quantum number is [0.5,0.5,0]
!that is [2,2,2]
	subroutine Symindex(outvec,dimen,QN)
		type(SymDimension),intent(in)::dimen
		integer,intent(out)::outvec(:)
		real*4,intent(in)::QN(:)
		integer::i,j
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,SymQN"
			stop
		end if
		outvec=-1
		do i=1,dimen%Dimsize
			do j=1,size(dimen%dimName(i)%QN)
				if(abs(dimen%dimName(i)%QN(j)-QN(i)).lt.zero_error_) then
					outvec(i)=j
					exit
				end if
			end do
			if(outvec(i).eq.-1)then
				write(*,*)"ERROR in Symindex"
				write(*,*)dimen%dimName(i)%QN
				write(*,*)QN(i)
				stop
			end if
		end do
		return
	end subroutine
	integer function Symindexi(dimen,ith,QN)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::ith
		real*4,intent(in)::QN
		integer::j
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,SymQN"
			stop
		end if
		do j=1,size(dimen%dimName(ith)%QN)
			if(abs(dimen%dimName(ith)%QN(j)-QN).lt.zero_error_) then
				Symindexi=j
				return
			end if
		end do
		Symindexi=0
		return
	end function
!the dimensionis[2,2,3] the quantum number are [(-0.5,0.5),(-0.5,0.5),(-1,0,1)]
!the degeneracy are [(1,1),(1,1),(1,2,1)]
!output the degeneracy the dimension [1,1,3],That is [1,1,1]
	subroutine SymDimDegeneracy1(outvec,dimen,inde)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::inde(:)
		integer,intent(out)::outvec(:)
		integer::i
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,SymDimDegeneracy"
			stop
		end if
		do i=1,dimen%Dimsize
			outvec(i)=dimen%dimName(i)%degeneracy(inde(i))
		end do
		return
	end subroutine
!the dimensionis[2,2,3] the quantum number are [(-0.5,0.5),(-0.5,0.5),(-1,0,1)]
!the degeneracy are [(1,1),(1,1),(1,2,1)]
!output the degeneracy the dimension [-0.5,-0.5,1],That is [1,1,1]
	subroutine SymDimDegeneracy2(outvec,dimen,QN)
		type(SymDimension),intent(in)::dimen
		integer,intent(out)::outvec(:)
		real*4,intent(in)::QN(:)
		integer::i,j
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,SymDimDegeneracy"
			stop
		end if
		do i=1,dimen%Dimsize
			do j=1,size(dimen%dimName(i)%QN)
				if(abs(dimen%dimName(i)%QN(j)-QN(i)).lt.zero_error_) then
					outvec(i)=dimen%dimName(i)%degeneracy(j)
					exit
				end if
			end do
		end do
		return
	end subroutine
	subroutine SymDimDegeneracy3(outvec,dimen,ith,QN)
		type(SymDimension),intent(in)::dimen
		integer,intent(out)::outvec(:)
		integer,intent(in)::ith(:)
		real*4,intent(in)::QN(:)
		integer::i,j
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,SymDimDegeneracy"
			stop
		end if
		do i=1,size(ith)
			do j=1,size(dimen%dimName(ith(i))%QN)
				if(abs(dimen%dimName(ith(i))%QN(j)-QN(ith(i))).lt.zero_error_) then
					outvec(i)=dimen%dimName(ith(i))%degeneracy(j)
					exit
				end if
			end do
		end do
		return
	end subroutine
!the dimensionis[2,2,3] the quantum number are [(-0.5,0.5),(-0.5,0.5),(-1,0,1)]
!the degeneracy are [(1,1),(1,1),(1,2,1)]
!output the degeneracy the dimension 3,That is [1,2,1]
	subroutine SymDimDegeneracy4(outvec,dimen,ith)
		type(SymDimension),intent(in)::dimen
		integer,intent(out)::outvec(:)
		integer,intent(in)::ith
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,SymDimDegeneracy"
			stop
		end if
		outvec=dimen%dimName(ith)%degeneracy
		return
	end subroutine
	subroutine symResetdeg1(dimen,ith,jth,deg)
		type(SymDimension),intent(inout)::dimen
		integer,intent(in)::ith,jth,deg
		dimen%dimName(ith)%degeneracy(jth)=deg
		return
	end subroutine
	subroutine symResetdeg2(dimen,ith,deg)
		type(SymDimension),intent(inout)::dimen
		integer,intent(in)::ith,deg
		dimen%dimName(ith)%degeneracy=deg
		return
	end subroutine
	subroutine symResetdeg3(dimen,ith,deg)
		type(SymDimension),intent(inout)::dimen
		integer,intent(in)::ith,deg(:)
		dimen%dimName(ith)%degeneracy=deg
		return
	end subroutine
	
!output the rule of the dimension
	subroutine SymRule1(outvec,dimen)
		type(SymDimension),intent(in)::dimen
		integer,intent(out)::outvec(:)
		integer::i
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,SymRule"
			stop
		end if
		do i=1,dimen%Dimsize
			outvec(i)=dimen%dimName(i)%rule
		end do
		return
	end subroutine
	subroutine SymRule2(outvec,dimen,ith)
		type(SymDimension),intent(in)::dimen
		integer,intent(out)::outvec
		integer,intent(in)::ith
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,SymRule"
			stop
		end if
		outvec=dimen%dimName(ith)%rule
		return
	end subroutine
!the dimensionis[2,2,3] the quantum number are [(-0.5,0.5),(-0.5,0.5),(-1,0,1)]
!output the quantum number of [1,2,2],that is [-0.5,0.5,0]
	subroutine SymQN1(outvec,dimen,inde)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::inde(:)
		real*4,intent(out)::outvec(:)
		integer::i
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,SymQN"
			stop
		end if
		do i=1,dimen%Dimsize
			outvec(i)=dimen%dimName(i)%QN(inde(i))
		end do
		return
	end subroutine
!output the max_QN
	subroutine SymQN2(outvec,dimen)
		type(SymDimension),intent(in)::dimen
		real*4,intent(out)::outvec(:)
		integer::i
		if(.not.Symif_original_dim(dimen))then
			write(*,*)"The SymDimension should be in its original SymDimension,SymQN"
			stop
		end if
		do i=1,dimen%Dimsize
			outvec(i)=dimen%dimName(i)%max_QN
		end do
		return
	end subroutine
!output the max_QN
	real*4 function SymQNum1(dimen,ith)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::ith
		SymQNum1=dimen%dimName(ith)%max_QN
		return
	end function
	real*4 function SymQNum2(dimen,ith,jth)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::ith,jth
		SymQNum2=dimen%dimName(ith)%QN(jth)
		return
	end function
	subroutine reverseDimrule1(dimen)
		type(SymDimension),intent(inout)::dimen
		integer::i
		do i=1,dimen%Dimsize
			dimen%dimName(i)%rule=-1*dimen%dimName(i)%rule
		end do
		return
	end subroutine
	subroutine reverseDimrule2(dimen,ith)
		type(SymDimension),intent(inout)::dimen
		integer,intent(in)::ith
			dimen%dimName(ith)%rule=-1*dimen%dimName(ith)%rule
		return
	end subroutine
!degeneracy(i,:) are the degeneracy of QN(i). QN are quantum number(max quantum number corrding to non-zero term)
	subroutine setDimQN1(dimen,QN,degeneracy,rule)
		type(SymDimension),intent(inout)::dimen
		real*4,intent(in)::QN(:)
		integer,intent(in)::rule(:),degeneracy(:,:)
		integer::i,lenQN
		lenQN=size(QN)
		do i=1,lenQN
			call set_QN1(dimen,i,QN(i),degeneracy(i,1:dimen.i.i),rule(i))
		end do
		return
	end subroutine
	subroutine setDimQN2(dimen,QN,rule)
		type(SymDimension),intent(inout)::dimen
		real*4,intent(in)::QN(:)
		integer,intent(in)::rule(:)
		integer::i,lenQN
		lenQN=size(QN)
		do i=1,lenQN
			call set_QN2(dimen,i,QN(i),dimen.i.i,rule(i))
		end do
		return
	end subroutine
	subroutine setDimQN3(dimen,ith,QN,degeneracy,rule)
		type(SymDimension),intent(inout)::dimen
		real*4,intent(in)::QN
		integer,intent(in)::rule,degeneracy(:),ith
		call set_QN1(dimen,ith,QN,degeneracy(1:dimen.i.ith),rule)
		return
	end subroutine
	
	type(SymDimension) function SymU1Dim1(QN,degeneracy,rule)result(dimen)
		real*4,intent(in)::QN(:)
		integer,intent(in)::rule(:),degeneracy(:,:)
		integer::i,lenQN
		integer,allocatable::maxinde(:)
		lenQN=size(QN)
		allocate(maxinde(lenQN))
		maxinde=2*QN+1
		dimen=maxinde
		do i=1,lenQN
			call set_QN1(dimen,i,QN(i),degeneracy(i,1:dimen.i.i),rule(i))
		end do
		return
	end function
!QN and Deg
!0.5:       1  1
!1  :      1  2  1
!1.5:     1  3  3  1
!2  :    1  4  6  4  1		
!2.5:   1 5 10 10  5  1 
!3  :  1 6 15 20 15 6  1
!3.5: 1 7 21 35 35 21 7 1
!These QNs or N*Deg give a good result on MPS
!input QN,Num
!output a random symTensor with degeneracy deg(i)=Num(i)*C_{N-1}^{i-1},i=0,N-1
!if no Num(:), all the degeneracy will be 1
	type(SymDimension) function SymU1Dim2(QN,rule,Num)result(dimen)
		real*4,intent(in)::QN(:)
		integer,intent(in)::rule(:)
		integer,optional,intent(in)::Num(:)
		integer::i,lenQN
		integer,allocatable::maxinde(:)
		lenQN=size(QN)
		allocate(maxinde(lenQN))
		maxinde=2*QN+1
		dimen=maxinde
		if(present(Num))then
			do i=1,lenQN
				call set_QN2(dimen,i,QN(i),dimen.i.i,rule(i),Num(i))
			end do
		else
			do i=1,lenQN
				call set_QN2(dimen,i,QN(i),dimen.i.i,rule(i))
			end do
		end if
		return
	end function
	type(SymDimension) function SymU1Dim3(QN,degeneracy,rule)result(dimen)
		real*4,intent(in)::QN
		integer,intent(in)::rule,degeneracy(:)
		integer,allocatable::maxinde(:)
		allocate(maxinde(1))
		maxinde=2*QN+1
		dimen=maxinde
		call set_QN1(dimen,1,QN,degeneracy(1:maxinde(1)),rule)
		return
	end function
	subroutine Dimtest()
		type(SymDimensionName)::dimenname,dimenname2
		type(SymDimension)::dimen
		!QN,degeneracy,rule)
		dimen=(/2,3,2/)
		call set_QN1(dimen,1,0.5,(/1,1/),1)
		call set_QN1(dimen,2,1.,(/1,2,1/),1)
		call set_QN1(dimen,3,0.5,(/1,1/),-1)
		call SymDprint2(dimen)
		call SymprintdimenQN(dimen,2)
	end subroutine
		
end module
!****************************************************
!*************** ENF OF SymDimension *************
!****************************************************
















