!************************************************************
!************* START OF Dimension *******************
!************************************************************
!*****************     worning    *********************************
!				1.In fortran,the data store in memory as:A_{1,1}->A_{2,1}->A_{3,1}
!			->...A{M,1}->A_{1,2}->A_{2,2}->..->A_{M,2}->..->A_{M-1,N}->A_{M,N}.
!			But in matlab or C,it is as:A{0,0}->A_{0,1}->A_{0,2}->..A{0,N}
!			->A{1,0}->..->A{M,N-1}->A_{M,N}.The code of DimConstract and Dimdecompose
!			is base on the rule of storing data.
!				2.A is a allocatable array and B is array too,when use A=B before
!			allocating storage space for A,A will get no data.This happen in 
!			fortran90 but do not in fortran77.The interface assignments is to  
!			solve this program.
!				3.The code is fortran90 version.
!******************   note   ********************
!				DimData :store the dimension data
!				boundary: store the boundary of every dimension,for
!			 example,the dimension D1=[2,2,3,4,5],after constract the
!			 first and second index ,it will be D2=[4,3,4,5],then
!			 the data of D1 :
!						boundarysize	=6
!						Dimsize			=5
!						boundary		=[0,1,2,3,4,5]
!						DimData			=[2,2,3,4,5]
!					which means the dimension is [2,2,3,4,5]
!			 the data of D2 :
!						boundarysize	=5
!						Dimsize			=4
!						boundary		=[0,2,3,4,5]
!						DimData			=[2,2,3,4,5]
!				which means the dimension is [4,3,4,5] (or[(2*2),3,2,4] )
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
!			Add a new type DimensionName into type(Dimension)  
!						DimData			=[2,2,3,4,5]
!						DimName		=[1,2,3,4,5]  
!                              [a,a,a,a,a] beause this is the same Tensor,indexname are the same
!    DimName are
!		[1,a],[2,a],[3,a],[4,a],[5,a]
!					which means the dimension is [2,2,3,4,5],the first index call "a.1" second call "a.2",
!		When permute if to 
!												[2,3,2,4,5]
!     The value of indexID will be
!												[1,3,2,4,5]
!		 DimName are
!		[1,a],[3,a],[2,a],[4,a],[5,a]
!		The name will change according to the permutation, but the value will not change. So is nameID.
!		nameID is the same in a dimension. When there are two dimesion. 
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
!	type DimensionData
!		CHARACTER(len=len_of_Name)::TensorName
!		CHARACTER(len=len_of_Name)::dimName
!	end type DimensionData
!
!	type Dimension
!		integer,private :: boundarysize=0
!		integer,private :: Dimsize=0
!		integer,allocatable,private :: boundary(:)
!		type(DimensionData),allocatable,private::dimData(:)
!		logical,private::nameflag=.false.
!	end type Dimension
!			But it cost all lot to modify in the files.
!***************
!
! original dimension means that the dimension do not do any contract 
! Can be fast by modifying Dimpermute , DimConstract, DimDecompose and DimDecomposeAll to subroutine 
!
include '/home/sjdong/Tensor/MPI/function.f90'
module Dimension_typede
	use usefull_function
	implicit none
	include "mpif.h"
	CHARACTER*1,save,private::indexsymbol='.'
	CHARACTER*1,save,private::intNamesymbol='_'!!0_0_1.0_1_2
	integer,private,save::len_of_intName=3
	integer,parameter::len_of_Name=20
	integer,parameter::len_of_intName_in_type_define=4!Cannot define using the value that is not parameter
	type DimensionName
		CHARACTER(len=len_of_Name)::TensorName
		CHARACTER(len=len_of_Name)::dimenName
	end type DimensionName
	type DimensionIntName
		integer::TensorName(len_of_intName_in_type_define)=0
		integer::dimenName(len_of_intName_in_type_define)=0
	end type DimensionIntName
	type Dimension
		integer,private :: boundarysize=0!=size(boundary)
		integer,private :: lenDimData=0!=size(DimData)
		integer,private :: Dimsize=0
		integer,allocatable,private :: boundary(:)
		integer,allocatable,private :: DimData(:)
		type(dimensionName),allocatable,private::DimName(:)
		type(DimensionIntName),allocatable,private::DimIntName(:)
		integer,private::nameflag=0!=0  means no name
		                           !=1  means there are names that is CHARACTER(DimName)
		                           !=2  means there are names that is integer(DimintName)
	end type Dimension
	integer,save::test_flag=0
	logical,private::Dim_check_flag=.false.
	private::getSubDimlen!return the length of the ith dimension
	private::Dimensioniniti
	
	interface Nameinit
		module procedure Nameinit1
		module procedure Nameinit2
		module procedure intNameinit1
		module procedure intNameinit2
		module procedure intNameinit3
		module procedure intNameinit4
	end interface	
	
	!find the index,whose name is w	,output the order of it in the dimension
	!If can not find , output 0
	interface Nameorder
		module procedure Nameorder2!(dimen,character)
		module procedure Nameorder3!(dimen,name)
		module procedure Nameorder4
		module procedure Nameorder5
		module procedure Nameorder6
		module procedure NameorderArray!(dimen,character(:)),output a integer(:)
		module procedure intNameorderArray1!(dimen,integer(:),integer(:)),output a integer(:)
		module procedure intNameorderArray2!(dimen,DimensionIntName(:)),output a integer(:)
	end interface
	
!Find all the DimensionNames, whose indexname is oldname
!change them to newname
!If Cannot find it , do nothing
!If input 'A' 'B', all name such as 'A.1','A.2','A.as' will change to 'B.1','B.2','B.as'
!if input 'A.1' 'B.1',the name 'A.1' will change to 'B.1'
	interface setDimName
		module procedure DimName!(dimen,w),set a name to the dimenison
		module procedure integerDimName!(dimen,int(:)),set a name to the dimenison
		module procedure setDimName1
		module procedure setDimName2
		module procedure setDimName3
		module procedure setDimName4
		module procedure setDimName5
		module procedure setDimName6
	end interface
	interface assignment(=)
		module procedure DimInitialization!dimension=vector
		module procedure DimInitialization2
		module procedure copyDimToVec
		module procedure copyName!name1=name2
		module procedure copyintName!intname1=intname2
		module procedure copyNameArray!name1(:)=name2(:),willnot allocate name1
		module procedure copyintNameArray!intname1(:)=intname2(:),willnot allocate intname1
		module procedure charName!name=character
		module procedure Namechar!character=name
		module procedure intNamechar!character=intname
		module procedure charNameArray!name(:)=character(:)
		module procedure NamecharArray!character(:)=name(:)
		module procedure intNamecharArray!character(:)=intname(:)
	end interface
	interface operator(+)
		module procedure Dimadd
		module procedure DimAdd2
		module procedure DimAdd3
	end interface
	interface operator(.sub.)
		module procedure getSubDim2
		module procedure getSubDim2_name
		module procedure getSubDim2_intname2
	end interface
	!getSubDim_intname(dimen,TensorName(:),DimenName(:))
	!!	return the inde  dimension	,outpout in a vector of the dimenison    
! If do the contract, onedimenison will have more than one value
! [2,3,4,5]	-->contract the 2,3 index -->[2,(3,4),5]!the dimension is 2,12,5
! then getSubDim(Dimen,2,dimenVec)==>dimenVec=[3,4]
!	subroutine getSubDim(Dimen,inde,dimenVec)
	interface operator(.i.)
		module procedure Dim_i
	end interface
	interface operator(.equ.)
		module procedure equal_of_array!If two array of one dimension are equal
		module procedure equal_of_dim
		module procedure  equal_name1
		module procedure  equal_name2
		module procedure  equal_name3
		module procedure  equal_name4
		module procedure  equal_intname
	end interface
!**********************************************************
!	Other Function or Subroutine:
!
!		cleanDimension:clean the data in type(Dimension)
!
!		cleanDimensionName:clean the name in the dimension
!
!		chartoName_log:input a character,output the two element of the 
!					type(DimensionName),example:input asd.ad output 'asd' and
!					 'ad', chartoName_log=.true.  ;  input asd output 'asd' ,
!					 '0'  and chartoName_log=.false.
!
!		copyNameAllocate(outName,inName):outName(:)=inName(:),will allocate outName
!
!		copyintNameAllocate(outName,inName):outName(:)=inName(:),will allocate outName
!
!		printName:print the name in the dimension
!
!		Dimensioniniti_name(boundarysize,Dimsize,boundary,DimData,DimName): initialation a dimension with a name
!
!		outNameTen(dimen,ith):output the TensorName of the ith dimension,as output as CHARACTER, no matter intName of charName
!
!		outIntNameTen(dimen,ith):output the TensorName of the ith dimension,output as integer(:)
!
!		outName(dimen,ith):output the index of the ith dimension,which means output TensorName"."dimenName, no matter intName of charName
!
!		outDimintName:output the Name of the ith dimension,output as integer(2,:)
!
!		getSubDimName:get the inde Name in the dimension,ouput in type(DimensionName)
!
!		getSubDimIntName:get the inde Name in the dimension,ouput in type(DimensionIntName)
!
!		RNDim(dimen) :input dimension [1,1,2,1,3,1,1,4,1],output dimenison [2,3,4],if input [1,1,1,1,1],ouput [1] without name
!
!		outDimData:output dimension%DimData
!
!		copydimension:copy dimension to a array,it will allocate array
!
!		print the type(dimension):
!			Dprint
!			Dprint0
!			Dprint2
!
!		Dimsize:return the Dimsize of type(Dimension)
!
!		outtotalData: product(DimData)
!
!		getSubDim(Dimen,inde,dimenVec)return the inde  dimension	,outpout in a vector of the dimenison   	
!
!		getSubDim_intname(dimen,TensorName(:),DimenName(:))
!			
!		DimConstract(dimen,inde,num):
!				1,2,3,4,..inde,inde+1,...rank-->1,2,3,4,..(inde*inde+1*..inde+num),...rank     	
!			if inde+num>boundarysize,it will constract all the index that larger than inde
!
!		DimDecompose(dimen,inde,midindde):
!				(inde,inde+1,..midindde,midindde+1,...rank)-->(inde,inde+1,..midindde),(midindde+1,...rank)	     		
!			if (inde,inde+1 ..midindde),	midindde larger then next elmemt,do nothing
!			for example:	
!				D=[2,2,(3,4,5),2,(3,4)],	boundary=[0,1,2,5,6,8]
!			DimDecompose(D,3,2)=[2,2,(3,4),5,2,(3,4)],boundary=[0,1,2,4,5,6,8]
!			DimDecompose(D,3,4) will do nothing
!
!		DimDecomposeAll:Decompose all dimension
!
!		Dimpermute:permutation of dimension	
!
!if flag=1,inv[1,2,3,4,...,ith,ith+1,...]->[ith,1,2,3,..,ith-1,ith+1,ith+2,...]
!if flag=0,inv[1,2,3,4,...,ith,ith+1,...]->[2,3,..,ith-1,ith,1,ith+1,ith+2,...]
!	permuteorder(outv,inv,ith,flag)
!
!
!
!
!		MPI function:
!			sent_Dimension(dim1,dim2,ID1,ID2,ierr):send the data of dim1 in ID1 to dim2 in ID2
!			BCAST_Dimension(dim1,ID,ierr):BCAST The dim1 in ID to every cpus
!
!			example
!				integer::ierr,proID,proNum
!				type(Dimension)::T1,T2
!				call mpi_init(ierr)
!				call mpi_comm_rank(mpi_comm_world,proID,ierr)
!				call mpi_comm_size(mpi_comm_world,proNum,ierr ) 
!				call sent_Dimension(T1,T2,0,1,ierr)  !T1 in cpu0 sent to cpu2 and store in T2
!				call BCAST_Dimension(T1,0,ierr) 		 !T1 in cpu 0 send to every cpus,store in T1 in other cpus
!
!**********************************************************	
	contains
   subroutine set_len_of_intName(lenN)
      integer,intent(in)::lenN
      if(lenN.gt.len_of_intName_in_type_define)then 
         write(*,*)"The len_of_intName is larger than the max length"
         stop
       end if
      len_of_intName=lenN
      return
   end subroutine
!*****************************************************
!  input A0 output 'A' , 0
!*******  function or subroutine for name   **************
	subroutine cleanDimensionName(Dimen)
		type(Dimension),intent(inout)::	Dimen
		if(allocated(Dimen%DimName))then
			deallocate(Dimen%DimName)
		end if
		if(allocated(Dimen%DimIntName))then
			deallocate(Dimen%DimIntName)
		end if
		Dimen%nameflag=0
		return
	end subroutine
!  input asd.ad output 'asd' and 'ad', chartoName_log=.true.
!	!  input asd output 'asd' , '0'  and chartoName_log=.false.
	logical function chartoName_log(TensorName,DimenName,w_)
		character(len=*),intent(in)::w_
		character(len=len(trim(adjustl(w_))))::w
		character(len=len_of_Name),intent(out)::DimenName,TensorName
		integer::lenw,lenw2
		w=trim(adjustl(w_))
		lenw=index(w,indexsymbol)
		lenw2=len(trim(w))
		if(lenw.eq.0)then!input asd
			chartoName_log=.false.
			TensorName=w
			DimenName='0'
			return
		end if
		chartoName_log=.true.
		!input 'A2.1'(if indexsymbol='.') 
		TensorName=w(1:lenw-1)
		DimenName=w(lenw+1:lenw2)
		return
	end function
!*******  function or subroutine for name   **************
	type(DimensionName) function Nameinit1(TensorName,DimenName)
		character(len=*),intent(in)::TensorName
		character(len=*),intent(in)::DimenName
		Nameinit1%TensorName=TensorName
		Nameinit1%DimenName=DimenName
		return
	end function
	type(DimensionName) function Nameinit2(TensorName,DimenName)
		character(len=*),intent(in)::TensorName
		integer,intent(in)::DimenName
		Nameinit2%TensorName=TensorName
		Nameinit2%DimenName=DimenName
		return
	end function
	type(DimensionIntName) function intNameinit1(TensorName,DimenName)
		integer,intent(in)::TensorName(:),DimenName(:)
		integer::lenname
		lenname=size(TensorName)
		intNameinit1%TensorName(1:lenname)=TensorName
		intNameinit1%DimenName(1:lenname)=DimenName
		return
	end function 
	type(DimensionIntName) function intNameinit2(intDimname)
		integer,intent(in)::intDimname(:,:)
		integer::lenname
		lenname=size(intDimname,2)
		intNameinit2%TensorName(1:lenname)=intDimname(1,:)
		intNameinit2%DimenName(1:lenname)=intDimname(2,:)
		return
	end function 
	type(DimensionIntName) function intNameinit3(TensorName,DimenName)
		integer,intent(in)::TensorName(:),DimenName
		integer::lenname
		lenname=size(TensorName)
		intNameinit3%TensorName(1:lenname)=TensorName
		intNameinit3%DimenName(1)=DimenName
		return
	end function 
	type(DimensionIntName) function intNameinit4(TensorName,DimenName)
		integer,intent(in)::TensorName,DimenName
		intNameinit4%TensorName(1)=TensorName
		intNameinit4%DimenName(1)=DimenName
		return
	end function 
		
!*******  function or subroutine for name   **************
	subroutine copyName(Nameout,Namein)
		type(DimensionName),intent(inout)::Nameout
		type(DimensionName),intent(in)::Namein
		Nameout%TensorName=Namein%TensorName
		Nameout%DimenName=Namein%DimenName
		return
	end subroutine		
	subroutine copyintName(Nameout,Namein)!the length of TensorName in Nameout or Namein are the same
		type(DimensionIntName),intent(inout)::Nameout
		type(DimensionIntName),intent(in)::Namein
		Nameout%TensorName=Namein%TensorName
		Nameout%DimenName=Namein%DimenName
		return
	end subroutine		
!*******  function or subroutine for name   **************
	subroutine copyNameArray(Nameout,Namein)
		type(DimensionName),intent(inout)::Nameout(:)
		type(DimensionName),intent(in)::Namein(:)
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
	subroutine copyintNameArray(Nameout,Namein)
		type(DimensionIntName),intent(inout)::Nameout(:)
		type(DimensionIntName),intent(in)::Namein(:)
		integer::Lenout,lenIN,i
		lenIN=size(Namein)
		if(size(Nameout).lt.lenIN)then
			write(*,*)"ERROR in assignment of two intName array "
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
!w='A','A1','A1.1'
	subroutine charName(outName,w)
		type(DimensionName),intent(inout)::outName
		character(len=*),intent(in)::w
		logical::flag
		flag=chartoName_log(outName%TensorName,outName%DimenName,w)
		return
	end subroutine
	subroutine Namechar(w,inName)
		type(DimensionName),intent(in)::inName
		character(len=*),intent(inout)::w
		w=inName%TensorName+indexsymbol+inName%DimenName
		return
	end subroutine
!example:inName is( [0,0,1][0,1,2])
!output character will be "0_0_1.0_1_2"
	subroutine intNamechar(w,inName)
		type(DimensionIntName),intent(in)::inName
		character(len=*),intent(inout)::w
		integer::i
		w=""
		do i=1,len_of_intName
			w=w+inName%TensorName(i)
			if(i.ne.len_of_intName) w=w+intNamesymbol
		end do
		w=w+indexsymbol
		do i=1,len_of_intName
			w=w+inName%DimenName(i)
			if(i.ne.len_of_intName) w=w+intNamesymbol
		end do
		return
	end subroutine
	
!	input array ,for example (/'A1.12','A1.2 ','A1.3 '/)
!	output a array of DimensionName
!	the array input should be the same length
!	'A1.12' length is 5, 'A1.2 ' is 5 too
!name(:)=w(:)
	 subroutine charNameArray(R,w)
	 	type(DimensionName),intent(inout)::R(:)
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
!character(:)=DimensionName(:)
	subroutine NamecharArray(w,inName)
	 	type(DimensionName),intent(in)::inName(:)
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
	subroutine intNamecharArray(w,inName)
	 	type(DimensionIntName),intent(in)::inName(:)
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
	subroutine copyNameAllocate(outName,inName)
		type(DimensionName),allocatable,intent(inout)::outName(:)
		type(DimensionName),intent(in)::inName(:)
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
	subroutine copyintNameAllocate(outName,inName)
		type(DimensionIntName),allocatable,intent(inout)::outName(:)
		type(DimensionIntName),intent(in)::inName(:)
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
	
!*******  function or subroutine for name   **************
!set a name to the diemsnion
!w may be something like 'A' or 'A1'
	subroutine DimName(dimen,w)
		type(Dimension),intent(inout)::dimen
		character(len=*),intent(in)::w
		character(len=len_of_Name)::TensorName,DimenName
		integer::lenD,i
		logical::flag
		flag=chartoName_log(TensorName,DimenName,w)
		if(flag)then
			write(*,*)"ERROR in Set the name to the Tensor"
			write(*,*)"input:",w
			write(*,*)"could not contains the indexsymbol:",indexsymbol
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"The dimension should be in its original dimension"
			stop
		end if
		if(dimen%nameflag.eq.2) then
			write(*,*)"The Dimension have the name of integer"
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
		if(dimen%nameflag.eq.1)then!There are name in the dimension,rename
			do i=1,lenD
				dimen%dimname(i)%TensorName=TensorName
			end do
		else
			do i=1,lenD
				dimen%DimName(i)=Nameinit(TensorName,i)
			end do
			dimen%nameflag=1
		end if
		return
	end subroutine
	subroutine integerDimName(dimen,TensorName)!set TensorName
		type(Dimension),intent(inout)::dimen
		integer,intent(in)::TensorName(:)
		integer::lenD,i
		if(.not.if_original_dim(dimen))then
			write(*,*)"The dimension should be in its original dimension"
			stop
		end if
		if(dimen%nameflag.eq.1) then
			write(*,*)"The Dimension have the name of character"
			stop
		end if
		lenD=size(dimen%DimData)
		if(allocated(dimen%DimintName)) then
			if(lenD.ne.size(dimen%DimintName)) then
				deallocate(dimen%DimintName)
				allocate(dimen%DimintName(lenD))
			end if
		else
			allocate(dimen%DimintName(lenD))
		end if
		if(dimen%nameflag.eq.2)then!There are name in the dimension,rename
			do i=1,lenD
				dimen%DimintName(i)%TensorName=TensorName
			end do
		else
			do i=1,lenD
				dimen%DimintName(i)=Nameinit(TensorName,i)
			end do
			dimen%nameflag=2
		end if
		return
	end subroutine

	
	subroutine printName(na)
		type(DimensionName),intent(in)::na
		character*1000::w
		w=na
		write(*,*)trim(adjustl(w))
		return
	end subroutine
	subroutine printIntName(na)
		type(DimensionName),intent(in)::na
		character*1000::w
		w=na
		write(*,*)trim(adjustl(w))
		return
	end subroutine
			
		
		
!***************   cleanDimension   *****************
	subroutine cleanDimension(Dimen)
		type(Dimension),intent(inout)::	Dimen
		Dimen%boundarysize=0
		Dimen%Dimsize=0
		Dimen%lenDimData=0
		if(allocated(Dimen%boundary)) then
			deallocate(Dimen%boundary)
		end if
		if(allocated(Dimen%DimData)) then
			deallocate(Dimen%DimData)
		end if
		call cleanDimensionName(Dimen)
		return
	end subroutine
	
	subroutine outDimData(DimData,dimen)
		type(Dimension),intent(in)::dimen
		integer,intent(inout)::DimData(:)
		DimData=dimen%DimData
		return
	end subroutine
!****************print the information of dimension ****************************************	
	subroutine Dprint1(Dimen)
		type(Dimension),intent(in) ::Dimen
		integer,allocatable :: dimenVec(:)
		CHARACTER*5000,allocatable::w(:)
		integer::i
		write(*,*) "***   START   ***"
		write(*,*) Dimen%DimData
		call copydimension(dimenVec,Dimen)
		if(Dimen%nameflag.eq.1)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		if(Dimen%nameflag.eq.2)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimintName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		
		write(*,*) "			--				"
		write(*,*)	dimenVec
		write(*,*) "***   END   ***"
	end subroutine
				
	subroutine Dprint0(Dimen)
		type(Dimension),intent(in) ::Dimen
		integer,allocatable :: dimenVec(:)
		CHARACTER*5000,allocatable::w(:)
		integer::i
		write(*,*) "***   Dimension Data    ***"
		write(*,*) Dimen%DimData
		call copydimension(dimenVec,Dimen)
		write(*,*)	dimenVec
		write(*,*) "***   Dimension END   ***"
		if(Dimen%nameflag.eq.1)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		if(Dimen%nameflag.eq.2)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimintName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
	end subroutine
	subroutine Dprint(Dimen)
		type(Dimension),intent(in) ::Dimen
		integer,allocatable :: dimenVec(:)
		CHARACTER*5000,allocatable::w(:)
		integer::i
		write(*,*) "***   Dimension Data    ***"
		call copydimension(dimenVec,Dimen)
		write(*,*)	dimenVec
		write(*,*) "***   Dimension END   ***"
		if(Dimen%nameflag.eq.1)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		if(Dimen%nameflag.eq.2)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimintName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		write(*,*) " "
	end subroutine
				
	subroutine Dprint2(Dimen)
		type(Dimension),intent(in) ::Dimen
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
		if(Dimen%nameflag.eq.1)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		if(Dimen%nameflag.eq.2)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimintName
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		call copydimension(dimenVec,Dimen)
		write(*,*) "			--				"
		write(*,*)	dimenVec
		write(*,*) "***   END   ***"
		write(*,*) " "
	end subroutine
	
	type(Dimension) function Dimensioniniti(boundarysize,Dimsize,boundary,DimData)
		integer,intent(in) :: boundarysize
		integer,intent(in) :: Dimsize
		integer,intent(in) :: boundary(:)
		integer,intent(in) :: DimData(:)
		allocate(Dimensioniniti%boundary(size(boundary,1)))
		allocate(Dimensioniniti%DimData(size(DimData,1)))
		Dimensioniniti%boundary=boundary
		Dimensioniniti%DimData=DimData
		Dimensioniniti%boundarysize=boundarysize
		Dimensioniniti%Dimsize=Dimsize
		Dimensioniniti%lenDimData=size(DimData,1)
		Dimensioniniti%Nameflag=0
		return
	end function
	subroutine DimInitialization(Dimen,DimData)
		type(Dimension),intent(inout) ::Dimen
		integer,intent(in) :: DimData(:)
		integer::i
		Dimen%Dimsize=size(DimData)
		Dimen%lenDimData=Dimen%Dimsize
		Dimen%boundarysize=Dimen%Dimsize+1
		call allocateCheck(Dimen%boundary,Dimen%boundarysize)
		do i=1,Dimen%boundarysize
			Dimen%boundary(i)=i-1
		end do
		call allocateCheck(Dimen%DimData,Dimen%lenDimData)
		Dimen%DimData(1:Dimen%lenDimData)=DimData
		Dimen%nameflag=0
		return
	end subroutine
	subroutine DimInitialization2(Dimen,Dimen2)
		type(Dimension),intent(inout) ::Dimen
		type(Dimension),intent(in) ::Dimen2
		integer::lenDim,lenBoun,i
		lenDim=Dimen2%lenDimData
		lenBoun=Dimen2%boundarysize
		
		Dimen%Dimsize=Dimen2%Dimsize
		Dimen%boundarysize=lenBoun
		Dimen%lenDimData=lenDim
		call allocateCheck(Dimen%DimData , lenDim)
		call allocateCheck(Dimen%boundary , lenBoun)
		Dimen%boundary(1:lenBoun)=Dimen2%boundary(1:lenBoun)
		Dimen%DimData(1:lenDim)=Dimen2%DimData(1:lenDim)
		Dimen%nameflag=Dimen2%nameflag
		if(Dimen%nameflag.eq.0) return
		if(Dimen2%nameflag.eq.1)then
			if(allocated(dimen%dimName))then
				if(size(dimen%dimName).lt.lenDim) then
					deallocate(dimen%dimName)
					allocate(dimen%dimName(lenDim))
				end if
			else
					allocate(dimen%dimName(lenDim))
			end if
			do i=1,lenDim
				dimen%dimName(i)=dimen2%dimName(i)
			end do
		end if
		if(Dimen2%nameflag.eq.2)then
			if(allocated(dimen%dimintName))then
				if(size(dimen%dimintName).lt.lenDim) then
					deallocate(dimen%dimintName)
					allocate(dimen%dimintName(lenDim))
				end if
			else
					allocate(dimen%dimintName(lenDim))
			end if
			do i=1,lenDim
				dimen%dimintName(i)=dimen2%dimintName(i)
			end do
		end if
		return
	end subroutine
	integer	function DimSize(Dimen)
		type(Dimension),intent(in) :: Dimen
		DimSize=Dimen%dimSize
		return
	end function
!*******  function or subroutine for name   **************		
!output the TensorName of the ith dimension
	CHARACTER(len=len_of_Name) function outNameTen(dimen,ith)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::ith
		integer::i
		if(Dimen%nameflag.eq.0)then
			write(*,*)"There is no CHARACTER name in the dimension"
			stop
		end if
		if(ith.gt.Dimen%lenDimData)then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,Dimen%lenDimData
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"NameID is use in original dimension"
			stop
		end if
		if(Dimen%nameflag.eq.1)then
			outNameTen=Dimen%DimName(ith)%TensorName
			return
		end if
		if(Dimen%nameflag.eq.2)then
			outNameTen=''
			do i=1,len_of_intName
				outNameTen=outNameTen+Dimen%DimintName(ith)%TensorName(i)
				if(i.ne.len_of_intName)outNameTen=outNameTen+intNamesymbol
			end do
			return
		end if
		write(*,*)"ERROR of nameFlag"
		stop
	end function
	function outIntNameTen(dimen,ith)
		integer::outIntNameTen(len_of_intName)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no integer name in the dimension"
			stop
		end if
		if(ith.gt.size(Dimen%DimintName))then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,size(Dimen%DimintName)
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"NameID is use in original dimension"
			stop
		end if
		outIntNameTen=Dimen%DimintName(ith)%TensorName(1:len_of_intName)
		return
	end function
!*******  function or subroutine for name   **************		
!output the index of the ith dimension
	CHARACTER(len=len_of_Name+len_of_Name) function outName(dimen,ith)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(Dimen%nameflag.eq.0)then
			write(*,*)"There is no name in the dimension"
			stop
		end if
		if(ith.gt.size(Dimen%DimName))then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,size(Dimen%DimName)
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"indexID is use in original dimension"
			stop
		end if
		if(Dimen%nameflag.eq.1)then
			outName=Dimen%DimName(ith)
			return
		end if
		if(Dimen%nameflag.eq.2)then
			outName=Dimen%DimintName(ith)
			return
		end if
		write(*,*)"ERROR of nameFlag"
		stop
	end function
	function outDimintName(dimen,ith)
		integer::outDimintName(2,len_of_intName)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no integer name in the dimension"
			stop
		end if
		if(ith.gt.size(Dimen%DimintName))then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,size(Dimen%DimintName)
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"NameID is use in original dimension"
			stop
		end if
		outDimintName(1,:)=Dimen%DimintName(ith)%TensorName(1:len_of_intName)
		outDimintName(2,:)=Dimen%DimintName(ith)%DimenName(1:len_of_intName)
		return
	end function
		
		
!*******  function or subroutine for name   **************		
!find the index,whose name is w	,output the order of it in the dimension
!If can not find , output 0
	integer function Nameorder2(dimen,w)
		type(Dimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w
		integer::i
		type(DimensionName)::nam
		if(Dimen%nameflag.ne.1)then
			write(*,*)"There is no char name in the dimension,Nameorder2"
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"Nameorder2 is use in original dimension"
			stop
		end if
		nam=w
		do i=1,dimen%lenDimData
			if(dimen%dimname(i).equ.nam)then
				Nameorder2=i
				return
			end if
		end do
		Nameorder2=0
		return
	end function
	integer function Nameorder3(dimen,dimname)
		type(Dimension),intent(in) :: Dimen
		integer::i
		type(DimensionName),intent(in)::dimname
		if(Dimen%nameflag.ne.1)then
			write(*,*)"There is no char name in the dimension,Nameorder3"
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"Nameorder3 is use in original dimension"
			stop
		end if
		do i=1,dimen%lenDimData
			if(dimen%dimname(i).equ.dimname)then
				Nameorder3=i
				return
			end if
		end do
		Nameorder3=0
		return
	end function
	integer function Nameorder4(dimen,TensorName,dimenName)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::TensorName(len_of_intName),dimenName(len_of_intName)
		integer::i
		type(DimensionIntName)::testName 
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,Nameorder4"
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"Nameorder4 is use in original dimension"
			stop
		end if
		testName=Nameinit(TensorName,DimenName)
		do i=1,dimen%lenDimData
			if(dimen%dimintname(i).equ.testName)then
				Nameorder4=i
				return
			end if
		end do
		Nameorder4=0
		return
	end function
	integer function Nameorder5(dimen,TensoDimrName)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::TensoDimrName(2,len_of_intName)
		integer::i
		type(DimensionIntName)::testName 
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,Nameorder5"
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"Nameorder5 is use in original dimension"
			stop
		end if
		testName=Nameinit(TensoDimrName)
		do i=1,dimen%lenDimData
			if(dimen%dimintname(i).equ.testName)then
				Nameorder5=i
				return
			end if
		end do
		Nameorder5=0
		return
	end function
	
	integer function Nameorder6(dimen,TensoDimName)
		type(Dimension),intent(in) :: Dimen
		type(DimensionIntName),intent(in)::TensoDimName
		integer::i
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,Nameorder5"
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"Nameorder5 is use in original dimension"
			stop
		end if
		do i=1,dimen%lenDimData
			if(dimen%dimintname(i).equ.TensoDimName)then
				Nameorder6=i
				return
			end if
		end do
		Nameorder6=0
		return
	end function
	
	function NameorderArray(dimen,w)
		type(Dimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w(:)
		integer::NameorderArray(size(w))
		integer::i,lenn
		if(Dimen%nameflag.ne.1)then
			write(*,*)"There is no char name in the dimension,NameorderArray"
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"NameorderArray is use in original dimension"
			stop
		end if
		lenn=size(w)
		if(lenn.gt.dimen%lenDimData)then
			write(*,*)"ERROR i NameorderArray"
			write(*,*)lenn,dimen%Dimsize
			stop
		end if
		do i=1,lenn
			NameorderArray(i)=Nameorder(dimen,w(i))
		end do
		return
	end function
	function intNameorderArray1(dimen,TensorName,dimenName)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::TensorName(:,:),dimenName(:,:)
		integer::intNameorderArray1(size(TensorName,1))
		integer::i,lenn
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,NameorderArray"
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"intNameorderArray is use in original dimension"
			stop
		end if
		lenn=size(TensorName,1)
		if(lenn.gt.dimen%lenDimData)then
			write(*,*)"ERROR i intNameorderArray"
			write(*,*)lenn,dimen%Dimsize
			stop
		end if
		do i=1,lenn
			intNameorderArray1(i)=Nameorder(dimen,TensorName(i,:),dimenName(i,:))
		end do
		return
	end function
	function intNameorderArray2(dimen,TensorintName)
		type(Dimension),intent(in) :: Dimen
		type(DimensionIntName),intent(in)::TensorintName(:)
		integer::intNameorderArray2(size(TensorintName))
		integer::i,lenn
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,NameorderArray"
			stop
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"intNameorderArray is use in original dimension"
			stop
		end if
		lenn=size(TensorintName)
		if(lenn.gt.dimen%lenDimData)then
			write(*,*)"ERROR i intNameorderArray"
			write(*,*)lenn,dimen%Dimsize
			stop
		end if
		do i=1,lenn
			intNameorderArray2(i)=Nameorder6(dimen,TensorintName(i))
		end do
		return
	end function
!*******  function or subroutine for name   **************			
	logical function equal_name2(dimname,dimname2)
		type(DimensionName),intent(in)::dimname,dimname2
		if(trim(adjustl(dimname%TensorName)).ne.trim(adjustl(dimname2%TensorName)))then
			equal_name2=.false.
			return
		end if
		if(trim(adjustl(dimname%dimenName)).ne.trim(adjustl(dimname2%dimenName)))then
			equal_name2=.false.
			return
		end if
		equal_name2=.true.
		return
	end function
	logical function equal_name1(dimname,w)
		type(DimensionName),intent(in)::dimname
		CHARACTER(len=*),intent(in)::w
		type(DimensionName)::tempName
		tempName=w
		equal_name1=equal_name2(dimname,tempName)
		return
	end function
	logical function equal_name3(w,dimname)
		type(DimensionName),intent(in)::dimname
		CHARACTER(len=*),intent(in)::w
		type(DimensionName)::tempName
		tempName=w
		equal_name3=equal_name2(dimname,tempName)
		return
	end function	
	logical function equal_name4(w1,w2)
		CHARACTER(len=*),intent(in)::w1
		CHARACTER(len=*),intent(in)::w2
		equal_name4=trim(adjustl(w1)).eq.trim(adjustl(w2))
		return
	end function	
	logical function equal_intname(name1,name2)
		type(DimensionIntName),intent(in)::name1,name2
		logical::lg1,lg2
		lg1=equal_of_array(name1%TensorName(:len_of_intName),name2%TensorName(:len_of_intName))
		lg2=equal_of_array(name1%dimenName(:len_of_intName),name2%dimenName(:len_of_intName))
		equal_intname=lg1.and.lg2
		return
	end function
	
!*******  function or subroutine for name   **************		
!Find all the DimensionNames, whose indexname is oldname
!change them to newname
!If Cannot find it , do nothing
	subroutine setDimName1(dimen,oldname,newname)
		type(Dimension),intent(inout) :: dimen
		CHARACTER(len=*),intent(in)::oldname,newname
		integer::i
		CHARACTER(len=len_of_Name)::oldTensorName,olddimenName
		CHARACTER(len=len_of_Name)::newTensorName,newdimenName
		logical::oldfullname,newfullname
		if(Dimen%nameflag.ne.1)then
			write(*,*)"There is no char name in the dimension,resetindexname"
			stop
		end if
		oldfullname=chartoName_log(oldTensorName,olddimenName,oldname)
		!(TensorName,DimenName,w_)
		newfullname=chartoName_log(newTensorName,newdimenName,newname)
		if(oldfullname.neqv.newfullname)then
			write(*,*)"ERROR in resetIndexname1"
			stop
		end if
		if(newfullname)then
			do i=1,size(dimen%DimName)
				if(dimen%dimname(i)%TensorName.equ.oldTensorName)then
					if(dimen%dimname(i)%dimenName.equ.olddimenName)then
						dimen%dimname(i)%TensorName=newTensorName
						dimen%dimname(i)%dimenName=newdimenName
					end if
				end if
			end do
		else
			do i=1,size(dimen%DimName)
				if(dimen%dimname(i)%TensorName.equ.oldTensorName)then
						dimen%dimname(i)%TensorName=newTensorName
				end if
			end do
		end if
		return
	end subroutine
!*******  function or subroutine for name   **************		
!set the NameID and Indexname of the ith index
	subroutine setDimName2(dimen,ith,newname)
		type(Dimension),intent(inout) :: dimen
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
		if(.not.if_original_dim(dimen))then
			write(*,*)"resetIndexNameID is use in original dimension"
			stop
		end if
		fullname=chartoName_log(newTensorName,newdimenName,newname)
		if(Dimen%nameflag.eq.0)then!There is no name in dimension
			lenName=size(dimen%DimData)
			allocate(dimen%dimname(lenName))
			do i=1,lenName
				dimen%dimname(i)%TensorName=newTensorName
				dimen%dimname(i)%dimenName=i
			end do
			if(fullname)then
				dimen%dimname(ith)%dimenName=newdimenName
			end if
			Dimen%nameflag=1
			return
		end if
		
		if(fullname)then
			dimen%dimname(ith)%TensorName=newTensorName
			dimen%dimname(ith)%dimenName=newdimenName
		else
			dimen%dimname(ith)%TensorName=newTensorName
		end if
		return
	end subroutine
!Find all the DimensionNames, whose indexname is oldname
!change them to newname
!If Cannot find it , do nothing
	subroutine setDimName3(dimen,oldnameTen,oldnamedim,newnameTen,newnameDim)
		type(Dimension),intent(inout) :: dimen
		integer,intent(in)::oldnameTen(:),oldnamedim(:)
		integer,intent(in)::newnameTen(:),newnameDim(:)
		integer::i
		type(DimensionIntName)::oldname,newName
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,resetindexname"
			stop
		end if
		oldname=Nameinit(oldnameTen,oldnamedim)
		newname=Nameinit(newnameTen,newnameDim)
		do i=1,dimen%lenDimData
			if(dimen%dimintname(i).equ.oldname)then
				dimen%dimintname(i)=newName
			end if
		end do
		return
	end subroutine
	subroutine setDimName4(dimen,oldnameTen,newnameTen)
		type(Dimension),intent(inout) :: dimen
		integer,intent(in)::oldnameTen(:)
		integer,intent(in)::newnameTen(:)
		integer::i
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,resetindexname"
			stop
		end if
		do i=1,dimen%lenDimData
			if(dimen%dimintname(i)%TensorName(:len_of_intName).equ.oldnameTen)then
				dimen%dimintname(i)%TensorName(:len_of_intName)=newnameTen
			end if
		end do
		return
	end subroutine
	subroutine setDimName5(dimen,ith,newnameTen)
		type(Dimension),intent(inout) :: dimen
		integer,intent(in)::ith
		integer,intent(in)::newnameTen(:)
		integer::i
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,resetindexname"
			stop
		end if
		dimen%dimintname(ith)%TensorName(:len_of_intName)=newnameTen
		return
	end subroutine
	subroutine setDimName6(dimen,ith,newnameTen,newnameDim)
		type(Dimension),intent(inout) :: dimen
		integer,intent(in)::ith
		integer,intent(in)::newnameTen(:),newnameDim(:)
		integer::i
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,resetindexname"
			stop
		end if
		dimen%dimintname(ith)%TensorName(:len_of_intName)=newnameTen
		dimen%dimintname(ith)%DimenName(:len_of_intName)=newnameDim
		return
	end subroutine
	
	integer function outtotalData(Dimen)
		type(Dimension),intent(in) :: Dimen
		outtotalData=product(Dimen%DimData)
		return
	end function				
!	return the inde dimension	,outpout in a integer
	integer function Dim_i(Dimen,inde)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		integer :: i,D1
		if(inde.gt.dimen%dimSize) then 
			write(*,*) "ERROR in Dim_i"
			write(*,*)"stop"
			stop
		end if
		D1=1
		do i=Dimen%boundary(inde) + 1,Dimen%boundary(inde+1)
			D1=D1*Dimen%DimData(i)
		end do
		Dim_i=D1
	end function 
!	return  all the  dimension	,outpout in a vector	
	subroutine copydimension(dimenVec,Dimen)
		integer,allocatable,intent(inout) :: dimenVec(:)
		type(Dimension),intent(in) :: Dimen
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
			dimenVec(i)=Dim_i(Dimen,i)
		end do
	end 	subroutine		
	subroutine copyDimToVec(dimenVec,Dimen)
		integer,intent(inout) :: dimenVec(:)
		type(Dimension),intent(in) :: Dimen
		integer :: i
		if(size(dimenVec).lt.Dimen%Dimsize)then
			write(*,*)"ERROR in assignment dimension to array "
			write(*,*)"array(:)=dimension(:),size(array)<size(dimension)"
			write(*,*)size(dimenVec),Dimen%Dimsize
			stop
		end if
		do i=1,Dimen%Dimsize
			dimenVec(i)=Dim_i(Dimen,i)
		end do
	end 	subroutine

!	return the inde  dimension	,outpout in a vector of the dimenison    
! If do the contract, onedimenison will have more than one value
! [2,3,4,5]	-->contract the 2,3 index -->[2,(3,4),5]!the dimension is 2,12,5
! then getSubDim(Dimen,2,dimenVec)==>dimenVec=[3,4]
	subroutine getSubDim(Dimen,inde,dimenVec)
		type(Dimension),intent(in) :: Dimen
		integer,allocatable,intent(inout) :: dimenVec(:)
		integer,intent(in) :: inde
		integer::i,j,Dlen
		Dlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
		if(inde.gt.dimen%dimSize) then 
			write(*,*) "ERROR in getSubDim"
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
	subroutine getSubDimName(Dimen,inde,dimenName)
		type(Dimension),intent(in) :: Dimen
		type(dimensionName),allocatable,intent(inout) :: dimenName(:)
		integer,intent(in) :: inde
		integer::i,j,Dlen
		Dlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
		if(inde.gt.dimen%dimSize) then 
			write(*,*) "ERROR in getSubDim"
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
	subroutine getSubDimIntName(Dimen,inde,dimenintName)
		type(Dimension),intent(in) :: Dimen
		type(dimensionIntName),allocatable,intent(inout) :: dimenintName(:)
		integer,intent(in) :: inde
		integer::i,j,Dlen
		Dlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
		if(inde.gt.dimen%dimSize) then 
			write(*,*) "ERROR in getSubDim"
			return
		end if
		if(allocated(dimenintName)) then
			if(Dlen.ne.size(dimenintName)) then 
				deallocate(dimenintName)
				allocate(dimenintName(Dlen))
			end if
		else
			allocate(dimenintName(Dlen))
		end if
		j=1
		do i=Dimen%boundary(inde) + 1,Dimen%boundary(inde+1)
			dimenintName(j)=Dimen%DimIntName(i)
			j=j+1
		end do
	end subroutine
!input dimension [1,1,2,1,3,1,1,4,1]
!output dimenison [2,3,4]
! if input [1,1,1,1,1]
!  ouput [1] without name
	type(Dimension) function RNDim(dimen) 
		type(Dimension),intent(inout) :: Dimen
		integer::i,lenD,lenNewD
		integer,allocatable::Dimindex(:)
		if(.not.if_original_dim(dimen)) then
			write(*,*)"ERROR IN RNDim,in Diemnsion.f90"
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
			RNDim=(/1/)
			return
		end if
			
		allocate(RNDim%DimData(lenNewD))
		do i=1,lenNewD
			RNDim%DimData(i)=dimen%DimData(Dimindex(i))
		end do
		RNDim%Dimsize=lenNewD
		RNDim%boundarysize=RNDim%Dimsize+1
		allocate(RNDim%boundary(RNDim%boundarysize))
		do i=1,RNDim%boundarysize
			RNDim%boundary(i)=i-1
		end do
		RNDim%NameFlag=dimen%NameFlag
		if(dimen%NameFlag.eq.1)then
			allocate(RNDim%DimName(lenNewD))
			do i=1,lenNewD
				RNDim%DimName(i)=dimen%DimName(Dimindex(i))
			end do
		end if
		if(dimen%NameFlag.eq.2)then
			allocate(RNDim%DimintName(lenNewD))
			do i=1,lenNewD
				RNDim%DimintName(i)=dimen%DimintName(Dimindex(i))
			end do
		end if
		return
	end function
!*******  function or subroutine for name   **************	
!	return the inde  dimension	,outpout in a type(dimenison)
	type(Dimension) function  getSubDim2(Dimen,inde)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		integer::boundarysize,boundary(2),Dimsize
		integer,allocatable :: Dimdata(:)
		integer::i,j,Dlen
		call getSubDim(Dimen,inde,Dimdata)
		boundary(1)=0
		boundary(2)=size(Dimdata)
		boundarysize=2
		getSubDim2=Dimensioniniti(boundarysize,1,boundary,DimData)
		getSubDim2%NameFlag=Dimen%NameFlag
		if(Dimen%nameFlag.eq.1)then
			call getSubDimName(Dimen,inde,getSubDim2%dimName)
		end if
		if(Dimen%nameFlag.eq.2)then
			call getSubDimIntName(Dimen,inde,getSubDim2%dimIntName)
		end if
		return
	end function
	type(Dimension) function  getSubDim2_name(Dimen,w)
		type(Dimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w
		integer::inde
		inde=Nameorder(Dimen,w)
		getSubDim2_name=getSubDim2(Dimen,inde)
		return
	end function
	type(Dimension) function  getSubDim_intname(Dimen,TensorName,DimenName)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::TensorName(:)
		integer,intent(in)::DimenName(:)
		integer::inde
		inde=Nameorder(Dimen,TensorName,DimenName)
		getSubDim_intname=getSubDim2(Dimen,inde)
		return
	end function
	type(Dimension) function  getSubDim2_intname2(Dimen,TensorDimName)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::TensorDimName(:,:)
		integer::inde
		inde=Nameorder(Dimen,TensorDimName)
		getSubDim2_intname2=getSubDim2(Dimen,inde)
		return
	end function
!*******  function or subroutine for name   **************	
! The type of input are both dimension	
!If dimen have a name but Dimen2 do not or Dimen2 have but not dimen
!then the name for the one with no name will be ['0',0,i]
	type(Dimension) function  DimAdd(Dimen,Dimen2)
		type(Dimension),intent(in) :: Dimen,Dimen2
		integer::i,j,l1,l2
		DimAdd%boundarysize=Dimen%boundarysize+Dimen2%boundarysize-1
		DimAdd%lendimData=Dimen%lendimData+Dimen2%lendimData
		allocate(DimAdd%boundary(DimAdd%boundarysize))
		do i=1,Dimen%boundarysize
			DimAdd%boundary(i)=Dimen%boundary(i)
		end do
		l1=Dimen%boundary(Dimen%boundarysize)
		do i=1,Dimen2%boundarysize-1
			DimAdd%boundary(i+Dimen%boundarysize)=Dimen2%boundary(i+1)+l1
		end do
		l1=size(Dimen%DimData)
		l2=size(Dimen2%DimData)
		DimAdd%lenDimData=l1+l2
		allocate(DimAdd%Dimdata(l1+l2))
		DimAdd%Dimdata(:l1)=Dimen%DimData
		DimAdd%Dimdata(l1+1:)=Dimen2%DimData
		DimAdd%Dimsize=Dimen%Dimsize+Dimen2%Dimsize
		if( (Dimen%nameFlag.eq.0) .and. (Dimen2%nameFlag.eq.0) ) then
			DimAdd%nameFlag=0
			return
		end if
		if((Dimen%nameFlag.eq.0) .and. (Dimen2%nameFlag.eq.1))then
			allocate(DimAdd%DimName(l1+l2))
			do i=1,l1
				DimAdd%DimName(i)=Nameinit('0',i)
			end do
			DimAdd%DimName(l1+1:)=Dimen2%DimName
			DimAdd%nameFlag=1
			return
		end if
		if((Dimen%nameFlag.eq.0) .and. (Dimen2%nameFlag.eq.2))then
			allocate(DimAdd%DimIntName(l1+l2))
			do i=1,l1
				DimAdd%DimIntName(i)=Nameinit(0,i)
			end do
			DimAdd%DimIntName(l1+1:)=Dimen2%DimIntName
			DimAdd%nameFlag=2
			return
		end if
		if((Dimen%nameFlag.eq.1) .and. (Dimen2%nameFlag.eq.1))then
			allocate(DimAdd%DimName(l1+l2))
			DimAdd%DimName(1:l1)=Dimen%DimName
			DimAdd%DimName(l1+1:)=Dimen2%DimName
			DimAdd%nameFlag=1
			return
		end if
		if((Dimen%nameFlag.eq.1) .and. (Dimen2%nameFlag.eq.0))then
			allocate(DimAdd%DimName(l1+l2))
			DimAdd%DimName(1:l1)=Dimen%DimName
			do i=l1+1,l1+l2
				DimAdd%DimName(i)=Nameinit('0',i-l1)
			end do
			DimAdd%nameFlag=1
			return
		end if
		if((Dimen%nameFlag.eq.2) .and. (Dimen2%nameFlag.eq.2))then
			allocate(DimAdd%DimIntName(l1+l2))
			DimAdd%DimIntName(1:l1)=Dimen%DimIntName
			DimAdd%DimIntName(l1+1:)=Dimen2%DimIntName
			DimAdd%nameFlag=2
			return
		end if
		if((Dimen%nameFlag.eq.2) .and. (Dimen2%nameFlag.eq.0))then
			allocate(DimAdd%DimIntName(l1+l2))
			DimAdd%DimIntName(1:l1)=Dimen%DimIntName
			do i=l1+1,l1+l2
				DimAdd%DimIntName(i)=Nameinit(0,i-l1)
			end do
			DimAdd%nameFlag=2
			return
		end if
		write(*,*)"ERROR in DimAdd,name type different"
		stop
	end function
! The type of input are one dimension and the other vector	
	type(Dimension) function  DimAdd2(Dimen,Dimenvec)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in) :: Dimenvec(:)
		!integer::boundarysize,Dimsize
		!integer,allocatable :: Dimdata(:),boundary(:)
		integer::i,j,l1,l2
		DimAdd2%boundarysize=Dimen%boundarysize+size(Dimenvec)
		allocate(DimAdd2%boundary(DimAdd2%boundarysize))
		do i=1,Dimen%boundarysize
			DimAdd2%boundary(i)=Dimen%boundary(i)
		end do
		l1=Dimen%boundary(Dimen%boundarysize)
		do i=1,size(Dimenvec)
			DimAdd2%boundary(i+Dimen%boundarysize)=i+l1
		end do
		l1=size(Dimen%DimData)
		l2=size(Dimenvec)
		DimAdd2%lenDimData=l1+l2
		allocate(DimAdd2%Dimdata(l1+l2))
		DimAdd2%Dimsize=Dimen%Dimsize+size(Dimenvec)
		DimAdd2%Dimdata(:l1)=Dimen%DimData
		DimAdd2%Dimdata(l1+1:)=Dimenvec
		DimAdd2%nameFlag=Dimen%nameFlag
		if(Dimen%nameFlag.eq.0) return
			
		if(Dimen%nameFlag.eq.1)then
			allocate(DimAdd2%DimName(l1+l2))
			DimAdd2%DimName(1:l1)=Dimen%DimName
			do i=l1+1,l1+l2
				DimAdd2%DimName(i)=Nameinit('0',i-l1)
			end do
			return
		end if
		
		if(Dimen%nameFlag.eq.2)then
			allocate(DimAdd2%DimIntName(l1+l2))
			DimAdd2%DimIntName(1:l1)=Dimen%DimIntName
			do i=l1+1,l1+l2
				DimAdd2%DimIntName(i)=Nameinit(0,i-l1)
			end do
			return
		end if
		
		write(*,*)"ERROR in DimAdd,name type different"
		stop
		
	end function
	type(Dimension) function  DimAdd3(Dimenvec,Dimen)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in) :: Dimenvec(:)
		type(Dimension)::Dimenvec_
		Dimenvec_=Dimenvec
		DimAdd3=DimAdd(Dimenvec_,Dimen)
	end function
	
	integer function getSubDimlen(Dimen,inde)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		getSubDimlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
	end function
!		1,2,3,4,..inde,inde+1,...rank-->1,2,3,4,..(inde*inde+1*..inde+num),...rank     	
!		if inde+num>boundarysize,it will constract all the index that larger than inde
	type(Dimension) function DimConstract(dimen,inde,num)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in):: inde,num
		!integer :: boundarysize
		integer :: i!,Dimsize
		!integer,allocatable :: boundary(:)
		if(inde.gt.dimen%Dimsize) then 
			write(*,*) "ERROR in DimConstract"
			write(*,*)inde,dimen%Dimsize
			stop
		end if
		if((num.le.0).or.(inde.eq.dimen%Dimsize)) then
			DimConstract=dimen
			return
		end if
		if(inde+num.ge.dimen%boundarysize) then
			DimConstract%boundarysize=inde+1
			DimConstract%Dimsize=inde
			allocate(DimConstract%boundary(DimConstract%boundarysize))
			DimConstract%boundary(:inde)=dimen%boundary(:inde)
			DimConstract%boundary(inde+1)=dimen%boundary(dimen%boundarysize)
			allocate(DimConstract%DimData(size(dimen%DimData)))
			DimConstract%LenDimData=size(dimen%DimData)
			DimConstract%DimData=dimen%DimData
			DimConstract%Nameflag=dimen%Nameflag
			if(dimen%Nameflag.eq.0) return
			if(dimen%Nameflag.eq.1)then
				call copyNameAllocate(DimConstract%DimName,dimen%DimName)
				return
			end if
			if(dimen%Nameflag.eq.2)then
				call copyIntNameAllocate(DimConstract%DimIntName,dimen%DimIntName)
				return
			end if
			write(*,*)"ERROR in DimConstract, name flag"
			stop
		end if
		DimConstract%boundarysize=dimen%boundarysize-num
		DimConstract%Dimsize=dimen%Dimsize-num
		allocate(DimConstract%boundary(DimConstract%boundarysize))
		DimConstract%LenDimData=size(dimen%DimData)
		allocate(DimConstract%DimData(size(dimen%DimData)))
		DimConstract%boundary(:inde)=dimen%boundary(:inde)
		do i=inde+1,DimConstract%boundarysize
			DimConstract%boundary(i)=dimen%boundary(i+num)
		end do
		DimConstract%DimData=dimen%DimData
		DimConstract%Nameflag=dimen%Nameflag
		if(dimen%Nameflag.eq.0) return
		if(dimen%Nameflag.eq.1)then
			call copyNameAllocate(DimConstract%DimName,dimen%DimName)
			return
		end if
		if(dimen%Nameflag.eq.2)then
			call copyIntNameAllocate(DimConstract%DimIntName,dimen%DimIntName)
			return
		end if
		write(*,*)"ERROR in DimConstract, name flag"
		stop
	end function
     	

     	
!		(inde,inde+1,..midindde,midindde+1,...rank)-->(inde,inde+1,..midindde),(midindde+1,...rank)	     		
!		if (inde,inde+1 ..midindde),	midindde larger then next elmemt,do nothing
!		for example:	
!			D=[2,2,(3,4,5),2,(3,4)],	boundary=[0,1,2,5,6,8]
!		DimDecompose(D,3,2)=[2,2,(3,4),5,2,(3,4)],boundary=[0,1,2,4,5,6,8]
!		DimDecompose(D,3,4) will do nothing
	type(Dimension) function DimDecompose(dimen,inde,midindde)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in):: inde,midindde
		integer :: i
		if(dimen%boundary(inde) +midindde .ge.dimen%boundary(inde+1)) then
			DimDecompose=dimen
			return
		end if
		DimDecompose%boundarysize=dimen%boundarysize+1
		DimDecompose%Dimsize=dimen%Dimsize+1
		allocate(DimDecompose%boundary(DimDecompose%boundarysize))
		DimDecompose%boundary(:inde)=dimen%boundary(:inde)
		DimDecompose%boundary(inde+1)=dimen%boundary(inde) + midindde
		do i=inde+1,DimDecompose%boundarysize-1
			DimDecompose%boundary(i+1)=dimen%boundary(i)
		end do
		DimDecompose%LenDimData=size(dimen%DimData)
		allocate(DimDecompose%DimData(size(dimen%DimData)))
		DimDecompose%DimData=dimen%DimData
		DimDecompose%Nameflag=dimen%Nameflag
		if(dimen%Nameflag.eq.0) return
		if(dimen%Nameflag.eq.1)then
			call copyNameAllocate(DimDecompose%DimName,dimen%DimName)
			return
		end if
		if(dimen%Nameflag.eq.2)then
			call copyIntNameAllocate(DimDecompose%DimIntName,dimen%DimIntName)
			return
		end if
		write(*,*)"ERROR in DimDecompose, name flag"
		stop
		end function
		
	type(Dimension) function DimDecomposeAll(dimen)
		type(Dimension),intent(in) :: Dimen
		integer :: i
		if(dimen%boundarysize.eq.size(dimen%DimData)+1) then
			DimDecomposeAll=dimen
			return
     	end if
		DimDecomposeAll%boundarysize=size(dimen%DimData)+1
		DimDecomposeAll%Dimsize=size(dimen%DimData)
		allocate(DimDecomposeAll%boundary(DimDecomposeAll%boundarysize))
		do i=1,DimDecomposeAll%boundarysize
			DimDecomposeAll%boundary(i)=i-1
		end do
		DimDecomposeAll%LenDimData=size(dimen%DimData)
		allocate(DimDecomposeAll%DimData(size(dimen%DimData)))
		DimDecomposeAll%DimData=dimen%DimData
		DimDecomposeAll%Nameflag=dimen%Nameflag
		if(dimen%Nameflag.eq.0) return
		if(dimen%Nameflag.eq.1)then
			call copyNameAllocate(DimDecomposeAll%DimName,dimen%DimName)
			return
		end if
		if(dimen%Nameflag.eq.2)then
			call copyIntNameAllocate(DimDecomposeAll%DimIntName,dimen%DimIntName)
			return
		end if
		write(*,*)"ERROR in DimDecomposeAll, name flag"
		stop
	end function
     
     			
	type(Dimension) function Dimpermute(dimen,v)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in):: v(Dimen%dimSize)
		type(DimensionName),allocatable::dimenName(:)
		type(DimensionIntName),allocatable::dimenIntName(:)
		integer::i,datalen,subDlen,k
		integer,allocatable :: dimenVec(:)
		Dimpermute%Dimsize=Dimen%Dimsize
		Dimpermute%boundarysize=Dimen%boundarysize
		datalen=size(Dimen%DimData)
		allocate(Dimpermute%boundary(Dimpermute%boundarysize))
		Dimpermute%LenDimData=datalen
		allocate(Dimpermute%DimData(datalen))
		if(dimen%NameFlag.eq.1)then
			allocate(Dimpermute%DimName(datalen))
		end if
		if(dimen%NameFlag.eq.2)then
			allocate(Dimpermute%DimIntName(datalen))
		end if
		Dimpermute%boundary(1)=0
		k=1
		Dimpermute%Nameflag=dimen%NameFlag
		do i=1,Dimen%dimSize
			subDlen=getsubDimlen(dimen,v(i))
			Dimpermute%boundary(i+1)=Dimpermute%boundary(i)+subDlen
			call getSubDim(Dimen,v(i),dimenVec)
			Dimpermute%DimData(k:subDlen+k-1)=dimenVec
			if(dimen%NameFlag.eq.1)then
				call getSubDimName(Dimen,v(i),dimenName)
				Dimpermute%DimName(k:subDlen+k-1)=dimenName
			end if
			if(dimen%NameFlag.eq.2)then
				call getSubDimIntName(Dimen,v(i),dimenIntName)
				Dimpermute%DimIntName(k:subDlen+k-1)=dimenIntName
			end if
			k=k+subDlen
		end do
	end function

	logical function equal_of_dim(dim1,dim2)
		type(Dimension),intent(in) :: dim1,dim2
		integer,allocatable :: dimenVec1(:)
		integer,allocatable :: dimenVec2(:)
		call copydimension(dimenVec1,dim1)
		call copydimension(dimenVec2,dim2)
		equal_of_dim=equal_of_array(dimenVec1,dimenVec2)
	end function
!************* equal_of_array *************
	logical function equal_of_array(a,b)
		integer,intent(in) :: a(:),b(:)
		integer :: la,lb,i,l
		la=size(a)
		lb=size(b)
		equal_of_array=.false.
		if(la.eq.lb) then
			l=count(abs(a-b).gt.0)
			if(l.eq.0) then
				equal_of_array=.true.
			end if
		end if
		return
	end function
!******************************************************************************
!if flag=1,inv[1,2,3,4,...,ith,ith+1,...]->[ith,1,2,3,..,ith-1,ith+1,ith+2,...]
!if flag=0,inv[1,2,3,4,...,ith,ith+1,...]->[2,3,..,ith-1,ith,1,ith+1,ith+2,...]
	subroutine permuteorder(outv,inv_,ith,flag)
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
			write(*,*)"ERROR in permuteorder"
			stop
		end if
		return
	end subroutine
	logical function if_original_dim(dimen)
		type(Dimension),intent(in)::dimen
		if(size(dimen%DimData).eq.dimen%Dimsize) then	
			if_original_dim=.true.
		else
			if_original_dim=.false.
		end if
		return
	end function
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!**************  there is no complex*16 type function ***************************************
! This part written by Wenyuan Liu
!
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!!! Given tensor T, in order to decompos T in a certan way, we can set T%TenDim using this function.
!  T is a  matrix with size 8 x 8.

! eg 1:  We can use this function to let T become a 4-index tensor.  
! Input parameter: Dimen=T%TenDim, inboundary = (/0,1,2,3,4/),indimdata = (/2,2,2,2/)

! eg 2: if we set Dimen=T%TenDim, inboundary = (/0,1,2,4/),indimdata = (/2,2,2,2/)
! T will become a 3-index tensor, with each dimension size is 2, 2, 2*2.
! Now we can apply operator .dc. to the 3-index T, then T will become a 4-index  tensor.

	subroutine SetDimension(Dimen,inboundary,indimdata)
		type(Dimension),intent(inout) ::Dimen
		integer,intent(in) :: inboundary(:),indimdata(:)
		call cleanDimension(Dimen)
		Dimen%boundarysize=size(inboundary)
		Dimen%Dimsize=Dimen%boundarysize-1
		allocate(Dimen%boundary(size(inboundary)))
		allocate(Dimen%DimData(size(indimdata)))
		Dimen%boundary=inboundary
		Dimen%DimData=indimdata
	end subroutine
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************




!**********************************************************************
!**********************************************************************
!	the code below is for MPI
!**********************************************************************
	subroutine sent_dimensionName(Name1,Name2,ID1,ID2,ierr)
		type(DimensionName),intent(in)::Name1
		type(DimensionName),intent(inout)::Name2
		integer,intent(in)::ID1,ID2,ierr
		integer::proID,proNum,tag,istatus(MPI_STATUS_SIZE)
		tag=1
		call mpi_comm_rank(mpi_comm_world,proID,ierr)
		call mpi_comm_size(mpi_comm_world,proNum,ierr )
		
		if(ID1.eq.ID2) return !The some cpu, do nothing
		
		if((proID.ne.ID1).and.(proID.ne.ID2)) return!The proID do not sent or recv, return
		
!************************   TensorName   *************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Name1%TensorName,len_of_Name,MPI_CHARACTER,ID2,tag,MPI_Comm_world,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Name2%TensorName,len_of_Name,MPI_CHARACTER,ID1,tag,MPI_Comm_world,istatus,ierr)
		end if
!*************************   DimenName   **************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Name1%DimenName,len_of_Name,MPI_CHARACTER,ID2,tag,MPI_Comm_world,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Name2%DimenName,len_of_Name,MPI_CHARACTER,ID1,tag,MPI_Comm_world,istatus,ierr)
		end if
		return
	end subroutine
	
	subroutine BCAST_DimensionName(DimenName,ID,ierr)
		type(DimensionName),intent(inout)::DimenName
		integer,intent(in)::ID,ierr
		integer::proID,proNum,tag,len1,istatus(MPI_STATUS_SIZE)
		tag=1
		call mpi_comm_rank(mpi_comm_world,proID,ierr)
		call mpi_comm_size(mpi_comm_world,proNum,ierr )
!************************   TensorName   *************************************************		
		call MPI_BCAST(DimenName%TensorName,len_of_Name,MPI_CHARACTER,ID,MPI_COMM_WORLD,ierr)	
!*************************   DimenName   **************************************************		
		call MPI_BCAST(DimenName%DimenName,len_of_Name,MPI_CHARACTER,ID,MPI_COMM_WORLD,ierr)	
		return
	end subroutine
	
	subroutine sent_Dimension(Dimen1,Dimen2,ID1,ID2,ierr)
		type(Dimension),intent(in)::Dimen1
		type(Dimension),intent(inout)::Dimen2
		integer,intent(in)::ID1,ID2,ierr
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),i
		tag=1
		len2=0
		call mpi_comm_rank(mpi_comm_world,proID,ierr)
		call mpi_comm_size(mpi_comm_world,proNum,ierr )
		
		if(ID1.eq.ID2) return !The some cpu, do nothing
		
		if((proID.ne.ID1).and.(proID.ne.ID2)) return!The proID do not sent or recv, return
		
		
!*********************   boundarysize   ******************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%boundarysize,1,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
			if(Dimen1%boundarysize.eq.0) return !No date,return
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%boundarysize,1,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
			if(Dimen2%boundarysize.eq.0) then
				call cleanDimension(Dimen2)
				return !No date,return
			end if
		end if
!***********************    Dimsize   ****************************************************				
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%Dimsize,1,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
			if(Dimen1%Dimsize.eq.0) return !No date,return
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%Dimsize,1,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
			if(Dimen2%Dimsize.eq.0) then
				call cleanDimension(Dimen2)
				return !No date,return
			end if
		end if
!*********************   boundary   *******************************************************		
		if(proID.eq.ID1) then
			len1=size(Dimen1%boundary)
			call mpi_send(len1,1,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(len1,1,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
			!if(allocated(Dimen2%boundary)) then
			!	deallocate(Dimen2%boundary)
			!end if
			if(allocated(Dimen2%boundary)) then
				if(len1.ne.size(Dimen2%boundary))then
					deallocate(Dimen2%boundary)
					allocate(Dimen2%boundary(len1))
				end if
			else
				allocate(Dimen2%boundary(len1))
			end if
		end if
		
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%boundary,len1,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%boundary,len1,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
		end if
!******************   DimData  ***********************************************************		
		if(proID.eq.ID1) then
			len2=size(Dimen1%DimData)
			call mpi_send(len2,1,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(len2,1,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
		!	if(allocated(Dimen2%DimData)) then
		!		deallocate(Dimen2%DimData)
		!	end if
			if(allocated(Dimen2%DimData)) then
				if(len2.ne.size(Dimen2%DimData))then
					deallocate(Dimen2%DimData)
					allocate(Dimen2%DimData(len2))
				end if
			else
				allocate(Dimen2%DimData(len2))
			end if
		end if
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%DimData,len2,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%DimData,len2,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
		end if
!******************   NameFlag  ***********************************************************
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%NameFlag,1,MPI_logical,ID2,tag,MPI_Comm_world,ierr)
			if(Dimen1%NameFlag.eq.0) return!No Name,return
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%NameFlag,1,MPI_logical,ID1,tag,MPI_Comm_world,istatus,ierr)
		end if
!******************   dimensionName  ***********************************************************
		if(proID.eq.ID2) then
			if(Dimen2%NameFlag.eq.0) then
				call cleanDimensionName(Dimen2)
				return!No Name,return
			end if
			
			if(allocated(Dimen2%DimName)) then
				if(len2.ne.size(Dimen2%DimName))then
					deallocate(Dimen2%DimName)
					allocate(Dimen2%DimName(len2))
				end if
			else
				allocate(Dimen2%DimName(len2))
			end if
		end if
		do i=1,len2
			call sent_dimensionName(Dimen1%DimName(i),Dimen2%DimName(i),ID1,ID2,ierr)
		end do
		return
	end subroutine
	
	subroutine BCAST_Dimension(Dimen1,ID,ierr)
		type(Dimension),intent(inout)::Dimen1
		integer,intent(in)::ID,ierr
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),i
		tag=1
		call mpi_comm_rank(mpi_comm_world,proID,ierr)
		call mpi_comm_size(mpi_comm_world,proNum,ierr )
!*******************   boundarysize   **********************************************************		
		call MPI_BCAST(Dimen1%boundarysize,1,MPI_integer,ID,MPI_COMM_WORLD,ierr)	
		if(Dimen1%boundarysize.eq.0) return !No date,return
!********************   Dimsize   ****************************************************	
		call MPI_BCAST(Dimen1%Dimsize,1,MPI_integer,ID,MPI_COMM_WORLD,ierr)	
		if(Dimen1%Dimsize.eq.0) return !No date,return
!*********************   boundary   ***************************************************		
		if(proID.eq.ID) then
			len1=size(Dimen1%boundary)
		end if
		call MPI_BCAST(len1,1,MPI_integer,ID,MPI_COMM_WORLD,ierr)	
		if(len1.eq.0) then
			write(*,*)"ERROR in BCAST_Dimension,len1=0"
			return
		end if	
		if(proID.ne.ID) then
		!	if(allocated(Dimen1%boundary)) then
		!		deallocate(Dimen1%boundary)
		!	end if
			if(allocated(Dimen1%boundary)) then
				if(len1.ne.size(Dimen1%boundary))then
					deallocate(Dimen1%boundary)
					allocate(Dimen1%boundary(len1))
				end if
			else
				allocate(Dimen1%boundary(len1))
			end if
		end if
		call MPI_BCAST(Dimen1%boundary,len1,MPI_integer,ID,MPI_COMM_WORLD,ierr)		
!*****************    DimData   *********************************************************		
		if(proID.eq.ID) then
			len2=size(Dimen1%DimData)
		end if
		call MPI_BCAST(len2,1,MPI_integer,ID,MPI_COMM_WORLD,ierr)		
		if(len2.eq.0) then
			write(*,*)"ERROR in BCAST_Dimension,len2=0"
			return
		end if
		if(proID.ne.ID) then
		!	if(allocated(Dimen1%DimData)) then
		!		deallocate(Dimen1%DimData)
		!	end if
			if(allocated(Dimen1%DimData)) then
				if(len2.ne.size(Dimen1%DimData))then
					deallocate(Dimen1%DimData)
					allocate(Dimen1%DimData(len2))
				end if
			else
				allocate(Dimen1%DimData(len2))
			end if
		end if
		call MPI_BCAST(Dimen1%DimData,len2,MPI_integer,ID,MPI_COMM_WORLD,ierr)	
!******************   NameFlag  ***********************************************************
		call MPI_BCAST(Dimen1%NameFlag,1,MPI_logical,ID,MPI_COMM_WORLD,ierr)	
!******************   dimensionName  ***********************************************************
		if(Dimen1%NameFlag.eq.1)then
			if(proID.ne.ID) then
				if(allocated(Dimen1%DimName)) then
					if(len2.ne.size(Dimen1%DimName))then
						deallocate(Dimen1%DimName)
						allocate(Dimen1%DimName(len2))
					end if
				else
					allocate(Dimen1%DimName(len2))
				end if
			end if
			
			do i=1,len2
				call BCAST_DimensionName(Dimen1%DimName(i),ID,ierr)
			end do
		else
				call cleanDimensionName(Dimen1)
		end if
		return
	end subroutine
end module
!****************************************************
!*************** ENF OF Dimension *************
!****************************************************























