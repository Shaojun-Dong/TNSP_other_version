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

!
!by default, the DimensionName:TensorName='0',dimenName=the order of the index
!            DimensionIntName:TensorName=[0,0,0],dimenName=[-i,0,0],i is the order of the index
!
module Dimension_typede
	use usefull_function
	implicit none
	CHARACTER*1,save,private::indexsymbol='.'
	CHARACTER*1,save,private::intNamesymbol='_'!!0_0_1.0_1_2
	integer,private,save::len_of_intTensorName=1
	integer,private,save::len_of_intDimenName=2
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
		integer,private :: boundarysize=0!The length of the boundary that is useful
		integer,private :: lenDimData=0!The length of the DimData that is useful,If it is 0, means no data, use it to tell if allocate dimension
		integer,private :: Dimsize=0!The Rank of The Tensor
		integer,allocatable,private :: boundary(:)!only boundary(1:boundarysize) are the usefull data
		integer,allocatable,private :: DimData(:)!only DimData(1:lenDimData) are the usefull data
		type(dimensionName),allocatable,private::DimName(:)!only DimName(1:lenDimData) are the usefull data
		type(DimensionIntName),allocatable,private::DimIntName(:)!only DimIntName(1:lenDimData) are the usefull data
		integer,private::nameflag=0!=0  means no name
		                           !=1  means there are names that is CHARACTER(DimName)
		                           !=2  means there are names that is integer(DimintName)
		logical,private::sample_dimension_flag=.true.!if true, only store DimData, DimName, DimIntName, nameflag and lenDimData
		                                             !When do not do the DimConstract, DimDecompose and
		                                             !no name, it gose faster.
	end type Dimension
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
	
	interface allocateCheckName
		module procedure  allocateCheckName1
		module procedure  allocateCheckNameint
	end interface
!**********************************************************
!	Other Function or Subroutine:
!
!		if_sample_dimension
!
!		emptyDimension
!
!		set_len_of_intName
!
!		set_sample_dimension_flag(dimen)
!
!		outTenNameLen
!
!		outDimNameLen
!
!		set_sample_dimension(dimen)
!
!		outIntNameDim
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
!				call mpi_comm_rank(mpi_comm_worl d,proID,ierr)
!				call mpi_comm_size(mpi_comm_worl d,proNum,ierr ) 
!				call sent_Dimension(T1,T2,0,1,ierr)  !T1 in cpu0 sent to cpu2 and store in T2
!				call BCAST_Dimension(T1,0,ierr) 		 !T1 in cpu 0 send to every cpus,store in T1 in other cpus
!
!**********************************************************	
	contains
	logical function if_sample_dimension(dimen)
		type(Dimension),intent(in)::dimen
		if_sample_dimension=dimen%sample_dimension_flag
		return
	end function
   subroutine set_len_of_intName(lenTenD,lendimN)
      integer,intent(in)::LenTenD,lendimN
      if(lenTenD.gt.len_of_intName_in_type_define)then 
         write(*,*)"The len_of_lenTenD is larger than the max length"
         call error_stop()
       end if
       if(lendimN.gt.len_of_intName_in_type_define)then 
         write(*,*)"The len_of_lenTenD is larger than the max length"
         call error_stop()
       end if
      len_of_intTensorName=lenTenD
      len_of_intDimenName=lendimN
      return
   end subroutine
   integer function outTenNameLen()
   	outTenNameLen=len_of_intTensorName
   	return
   end function
   integer function outDimNameLen()
   	outDimNameLen=len_of_intDimenName
   	return
   end function
   subroutine set_sample_dimension_flag(dimen)!dimen%lenDimData,Dimen%DimData are already in dimen
   	type(Dimension),intent(inout)::dimen
   	integer::i
   	if(dimen%sample_dimension_flag)then
			dimen%sample_dimension_flag=.false.
			dimen%Dimsize=dimen%lenDimData
			dimen%boundarysize=dimen%lenDimData+1
			call allocateCheck(dimen%boundary,dimen%boundarysize)
			do i=1,Dimen%boundarysize
				Dimen%boundary(i)=i-1
			end do
			Dimen%nameflag=0
		end if
		return
	end subroutine
	subroutine set_sample_dimension(dimen)!make the dimension as a sample_dimension, and kill dimensionName
   	type(Dimension),intent(inout)::dimen
   	integer::i
   	if(dimen%sample_dimension_flag) return
   	
   	do i=1,dimen%Dimsize
   		dimen%DimData(i)=dimen.i.i !Will not gose wrong,as the index of dimen.i.i will larger then dimen%DimData(i)
   	end do
   	dimen%LenDimData=dimen%Dimsize
		dimen%sample_dimension_flag=.true.
		dimen%Dimsize=0
		dimen%boundarysize=0
		Dimen%nameflag=0
		return
	end subroutine
	subroutine outDimension_boundry(dimen1,boundary,boundarysize,Dimsize)
   	type(Dimension),intent(in)::dimen1
   	integer,intent(inout)::boundarysize,Dimsize
   	integer,allocatable,intent(inout)::boundary(:)
   	integer::i
   	if(dimen1%sample_dimension_flag)then
			Dimsize=dimen1%lenDimData
			boundarysize=dimen1%lenDimData+1
			call allocateCheck(boundary,boundarysize)
			do i=1,boundarysize
				boundary(i)=i-1
			end do
		else
			Dimsize=dimen1%Dimsize
			boundarysize=dimen1%boundarysize
			call allocateCheck(boundary,boundarysize)
			do i=1,Dimen1%boundarysize
				boundary(i)=Dimen1%boundary(i)
			end do
		end if
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
		intNameinit1%TensorName(1:len_of_intTensorName)=TensorName(1:len_of_intTensorName)
		intNameinit1%DimenName(1:len_of_intDimenName)=DimenName(1:len_of_intDimenName)
		return
	end function 
	type(DimensionIntName) function intNameinit2(intDimname)
		integer,intent(in)::intDimname(:,:)
		intNameinit2%TensorName(1:len_of_intTensorName)=intDimname(1,1:len_of_intTensorName)
		intNameinit2%DimenName(1:len_of_intDimenName)=intDimname(2,1:len_of_intDimenName)
		return
	end function 
	type(DimensionIntName) function intNameinit3(TensorName,DimenName)
		integer,intent(in)::TensorName(:),DimenName
		intNameinit3%TensorName(1:len_of_intTensorName)=TensorName(1:len_of_intTensorName)
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
			call error_stop()
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
			call error_stop()
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
		do i=1,len_of_intTensorName
			w=w+inName%TensorName(i)
			if(i.ne.len_of_intTensorName) w=w+intNamesymbol
		end do
		w=w+indexsymbol
		do i=1,len_of_intDimenName
			w=w+inName%DimenName(i)
			if(i.ne.len_of_intDimenName) w=w+intNamesymbol
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
			call error_stop()
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
			call error_stop()
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
			call error_stop()
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
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"The dimension should be in its original dimension"
			call error_stop()
		end if
		if(dimen%nameflag.eq.2) then
			write(*,*)"The Dimension have the name of integer"
			call error_stop()
		end if
		lenD=dimen%LenDimData
		call allocateCheckName(dimen%dimname,lenD)
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
			call error_stop()
		end if
		if(dimen%nameflag.eq.1) then
			write(*,*)"The Dimension have the name of character"
			call error_stop()
		end if
		lenD=dimen%LenDimData
		call allocateCheckName(dimen%DimintName,lenD)
		if(dimen%nameflag.eq.2)then!There are name in the dimension,rename
			do i=1,lenD
				dimen%DimintName(i)%TensorName=TensorName
			end do
		else
			do i=1,lenD
				dimen%DimintName(i)=Nameinit(TensorName,-i)
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
		Dimen%nameflag=0
		if(allocated(Dimen%boundary)) then
			deallocate(Dimen%boundary)
		end if
		if(allocated(Dimen%DimData)) then
			deallocate(Dimen%DimData)
		end if
		call cleanDimensionName(Dimen)
		Dimen%sample_dimension_flag=.true.
		return
	end subroutine
		subroutine emptyDimension(Dimen)!make dimension to a a empty dimension.but do not deallocate
		type(Dimension),intent(inout)::	Dimen
		Dimen%boundarysize=0
		Dimen%Dimsize=0
		Dimen%lenDimData=0
		Dimen%nameflag=0
		Dimen%sample_dimension_flag=.true.
		return
	end subroutine
	
	subroutine outDimData(DimData,dimen)
		type(Dimension),intent(in)::dimen
		integer,intent(inout)::DimData(:)
		DimData=dimen%DimData(1:dimen%lenDimData)
		return
	end subroutine
!****************print the information of dimension ****************************************	
	subroutine Dprint1(Dimen)
		type(Dimension),intent(in) ::Dimen
		integer,allocatable :: dimenVec(:)
		CHARACTER*5000,allocatable::w(:)
		integer::i
		write(*,*) "***   START   ***"
		write(*,*) Dimen%DimData(1:dimen%lenDimData)
		if(.not.Dimen%sample_dimension_flag) then
			if(Dimen%nameflag.eq.1)then
				allocate(w(Dimen%lenDimData))
				w=Dimen%DimName(1:Dimen%lenDimData)
				write(*,*)"index Name are"
				write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
			end if
			if(Dimen%nameflag.eq.2)then
				allocate(w(Dimen%lenDimData))
				w=Dimen%DimintName(1:Dimen%lenDimData)
				write(*,*)"index Name are"
				write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
			end if
			call copydimension(dimenVec,Dimen)
			write(*,*) "			--				"
			write(*,*)	dimenVec
		end if
		write(*,*) "***   END   ***"
	end subroutine
				
	subroutine Dprint0(Dimen)
		type(Dimension),intent(in) ::Dimen
		integer,allocatable :: dimenVec(:)
		CHARACTER*5000,allocatable::w(:)
		integer::i
		write(*,*) "***   Dimension Data    ***"
		write(*,*) Dimen%DimData(1:dimen%lenDimData)
		if(.not.Dimen%sample_dimension_flag)then
			call copydimension(dimenVec,Dimen)
			write(*,*)	dimenVec
			write(*,*) "***   Dimension END   ***"
		else
			write(*,*) "***   Dimension END   ***"
		end if
		if(Dimen%nameflag.eq.1)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName(1:Dimen%lenDimData)
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		if(Dimen%nameflag.eq.2)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimintName(1:Dimen%lenDimData)
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
		if(Dimen%sample_dimension_flag)then
			write(*,*)Dimen%DimData(1:Dimen%LenDimData)
		else
			call copydimension(dimenVec,Dimen)
			write(*,*)	dimenVec
			write(*,*)"It is not original dimension "
		end if
		write(*,*) "***   Dimension END   ***"
		if(Dimen%nameflag.eq.1)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName(1:Dimen%lenDimData)
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		if(Dimen%nameflag.eq.2)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimintName(1:Dimen%lenDimData)
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
		write(*,*) "***   Dimension Data    ***"
		write(*,*)Dimen%DimData(1:Dimen%LenDimData)
		if(.not.Dimen%sample_dimension_flag)then
			write(*,*) "***   START   ***"
			write(*,*) "Dimsize"
			write(*,*) Dimen%Dimsize
			write(*,*) "boundarysize"
			write(*,*) Dimen%boundarysize
			write(*,*) "boundary"
			write(*,*) Dimen%boundary(1:Dimen%boundarysize)
			call copydimension(dimenVec,Dimen)
			write(*,*) "		Dimension   	"
			write(*,*)	dimenVec
		end if
		if(Dimen%nameflag.eq.1)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName(1:Dimen%lenDimData)
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		if(Dimen%nameflag.eq.2)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimintName(1:Dimen%lenDimData)
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		
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
		Dimensioniniti%sample_dimension_flag=.false.
		return
	end function
	subroutine DimInitialization(Dimen,DimData)
		type(Dimension),intent(inout) ::Dimen
		integer,intent(in) :: DimData(:)
		integer::i
		Dimen%lenDimData=size(DimData)
		call allocateCheck(Dimen%DimData,Dimen%lenDimData)
		Dimen%DimData(1:Dimen%lenDimData)=DimData
		Dimen%sample_dimension_flag=.true.
		Dimen%nameflag=0
		return
	end subroutine
	subroutine DimInitialization2(Dimen,Dimen2)
		type(Dimension),intent(inout) ::Dimen
		type(Dimension),intent(in) ::Dimen2
		integer::lenDim,lenBoun,i
		lenDim=Dimen2%lenDimData
		if(lenDim.eq.0) then !no Data
			Dimen%Dimsize=0
			Dimen%boundarysize=0
			Dimen%lenDimData=0
			Dimen%nameflag=0
			Dimen%sample_dimension_flag=.true.
			return
		end if
		Dimen%lenDimData=lenDim
		call allocateCheck(Dimen%DimData , lenDim)
		Dimen%DimData(1:lenDim)=Dimen2%DimData(1:lenDim)
		Dimen%nameflag=Dimen2%nameflag
		Dimen%sample_dimension_flag=Dimen2%sample_dimension_flag
		if((Dimen%nameflag.eq.0).and.Dimen%sample_dimension_flag) return
		
		if(Dimen2%nameflag.eq.1)then
			call allocateCheckName(dimen%dimName,lenDim)
			do i=1,lenDim
				dimen%dimName(i)=dimen2%dimName(i)
			end do
		end if
		if(Dimen2%nameflag.eq.2)then
			call allocateCheckName(dimen%dimintName,lenDim)
			do i=1,lenDim
				dimen%dimintName(i)=dimen2%dimintName(i)
			end do
		end if
		if(.not.Dimen2%sample_dimension_flag) then
			lenBoun=Dimen2%boundarysize
			Dimen%boundarysize=lenBoun
			call allocateCheck(Dimen%boundary , lenBoun)
			Dimen%boundary(1:lenBoun)=Dimen2%boundary(1:lenBoun)
			Dimen%Dimsize=Dimen2%Dimsize
		end if
		return
	end subroutine
	integer	function DimSize(Dimen)
		type(Dimension),intent(in) :: Dimen
		if(Dimen%sample_dimension_flag)then
			DimSize=Dimen%LenDimData
		else
			DimSize=Dimen%dimSize
		end if
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
			call error_stop()
		end if
		if(ith.gt.Dimen%lenDimData)then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,Dimen%lenDimData
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"NameID is use in original dimension"
			call error_stop()
		end if
		if(Dimen%nameflag.eq.1)then
			outNameTen=Dimen%DimName(ith)%TensorName
			return
		end if
		if(Dimen%nameflag.eq.2)then
			outNameTen=''
			do i=1,len_of_intTensorName
				outNameTen=outNameTen+Dimen%DimintName(ith)%TensorName(i)
				if(i.ne.len_of_intTensorName)outNameTen=outNameTen+intNamesymbol
			end do
			return
		end if
		write(*,*)"ERROR of nameFlag"
		call error_stop()
	end function
	function outIntNameTen(dimen,ith)
		integer::outIntNameTen(len_of_intTensorName)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no integer name in the dimension"
			call error_stop()
		end if
		if(ith.gt.Dimen%LenDimData)then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,Dimen%LenDimData
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"NameID is use in original dimension"
			call error_stop()
		end if
		outIntNameTen=Dimen%DimintName(ith)%TensorName(1:len_of_intTensorName)
		return
	end function
	function outIntNameDim(dimen,ith)
		integer::outIntNameDim(len_of_intDimenName)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no integer name in the dimension"
			call error_stop()
		end if
		if(ith.gt.Dimen%LenDimData)then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,Dimen%LenDimData
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"NameID is use in original dimension"
			call error_stop()
		end if
		outIntNameDim=Dimen%DimintName(ith)%DimenName(1:len_of_intDimenName)
		return
	end function
!*******  function or subroutine for name   **************		
!output the index of the ith dimension
	CHARACTER(len=len_of_Name+len_of_Name) function outName(dimen,ith)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(Dimen%nameflag.eq.0)then
			write(*,*)"There is no name in the dimension"
			call error_stop()
		end if
		if(ith.gt.size(Dimen%DimName))then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,size(Dimen%DimName)
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"indexID is use in original dimension"
			call error_stop()
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
		call error_stop()
	end function
	function outDimintName(dimen,ith,TensorNameflag)
		integer,allocatable::outDimintName(:)
		type(Dimension),intent(in) :: Dimen
		logical,intent(in)::TensorNameflag
		integer,intent(in)::ith
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no integer name in the dimension"
			call error_stop()
		end if
		if(ith.gt.size(Dimen%DimintName))then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,size(Dimen%DimintName)
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"NameID is use in original dimension"
			call error_stop()
		end if
		if(TensorNameflag)then
			allocate(outDimintName(len_of_intTensorName))
			outDimintName=Dimen%DimintName(ith)%TensorName(1:len_of_intTensorName)
		else
			allocate(outDimintName(len_of_intDimenName))
			outDimintName=Dimen%DimintName(ith)%DimenName(1:len_of_intDimenName)
		end if
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
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"Nameorder2 is use in original dimension"
			call error_stop()
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
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"Nameorder3 is use in original dimension"
			call error_stop()
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
		integer,intent(in)::TensorName(:),dimenName(:)
		integer::i
		type(DimensionIntName)::testName 
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,Nameorder4"
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"Nameorder4 is use in original dimension"
			call error_stop()
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
		integer,intent(in)::TensoDimrName(:,:)
		integer::i
		type(DimensionIntName)::testName 
		if(Dimen%nameflag.ne.2)then
			write(*,*)"There is no int name in the dimension,Nameorder5"
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"Nameorder5 is use in original dimension"
			call error_stop()
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
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"Nameorder5 is use in original dimension"
			call error_stop()
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
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"NameorderArray is use in original dimension"
			call error_stop()
		end if
		lenn=size(w)
		if(lenn.gt.dimen%lenDimData)then
			write(*,*)"ERROR i NameorderArray"
			write(*,*)lenn,dimen%Dimsize
			call error_stop()
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
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"intNameorderArray is use in original dimension"
			call error_stop()
		end if
		lenn=size(TensorName,1)
		if(lenn.gt.dimen%lenDimData)then
			write(*,*)"ERROR i intNameorderArray"
			write(*,*)lenn,dimen%Dimsize
			call error_stop()
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
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"intNameorderArray is use in original dimension"
			call error_stop()
		end if
		lenn=size(TensorintName)
		if(lenn.gt.dimen%lenDimData)then
			write(*,*)"ERROR i intNameorderArray"
			write(*,*)lenn,dimen%Dimsize
			call error_stop()
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
		lg1=equal_of_array(name1%TensorName(:len_of_intTensorName),name2%TensorName(:len_of_intTensorName))
		lg2=equal_of_array(name1%dimenName(:len_of_intDimenName),name2%dimenName(:len_of_intDimenName))
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
			call error_stop()
		end if
		oldfullname=chartoName_log(oldTensorName,olddimenName,oldname)
		!(TensorName,DimenName,w_)
		newfullname=chartoName_log(newTensorName,newdimenName,newname)
		if(oldfullname.neqv.newfullname)then
			write(*,*)"ERROR in resetIndexname1"
			call error_stop()
		end if
		if(newfullname)then
			do i=1,Dimen%lenDimData
				if(dimen%dimname(i)%TensorName.equ.oldTensorName)then
					if(dimen%dimname(i)%dimenName.equ.olddimenName)then
						dimen%dimname(i)%TensorName=newTensorName
						dimen%dimname(i)%dimenName=newdimenName
					end if
				end if
			end do
		else
			do i=1,Dimen%lenDimData
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
		if(ith.gt.dimen%LenDimData)then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,dimen%LenDimData
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"resetIndexNameID is use in original dimension"
			call error_stop()
		end if
		fullname=chartoName_log(newTensorName,newdimenName,newname)
		if(Dimen%nameflag.eq.0)then!There is no name in dimension
			lenName=dimen%LenDimData
			!allocate(dimen%dimname(lenName))
			call allocateCheckName(dimen%dimname,lenName)
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
			call error_stop()
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
			call error_stop()
		end if
		do i=1,dimen%lenDimData
			if(dimen%dimintname(i)%TensorName(:len_of_intTensorName).equ.oldnameTen)then
				dimen%dimintname(i)%TensorName(:len_of_intTensorName)=newnameTen
			end if
		end do
		return
	end subroutine
	subroutine setDimName5(dimen,ith,newnameTen) 
		type(Dimension),intent(inout) :: dimen
		integer,intent(in)::ith
		integer,intent(in)::newnameTen(:)
		integer::i
		if(Dimen%nameflag.eq.1)then
			write(*,*)"There is cahr name in the dimension,resetindexname"
			call error_stop()
		end if
		if(ith.gt.dimen%LenDimData)then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,dimen%LenDimData
			call error_stop()
		end if
		if(Dimen%nameflag.eq.2) then
			dimen%dimintname(ith)%TensorName(:len_of_intTensorName)=newnameTen
			return
		end if
		if(Dimen%nameflag.eq.0)then
			call allocateCheckName(dimen%dimintname,dimen%lenDimData)
			do i=1,dimen%lenDimData
				dimen%dimintname(i)=Nameinit(newnameTen,-i)
			end do
			Dimen%nameflag=2
		end if
		return
	end subroutine
	subroutine setDimName6(dimen,ith,newnameTen,newnameDim)
		type(Dimension),intent(inout) :: dimen
		integer,intent(in)::ith
		integer,intent(in)::newnameTen(:),newnameDim(:)
		integer::i
		if(Dimen%nameflag.eq.1)then
			write(*,*)"There is char name in the dimension,resetindexname"
			call error_stop()
		end if
		if(ith.gt.dimen%LenDimData)then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,dimen%LenDimData
			call error_stop()
		end if
		if(Dimen%nameflag.eq.2) then
			dimen%dimintname(ith)%TensorName(:len_of_intTensorName)=newnameTen
			dimen%dimintname(ith)%DimenName(:len_of_intDimenName)=newnameDim
			return
		end if
		if(Dimen%nameflag.eq.0)then
			call allocateCheckName(dimen%dimintname,dimen%lenDimData)
			do i=1,dimen%lenDimData
				if(i.eq.ith) then
					dimen%dimintname(i)=Nameinit(newnameTen,newnameDim)
				else
					dimen%dimintname(i)=Nameinit(newnameTen,-i)
				end if
			end do
			Dimen%nameflag=2
		end if
		return
	end subroutine
	
	integer function outtotalData(Dimen)
		type(Dimension),intent(in) :: Dimen
		outtotalData=product(Dimen%DimData(1:dimen%lenDimData))
		return
	end function				
!	return the inde dimension	,outpout in a integer
	integer function Dim_i(Dimen,inde)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		integer :: i,D1
		if(Dimen%sample_dimension_flag)then
			if(inde.gt.dimen%LenDimData) then 
				write(*,*) "ERROR in Dim_i"
				write(*,*)"stop"
				call error_stop()
			end if
			Dim_i=Dimen%DimData(inde)
			return
		end if
		
		if(inde.gt.dimen%Dimsize) then 
			write(*,*) "ERROR in Dim_i"
			write(*,*)"stop"
			call error_stop()
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
		if(Dimen%sample_dimension_flag)then
			if(allocated(dimenVec)) then
				if(dimen%LenDimData.ne.size(dimenVec)) then 
					deallocate(dimenVec)
					allocate(dimenVec(dimen%LenDimData))
				end if
			else
				allocate(dimenVec(dimen%LenDimData))
			end if
			dimenVec=dimen%DimData
			return
		end if
		
		if(allocated(dimenVec)) then
			if(dimen%DimSize.ne.size(dimenVec)) then 
				deallocate(dimenVec)
				allocate(dimenVec(dimen%DimSize))
			end if
		else
			allocate(dimenVec(dimen%DimSize))
		end if
		do i=1,Dimen%Dimsize
			dimenVec(i)=Dim_i(Dimen,i)
		end do
		return
	end 	subroutine		
	subroutine copyDimToVec(dimenVec,Dimen)
		integer,intent(inout) :: dimenVec(:)
		type(Dimension),intent(in) :: Dimen
		integer :: i
		if(Dimen%sample_dimension_flag)then
			if(size(dimenVec).lt.dimen%LenDimData)then
				write(*,*)"ERROR in assignment dimension to array "
				write(*,*)"array(:)=dimension(:),size(array)<size(dimension)"
				write(*,*)size(dimenVec),Dimen%LenDimData
				call error_stop()
			end if
			dimenVec=dimen%DimData
			return
		end if
		
		if(size(dimenVec).lt.Dimen%Dimsize)then
			write(*,*)"ERROR in assignment dimension to array "
			write(*,*)"array(:)=dimension(:),size(array)<size(dimension)"
			write(*,*)size(dimenVec),Dimen%Dimsize
			call error_stop()
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
		if(Dimen%sample_dimension_flag)then
			if(allocated(dimenVec)) then
				if(size(dimenVec).ne.1) then 
					deallocate(dimenVec)
					allocate(dimenVec(1))
				end if
			else
				allocate(dimenVec(1))
			end if
			dimenVec(1)=Dimen%DimData(inde)
			return
		end if
		
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
		if(Dimen%sample_dimension_flag)then
			if(inde.gt.dimen%LenDimData) then 
				write(*,*) "ERROR in getSubDimName"
				return
			end if
			if(allocated(dimenName)) then
				if(1.ne.size(dimenName)) then 
					deallocate(dimenName)
					allocate(dimenName(1))
				end if
			else
				allocate(dimenName(1))
			end if
			dimenName(1)=Dimen%DimName(inde)
			return
		end if
		
		Dlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
		if(inde.gt.dimen%dimSize) then 
			write(*,*) "ERROR in getSubDimName"
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
		if(Dimen%sample_dimension_flag)then
			if(inde.gt.dimen%LenDimData) then 
				write(*,*) "ERROR in getSubDimIntName"
				return
			end if
			if(allocated(dimenintName)) then
				if(1.ne.size(dimenintName)) then 
					deallocate(dimenintName)
					allocate(dimenintName(1))
				end if
			else
				allocate(dimenintName(1))
			end if
			dimenintName(1)=Dimen%DimIntName(inde)
			return
		end if
		
		Dlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
		if(inde.gt.dimen%dimSize) then 
			write(*,*) "ERROR in getSubDimIntName"
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
			call error_stop()
		end if
		lenD=dimen%LenDimData
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
		if(.not.Dimen%sample_dimension_flag)then
			RNDim%boundarysize=RNDim%Dimsize+1
			allocate(RNDim%boundary(RNDim%boundarysize))
			do i=1,RNDim%boundarysize
				RNDim%boundary(i)=i-1
			end do
		end if
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
		Dlen=size(Dimdata)
		if(Dlen.eq.1) then
			allocate(getSubDim2%DimData(Dlen))
			getSubDim2%DimData=DimData
			getSubDim2%lenDimData=Dlen
			getSubDim2%sample_dimension_flag=.true.
		else
			boundary(1)=0
			boundary(2)=Dlen
			boundarysize=2
			getSubDim2=Dimensioniniti(boundarysize,1,boundary,DimData)
		end if
		getSubDim2%NameFlag=Dimen%NameFlag
		if(getSubDim2%NameFlag.eq.0) return
		if(Dimen%nameFlag.eq.1)then
			call getSubDimName(Dimen,inde,getSubDim2%dimName)
			return
		end if
		if(Dimen%nameFlag.eq.2)then
			call getSubDimIntName(Dimen,inde,getSubDim2%dimIntName)
			return
		end if
		write(*,*)"ERROR in getSubDim2,nameFlag"
		call error_stop()
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
		integer,allocatable::boundary1(:),boundary2(:)
		integer::boundarysize1,boundarysize2,Dimsize1,Dimsize2
		integer::i,j,l1,l2,templen
		l1=dimen%LenDimData
		l2=dimen2%LenDimData
		DimAdd%lenDimData=l1+l2
		allocate(DimAdd%Dimdata(l1+l2))
		DimAdd%Dimdata(:l1)=Dimen%DimData(1:Dimen%lenDimData)
		DimAdd%Dimdata(l1+1:)=Dimen2%DimData(1:Dimen2%lenDimData)
		if( (Dimen%nameFlag.eq.0) .and. (Dimen2%nameFlag.eq.0) ) then
			DimAdd%nameFlag=0
			if(Dimen%sample_dimension_flag.and.Dimen2%sample_dimension_flag)then
				DimAdd%sample_dimension_flag=.true.
				return
			end if
		end if
		if((Dimen%nameFlag.eq.0) .and. (Dimen2%nameFlag.eq.1))then
			allocate(DimAdd%DimName(l1+l2))
			do i=1,l1
				DimAdd%DimName(i)=Nameinit('0',i)
			end do
			DimAdd%DimName(l1+1:)=Dimen2%DimName(1:l2)
			DimAdd%nameFlag=1
		end if
		if((Dimen%nameFlag.eq.0) .and. (Dimen2%nameFlag.eq.2))then
			allocate(DimAdd%DimIntName(l1+l2))
			do i=1,l1
				DimAdd%DimIntName(i)=Nameinit(0,-i)
			end do
			DimAdd%DimIntName(l1+1:)=Dimen2%DimIntName(1:l2)
			DimAdd%nameFlag=2
		end if
		if((Dimen%nameFlag.eq.1) .and. (Dimen2%nameFlag.eq.1))then
			allocate(DimAdd%DimName(l1+l2))
			DimAdd%DimName(1:l1)=Dimen%DimName(1:l1)
			DimAdd%DimName(l1+1:)=Dimen2%DimName(1:l2)
			DimAdd%nameFlag=1
		end if
		if((Dimen%nameFlag.eq.1) .and. (Dimen2%nameFlag.eq.0))then
			allocate(DimAdd%DimName(l1+l2))
			DimAdd%DimName(1:l1)=Dimen%DimName(1:l1)
			do i=l1+1,l1+l2
				DimAdd%DimName(i)=Nameinit('0',i)
			end do
			DimAdd%nameFlag=1
		end if
		if((Dimen%nameFlag.eq.2) .and. (Dimen2%nameFlag.eq.2))then
			allocate(DimAdd%DimIntName(l1+l2))
			DimAdd%DimIntName(1:l1)=Dimen%DimIntName(1:l1)
			DimAdd%DimIntName(l1+1:)=Dimen2%DimIntName(1:l2)
			DimAdd%nameFlag=2
		end if
		if((Dimen%nameFlag.eq.2) .and. (Dimen2%nameFlag.eq.0))then
			allocate(DimAdd%DimIntName(l1+l2))
			DimAdd%DimIntName(1:l1)=Dimen%DimIntName(1:l1)
			do i=l1+1,l1+l2
				DimAdd%DimIntName(i)=Nameinit(0,-i)
			end do
			DimAdd%nameFlag=2
		end if
		if(Dimen%sample_dimension_flag.and.Dimen2%sample_dimension_flag)then
			DimAdd%sample_dimension_flag=.true.
			return
		end if
		DimAdd%sample_dimension_flag=.false.
		call outDimension_boundry(dimen,boundary1,boundarysize1,Dimsize1)
		call outDimension_boundry(dimen2,boundary2,boundarysize2,Dimsize2)
		DimAdd%boundarysize=boundarysize1+boundarysize2-1
		allocate(DimAdd%boundary(DimAdd%boundarysize))
		do i=1,boundarysize1
			DimAdd%boundary(i)=boundary1(i)
		end do
		templen=boundary1(boundarysize1)
		do i=1,boundarysize2-1
			DimAdd%boundary(i+boundarysize1)=boundary2(i+1)+templen
		end do
		DimAdd%Dimsize=Dimsize1+Dimsize2
		return
	end function
! The type of input are one dimension and the other vector	
	type(Dimension) function  DimAdd2(Dimen,Dimenvec)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in) :: Dimenvec(:)
		integer,allocatable::boundary(:)
		integer::boundarysize,Dimsize
		integer::i,j,l1,l2,tempLen
		l1=dimen%LenDimData
		l2=size(Dimenvec)
		DimAdd2%lenDimData=l1+l2
		allocate(DimAdd2%Dimdata(l1+l2))
		DimAdd2%Dimdata(:l1)=Dimen%DimData(1:dimen%lenDimData)
		DimAdd2%Dimdata(l1+1:)=Dimenvec
		DimAdd2%nameFlag=Dimen%nameFlag
		DimAdd2%sample_dimension_flag=Dimen%sample_dimension_flag
		if((DimAdd2%nameFlag.eq.0).and.(DimAdd2%sample_dimension_flag)) return
				
			
		if(Dimen%nameFlag.eq.1)then
			allocate(DimAdd2%DimName(l1+l2))
			DimAdd2%DimName(1:l1)=Dimen%DimName(1:l1)
			do i=l1+1,l1+l2
				DimAdd2%DimName(i)=Nameinit('0',i)
			end do
		end if
		
		if(Dimen%nameFlag.eq.2)then
			allocate(DimAdd2%DimIntName(l1+l2))
			DimAdd2%DimIntName(1:l1)=Dimen%DimIntName(1:l1)
			do i=l1+1,l1+l2
				DimAdd2%DimIntName(i)=Nameinit(0,-i)
			end do
		end if
		
		if(DimAdd2%sample_dimension_flag) return
		
		call outDimension_boundry(dimen,boundary,boundarysize,Dimsize)
		
		DimAdd2%boundarysize=boundarysize+size(Dimenvec)
		allocate(DimAdd2%boundary(DimAdd2%boundarysize))
		do i=1,boundarysize
			DimAdd2%boundary(i)=boundary(i)
		end do
		tempLen=boundary(boundarysize)
		do i=1,size(Dimenvec)
			DimAdd2%boundary(i+boundarysize)=i+tempLen
		end do
		DimAdd2%Dimsize=Dimsize+size(Dimenvec)
		return
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
		if(Dimen%sample_dimension_flag)then
			getSubDimlen=1
			return
		end if
		getSubDimlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
	end function
!		1,2,3,4,..inde,inde+1,...rank-->1,2,3,4,..(inde*inde+1*..inde+num),...rank     	
!		if inde+num>boundarysize,it will constract all the index that larger than inde
	type(Dimension) function DimConstract(dimen,inde,num)
		type(Dimension)::Dimen
		integer,intent(in):: inde,num
		integer :: i,boundarysize,Dimsize
		integer,allocatable::boundary(:)
		if(dimen%sample_dimension_flag)then
			if(inde.gt.dimen%LenDimData) then 
				write(*,*) "ERROR in DimConstract"
				write(*,*)inde,dimen%LenDimData
				call error_stop()
			end if
			if((num.le.0).or.(inde.eq.dimen%LenDimData)) then
				DimConstract=dimen
				return
			end if
		else
			if(inde.gt.dimen%Dimsize) then 
				write(*,*) "ERROR in DimConstract"
				write(*,*)inde,dimen%Dimsize
				call error_stop()
			end if
			if((num.le.0).or.(inde.eq.dimen%Dimsize)) then
				DimConstract=dimen
				return
			end if
		end if
		call outDimension_boundry(dimen,boundary,boundarysize,Dimsize)
		DimConstract%sample_dimension_flag=.false.
		if(inde+num.ge.boundarysize) then
			DimConstract%boundarysize=inde+1
			DimConstract%Dimsize=inde
			allocate(DimConstract%boundary(DimConstract%boundarysize))
			DimConstract%boundary(:inde)=boundary(:inde)
			DimConstract%boundary(inde+1)=boundary(boundarysize)
			DimConstract%LenDimData=dimen%LenDimData
			allocate(DimConstract%DimData(DimConstract%LenDimData))
			DimConstract%DimData=dimen%DimData(1:dimen%LenDimData)
			DimConstract%Nameflag=dimen%Nameflag
			if(dimen%Nameflag.eq.0) return
			if(dimen%Nameflag.eq.1)then
				call copyNameAllocate(DimConstract%DimName,dimen%DimName(1:Dimen%lenDimData))
				return
			end if
			if(dimen%Nameflag.eq.2)then
				call copyIntNameAllocate(DimConstract%DimIntName,dimen%DimIntName(1:Dimen%lenDimData))
				return
			end if
			write(*,*)"ERROR in DimConstract, name flag"
			call error_stop()
		end if
		DimConstract%boundarysize=boundarysize-num
		DimConstract%Dimsize=Dimsize-num
		allocate(DimConstract%boundary(DimConstract%boundarysize))
		DimConstract%LenDimData=dimen%LenDimData
		allocate(DimConstract%DimData(dimen%LenDimData))
		DimConstract%boundary(:inde)=boundary(:inde)
		do i=inde+1,DimConstract%boundarysize
			DimConstract%boundary(i)=boundary(i+num)
		end do
		DimConstract%DimData=dimen%DimData(1:dimen%LenDimData)
		DimConstract%Nameflag=dimen%Nameflag
		if(dimen%Nameflag.eq.0) return
		if(dimen%Nameflag.eq.1)then
			call copyNameAllocate(DimConstract%DimName,dimen%DimName(1:Dimen%lenDimData))
			return
		end if
		if(dimen%Nameflag.eq.2)then
			call copyIntNameAllocate(DimConstract%DimIntName,dimen%DimIntName(1:Dimen%lenDimData))
			return
		end if
		write(*,*)"ERROR in DimConstract, name flag"
		call error_stop()
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
		if(dimen%sample_dimension_flag)then
			write(*,*)"ERROR in DimDecompose,1"
			call error_stop()
		end if
		if(inde.gt.dimen%DimSize) then
			write(*,*)"ERROR in DimDecompose"
			call error_stop()
		end if
		if(dimen%boundary(inde) +midindde .ge.dimen%boundary(inde+1)) then
			DimDecompose=dimen
			return
		end if
		DimDecompose%LenDimData=dimen%LenDimData
		allocate(DimDecompose%DimData(dimen%LenDimData))
		DimDecompose%DimData=dimen%DimData(1:dimen%LenDimData)
		DimDecompose%Dimsize=dimen%Dimsize+1
		if(DimDecompose%LenDimData.eq.DimDecompose%Dimsize) then !It is sample dimension
			DimDecompose%sample_dimension_flag=.true.
			DimDecompose%Dimsize=0
			DimDecompose%boundarysize=0
		else
			DimDecompose%sample_dimension_flag=.false.
			DimDecompose%boundarysize=dimen%boundarysize+1
			allocate(DimDecompose%boundary(DimDecompose%boundarysize))
			DimDecompose%boundary(:inde)=dimen%boundary(:inde)
			DimDecompose%boundary(inde+1)=dimen%boundary(inde) + midindde
			do i=inde+1,DimDecompose%boundarysize-1
				DimDecompose%boundary(i+1)=dimen%boundary(i)
			end do
		end if
		DimDecompose%Nameflag=dimen%Nameflag
		if(dimen%Nameflag.eq.0) return
		if(dimen%Nameflag.eq.1)then
			call copyNameAllocate(DimDecompose%DimName,dimen%DimName(1:Dimen%lenDimData))
			return
		end if
		if(dimen%Nameflag.eq.2)then
			call copyIntNameAllocate(DimDecompose%DimIntName,dimen%DimIntName(1:Dimen%lenDimData))
			return
		end if
		write(*,*)"ERROR in DimDecompose, name flag"
		call error_stop()
		end function
		
		
!		Will not store boundary as it is sample dimension	
!		If dimen is already the simple type, do nothing	
	type(Dimension) function DimDecomposeAll(dimen)
		type(Dimension),intent(in) :: Dimen
		integer :: i
		if(dimen%sample_dimension_flag)then
			DimDecomposeAll=dimen
			return
		end if
     	DimDecomposeAll%sample_dimension_flag=.true.
     	DimDecomposeAll%LenDimData=dimen%LenDimData
     	allocate(DimDecomposeAll%DimData(dimen%LenDimData))
		DimDecomposeAll%DimData=dimen%DimData(1:dimen%LenDimData)
		DimDecomposeAll%boundarysize=0
		DimDecomposeAll%Dimsize=0
		DimDecomposeAll%Nameflag=dimen%Nameflag
		if(dimen%Nameflag.eq.0) return
		if(dimen%Nameflag.eq.1)then
			call copyNameAllocate(DimDecomposeAll%DimName,dimen%DimName(1:Dimen%lenDimData))
			return
		end if
		if(dimen%Nameflag.eq.2)then
			call copyIntNameAllocate(DimDecomposeAll%DimIntName,dimen%DimIntName(1:Dimen%lenDimData))
			return
		end if
		write(*,*)"ERROR in DimDecomposeAll, name flag"
		call error_stop()
	end function
     
     			
	type(Dimension) function Dimpermute(dimen,v)
		type(Dimension),intent(in) :: Dimen
		integer,intent(in):: v(Dimen%dimSize)
		type(DimensionName),allocatable::dimenName(:)
		type(DimensionIntName),allocatable::dimenIntName(:)
		integer::i,datalen,subDlen,k
		integer,allocatable :: dimenVec(:)
		datalen=dimen%LenDimData
		if(datalen.eq.0) then
			write(*,*)"There is no data in Dimension"
			call error_stop()
		end if
		Dimpermute%LenDimData=datalen
		allocate(Dimpermute%DimData(datalen))
		if(dimen%NameFlag.eq.1)then
			allocate(Dimpermute%DimName(datalen))
		end if
		if(dimen%NameFlag.eq.2)then
			allocate(Dimpermute%DimIntName(datalen))
		end if
		Dimpermute%Nameflag=dimen%NameFlag
		Dimpermute%sample_dimension_flag=dimen%sample_dimension_flag
		if(dimen%sample_dimension_flag)then
			do i=1,datalen
				if(dimen%NameFlag.eq.1)then
					Dimpermute%DimName(i)=dimen%DimName(v(i))
				end if
				if(dimen%NameFlag.eq.2)then
					Dimpermute%DimIntName(i)=dimen%DimIntName(v(i))
				end if
				Dimpermute%DimData(i)=dimen%DimData(v(i))
			end do
			return
		end if
		
		Dimpermute%Dimsize=Dimen%Dimsize
		Dimpermute%boundarysize=Dimen%boundarysize
		allocate(Dimpermute%boundary(Dimpermute%boundarysize))
		Dimpermute%boundary(1)=0
		k=1
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
			call error_stop()
		end if
		return
	end subroutine
	logical function if_original_dim(dimen)
		type(Dimension),intent(in)::dimen
		if_original_dim=.true.
		if(dimen%sample_dimension_flag) return
		if(dimen%LenDimData.eq.dimen%Dimsize) return	
		if_original_dim=.false.
		return
	end function
	
	
	subroutine allocateCheckName1(A,lenA)
		type(DimensionName),allocatable,intent(inout)::A(:)
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
	subroutine allocateCheckNameint(A,lenA)
		type(DimensionintName),allocatable,intent(inout)::A(:)
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
	subroutine sent_dimensionName(Name1,Name2,ID1,ID2,ierr,MPIcommon)
		type(DimensionName),intent(in)::Name1
		type(DimensionName),intent(inout)::Name2
		integer,optional,intent(in)::MPIcommon
		integer,intent(in)::ID1,ID2,ierr
		integer::proID,proNum,tag,istatus(MPI_STATUS_SIZE),mpi_comm
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		
		if(ID1.eq.ID2) return !The some cpu, do nothing
		
		if((proID.ne.ID1).and.(proID.ne.ID2)) return!The proID do not sent or recv, return
		
!************************   TensorName   *************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Name1%TensorName,len_of_Name,MPI_CHARACTER,ID2,tag,mpi_comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Name2%TensorName,len_of_Name,MPI_CHARACTER,ID1,tag,mpi_comm,istatus,ierr)
		end if
!*************************   DimenName   **************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Name1%DimenName,len_of_Name,MPI_CHARACTER,ID2,tag,mpi_comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Name2%DimenName,len_of_Name,MPI_CHARACTER,ID1,tag,mpi_comm,istatus,ierr)
		end if
		return
	end subroutine
	
	subroutine BCAST_DimensionName(DimenName,ID,ierr,MPIcommon)
		type(DimensionName),intent(inout)::DimenName
		integer,intent(in)::ID,ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,istatus(MPI_STATUS_SIZE),mpi_comm
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
!************************   TensorName   *************************************************		
		call MPI_BCAST(DimenName%TensorName,len_of_Name,MPI_CHARACTER,ID,mpi_comm,ierr)	
!*************************   DimenName   **************************************************		
		call MPI_BCAST(DimenName%DimenName,len_of_Name,MPI_CHARACTER,ID,mpi_comm,ierr)	
		return
	end subroutine
	
	subroutine sent_dimensionIntName(Name1,Name2,ID1,ID2,ierr,MPIcommon)
		type(dimensionIntName),intent(in)::Name1
		type(dimensionIntName),intent(inout)::Name2
		integer,intent(in)::ID1,ID2,ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,istatus(MPI_STATUS_SIZE),mpi_comm
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		
		if(ID1.eq.ID2) return !The some cpu, do nothing
		
		if((proID.ne.ID1).and.(proID.ne.ID2)) return!The proID do not sent or recv, return
		
!************************   TensorName   *************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Name1%TensorName,len_of_intTensorName,MPI_integer,ID2,tag,mpi_comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Name2%TensorName,len_of_intTensorName,MPI_integer,ID1,tag,mpi_comm,istatus,ierr)
		end if
!*************************   DimenName   **************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Name1%DimenName,len_of_intDimenName,MPI_integer,ID2,tag,mpi_comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Name2%DimenName,len_of_intDimenName,MPI_integer,ID1,tag,mpi_comm,istatus,ierr)
		end if
		return
	end subroutine
	
	subroutine BCAST_dimensionIntName(DimenName,ID,ierr,MPIcommon)
		type(dimensionIntName),intent(inout)::DimenName
		integer,intent(in)::ID,ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,istatus(MPI_STATUS_SIZE),mpi_comm
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
!************************   TensorName   *************************************************		
		call MPI_BCAST(DimenName%TensorName,len_of_intTensorName,MPI_integer,ID,mpi_comm,ierr)	
!*************************   DimenName   **************************************************		
		call MPI_BCAST(DimenName%DimenName,len_of_intDimenName,MPI_integer,ID,mpi_comm,ierr)	
		return
	end subroutine
	
	subroutine sent_Dimension(Dimen1,Dimen2,ID1,ID2,ierr,MPIcommon)
		type(Dimension),intent(in)::Dimen1
		type(Dimension),intent(inout)::Dimen2
		integer,intent(in)::ID1,ID2,ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,lenDimData,istatus(MPI_STATUS_SIZE),i,nameflag,mpi_comm
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		lenDimData=0
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		
		if(ID1.eq.ID2) return !The some cpu, do nothing
		
		if((proID.ne.ID1).and.(proID.ne.ID2)) return!The proID do not sent or recv, return
		
!*********************   LenDimData   ******************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%LenDimData,1,MPI_integer,ID2,tag,mpi_comm,ierr)
			if(Dimen1%LenDimData.eq.0) return !No date,return
			lenDimData=Dimen1%LenDimData
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%LenDimData,1,MPI_integer,ID1,tag,mpi_comm,istatus,ierr)
			if(Dimen2%LenDimData.eq.0) then
				call emptyDimension(Dimen2)
				return !No date,return
			else
				call allocateCheck(Dimen2%DimData,Dimen2%LenDimData)
			end if
			lenDimData=Dimen2%LenDimData
		end if
!******************   DimData  ***********************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%DimData(1:dimen1%lenDimData),dimen1%lenDimData,MPI_integer,ID2,tag,mpi_comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%DimData(1:dimen2%lenDimData),dimen2%lenDimData,MPI_integer,ID1,tag,mpi_comm,istatus,ierr)
		end if		
!******************   NameFlag  ***********************************************************
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%NameFlag,1,MPI_integer,ID2,tag,mpi_comm,ierr)
			nameflag=Dimen1%NameFlag
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%NameFlag,1,MPI_integer,ID1,tag,mpi_comm,istatus,ierr)
			nameflag=Dimen2%NameFlag
			if(Dimen2%NameFlag.eq.1) then
				call allocateCheckName(Dimen2%DimName,dimen2%lenDimData)
			end if
			if(Dimen2%NameFlag.eq.2) then
				call allocateCheckName(Dimen2%DimintName,dimen2%lenDimData)
			end if
		end if
!******************   dimensionName  ***********************************************************
		if(nameflag.eq.1)then
			do i=1,lenDimData
				call sent_dimensionName(Dimen1%DimName(i),Dimen2%DimName(i),ID1,ID2,ierr)
			end do
		else if(nameflag.eq.2)then
			do i=1,lenDimData
				call sent_dimensionIntName(Dimen1%DimIntName(i),Dimen2%DimIntName(i),ID1,ID2,ierr)
			end do
		end if
!******************   sample_dimension_flag  ***********************************************************
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%sample_dimension_flag,1,MPI_logical,ID2,tag,mpi_comm,ierr)
			if(Dimen1%sample_dimension_flag)return
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%sample_dimension_flag,1,MPI_logical,ID1,tag,mpi_comm,istatus,ierr)
			if(Dimen2%sample_dimension_flag)return
		end if		

!***********************    Dimsize   ****************************************************				
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%Dimsize,1,MPI_integer,ID2,tag,mpi_comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%Dimsize,1,MPI_integer,ID1,tag,mpi_comm,istatus,ierr)
		end if
!*********************   boundarysize   ******************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%boundarysize,1,MPI_integer,ID2,tag,mpi_comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%boundarysize,1,MPI_integer,ID1,tag,mpi_comm,istatus,ierr)
			call allocateCheck(Dimen2%boundary,Dimen2%boundarysize)
		end if
!*********************   boundary   *******************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%boundary(1:Dimen1%boundarysize),Dimen1%boundarysize,MPI_integer,ID2,tag,mpi_comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%boundary(1:Dimen2%boundarysize),Dimen2%boundarysize,MPI_integer,ID1,tag,mpi_comm,istatus,ierr)
		end if

		return
	end subroutine
	
	subroutine BCAST_Dimension(Dimen1,ID,ierr,MPIcommon)
		type(Dimension),intent(inout)::Dimen1
		integer,intent(in)::ID,ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,istatus(MPI_STATUS_SIZE),i,mpi_comm
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
!*******************   LenDimData   **********************************************************		
		call MPI_BCAST(Dimen1%LenDimData,1,MPI_integer,ID,mpi_comm,ierr)	
		if(Dimen1%LenDimData.eq.0) then
			call emptyDimension(Dimen1)
			return !No date,return
		end if
		if(proID.ne.ID) then
			call allocateCheck(Dimen1%DimData,Dimen1%LenDimData)
		end if
		call MPI_BCAST(Dimen1%DimData(1:dimen1%lenDimData),Dimen1%LenDimData,MPI_integer,ID,mpi_comm,ierr)	
!******************   NameFlag  ***********************************************************
		call MPI_BCAST(Dimen1%NameFlag,1,MPI_logical,ID,mpi_comm,ierr)	
!******************   dimensionName  ***********************************************************
		if(Dimen1%NameFlag.eq.1)then
			if(proID.ne.ID) then
				call allocateCheckName(Dimen1%DimName,Dimen1%LenDimData)
			end if
			do i=1,Dimen1%LenDimData
				call BCAST_DimensionName(Dimen1%DimName(i),ID,ierr)
			end do
		end if
		if(Dimen1%NameFlag.eq.2)then
			if(proID.ne.ID) then
				call allocateCheckName(Dimen1%DimIntName,Dimen1%LenDimData)
			end if
			do i=1,Dimen1%LenDimData
				call BCAST_DimensionIntName(Dimen1%DimIntName(i),ID,ierr)
			end do
		end if
!******************   sample_dimension_flag  ***********************************************************
		call MPI_BCAST(Dimen1%sample_dimension_flag,1,MPI_logical,ID,mpi_comm,ierr)
		if(Dimen1%sample_dimension_flag) return
!*******************   boundarysize   **********************************************************		
		call MPI_BCAST(Dimen1%boundarysize,1,MPI_integer,ID,mpi_comm,ierr)	
!*********************   boundary   ***************************************************		
		if(proID.ne.ID) then
			call allocateCheck(Dimen1%boundary,Dimen1%boundarysize)
		end if
		call MPI_BCAST(Dimen1%boundary(1:Dimen1%boundarysize),Dimen1%boundarysize,MPI_integer,ID,mpi_comm,ierr)		
!********************   Dimsize   ****************************************************	
		call MPI_BCAST(Dimen1%Dimsize,1,MPI_integer,ID,mpi_comm,ierr)	
		return
	end subroutine
end module
!****************************************************
!*************** ENF OF Dimension *************
!****************************************************























