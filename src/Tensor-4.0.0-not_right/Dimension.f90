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
! Report bugs of the package to sj.dong@outlook.com
!************************************************************
!************* START OF Dimension *******************
!************************************************************


module Dimension_typede
	use Tools
	use mpi
	use memory_type
	implicit none
	private
	type(memory),private::WorkingMemory
	public::Dimension
	type Dimension
		integer,private :: boundarysize=0!The length of the boundary that is useful
		integer,private :: lenDimData=0!The length of the DimData that is useful,If it is 0, means no data, use it to tell if allocate dimension
		integer,private :: Dimsize=0!The Rank of The Tensor
		integer,allocatable,private :: boundary(:)!only boundary(1:boundarysize) are the usefull data
		integer,allocatable,private :: DimData(:)!only DimData(1:lenDimData) are the usefull data
		character(len=len_of_Name),allocatable,private::DimName(:)!only DimName(1:lenDimData) are the usefull data
		logical,private::nameflag=.false.!=.false.  means no name
		                           		 !=.true.  means there are names 
		logical,private::sample_dimension_flag=.true.!if true, only store DimData, DimName, DimIntName, nameflag and lenDimData
		                                             !When do not do the fuse, split and
		                                             !no name, it gose faster.
	contains 
		procedure::Nameorder2,NameorderArray
		generic,public::Nameorder=>Nameorder2,NameorderArray
		procedure::Nameorder2_Check,NameorderArray_Check
		generic,public::FindOrder=>Nameorder2_Check,NameorderArray_Check
		procedure::setDimName0,setDimName1,setDimName2
		generic,public::setName=>setDimName0,setDimName1,setDimName2
		procedure::getName1,outNameAll,outNameTenChar
		procedure,public::outSomedimensionName=>outAllNameTenChar
		generic,public::getName=>getName1,outNameAll,outNameTenChar
		procedure,public::Dimpermute
		procedure::DimName_i,Dim_i,AllDim
		procedure,public::outDimData
		procedure::emptyDimension,cleanDimension
		generic,public::empty => emptyDimension
		generic,public::deallocate=>cleanDimension
		procedure,public::emptyName=>emptyDimensionName
		procedure,public::deallocateName=>cleanDimensionName
		generic,public::dim=>DimName_i,Dim_i,AllDim
		generic,public::subdim=>getSubDim2,getSubDim3,getSubDim2_name,getAllDim!output class(dimension)
		procedure,public::getRank=>getDimSize
		procedure,public::size=>productdim
		procedure,public::outlenDimData
		procedure::Dprint,Dprint2,readdimension,Dprint_File,Dprint2_file
		generic,public::print=>Dprint
		generic,public::print_file=>Dprint_File
		generic,public::info=>Dprint2
		generic,public::info_file=>Dprint2_file
		generic,public::read=>readdimension
		procedure,public::fuse=>fuseDimension
		procedure,public::fuseIndex=>fuseIndexFunc!no longer use
		procedure,public::split =>splitDimension
		procedure::resetDimension1,resetDimension2
		generic,public::resetDim =>resetDimension1,resetDimension2

		procedure::RNDimRoutine,RNDimRoutineint,RNDimRoutinechar
		generic,public::killLeg=>RNDimRoutine,RNDimRoutineint,RNDimRoutinechar
		procedure::RNDimFun,RNDimFunInt,RNDimFunChar
		generic,public::getkillLeg=>RNDimFun,RNDimFunInt,RNDimFunChar
		procedure,public::out_simple_flag=>out_simple_dimension_flag
		procedure,public::if_simple_dimension=>out_simple_dimension_flag
		procedure,public::getSubDimlen
		procedure,public::if_original_dim
		procedure,public::outNameFlag
		procedure,public::getNameFlag=>outNameFlag
		procedure,public::Same_Name_Order_forwards
      	procedure,public::Same_Name_Order_backwards
      	procedure,public::ifName
      	procedure::send_Dimension,BCAST_Dimension
      	procedure,public::MPI_send_Dimension=>send_Dimension
      	generic,public::Send=>send_Dimension
      	procedure,public::MPI_BCAST_Dimension=>BCAST_Dimension
      	generic,public::BCAST=>BCAST_Dimension
		!********* assignment **************************
		procedure::DimInitialization!dimension=vector
		procedure::DimInitialization2
		procedure,pass(Dimen)::copyDimToVec
		generic,public :: assignment(=) => DimInitialization,DimInitialization2,copyDimToVec
		!********* operators **************************
		procedure::fuseDimension_val,fuseDimension_vec
		generic,public :: operator(.fuse.)=>fuseDimension_val,fuseDimension_vec
		procedure::splitDimension2,splitDimension3,splitDimensionAll
		generic,public :: operator(.split.)=>splitDimension2,splitDimension3,splitDimensionAll
		procedure::AddFunc,AddFunc2
		procedure,pass(Dimen)::AddFunc3
		generic,public :: operator(+) =>AddFunc,AddFunc2,AddFunc3
		procedure::getSubDim2,getSubDim3,getSubDim2_name,getAllDim
		generic,public :: operator(.subdim.)=>getSubDim2,getSubDim3,getSubDim2_name,getAllDim
		generic,public :: operator(.dim.)=>DimName_i,Dim_i,AllDim
		procedure::equal_of_dim
		generic,public :: operator(.equ.)=>equal_of_dim
		generic,public :: operator(.eq.)=>equal_of_dim
	end type Dimension

	public::dimpermute_forwards,dimpermute_backwards,dimpermute_forwards_index,dimpermute_backwards_index
	public::writemess,if_long_Name,copydimension
	interface writemess
		module procedure writemess_dimension
	end interface
	public::Dimension_memory_report,Dimension_memory_length,dellocate_Dimension_memory
contains
	subroutine dellocate_Dimension_memory()
		call WorkingMemory%deallocate()
		return
	end subroutine
	
	subroutine Dimension_memory_report()
		call writemess('The memory used in Dimension are:')
		call WorkingMemory%print()
		call writemess(' ')
	end subroutine	
	subroutine Dimension_memory_length(length)
		integer,intent(inout)::length(:)
		call WorkingMemory%getlength(length)
	end subroutine	

	subroutine writemess_dimension(Dimen,cpu_number)!overwrite writemess
		type(dimension),intent(in)::Dimen
		integer,optional,intent(in)::cpu_number
		character(len=max_len_of_char)::w
		CHARACTER*5000,allocatable::ws(:)
		integer,allocatable :: dimenVec(:)
		integer::i,maxdim,totoal
		CHARACTER*100::wlen,w2
		totoal=Dimen%LenDimData
		if(totoal.eq.0)then
			call writemess('There is no data in the dimension',cpu_number)
			return
		end if
		if(Dimen%sample_dimension_flag)then
			w='dimension:('
			do i=1,Dimen%LenDimData-1
				w=w+Dimen%DimData(i)+','
			end do
			w=w+Dimen%DimData(Dimen%LenDimData)+')'
			call writemess(w,cpu_number)
		else
			call copydimension(dimenVec,Dimen)
			w='dimension:('
			do i=1,size(dimenVec)-1
				w=w+dimenVec(i)+','
			end do
			w=w+dimenVec(size(dimenVec))+'),It is not original dimension'
			call writemess(w,cpu_number)
			w='original dimension:('
			do i=1,Dimen%LenDimData-1
				w=w+Dimen%DimData(i)+','
			end do
			w=w+Dimen%DimData(Dimen%LenDimData)+')'
			call writemess(w,cpu_number)
		end if
		if(Dimen%nameflag)then
			allocate(ws(Dimen%lenDimData))
			ws=Dimen%DimName(1:Dimen%lenDimData)
			w="index Name are:"
			do i=1,Dimen%lenDimData
				w=w+','+ws(i)
			end do
			w=w+'.'
			call writemess(w,cpu_number)
		end if
		return
	end subroutine

	function Dimensioniniti(boundarysize,Dimsize,boundary,DimData)
		class(Dimension),allocatable::Dimensioniniti
		integer,intent(in) :: boundarysize
		integer,intent(in) :: Dimsize
		integer,intent(in) :: boundary(:)
		integer,intent(in) :: DimData(:)
		allocate(Dimension::Dimensioniniti)
		allocate(Dimensioniniti%boundary(size(boundary,1)))
		allocate(Dimensioniniti%DimData(size(DimData,1)))
		Dimensioniniti%boundary=boundary
		Dimensioniniti%DimData=DimData
		Dimensioniniti%boundarysize=boundarysize
		Dimensioniniti%Dimsize=Dimsize
		Dimensioniniti%lenDimData=size(DimData,1)
		Dimensioniniti%Nameflag=.false.
		Dimensioniniti%sample_dimension_flag=.false.
		return
	end function

	subroutine resetDimension1(dimen,DimData)
		class(Dimension),intent(inout) ::Dimen
		integer,intent(in) :: DimData(:)
		if(size(DimData).ne.dimen%outlenDimData())then
			call writemess("Can not reset the dimension in type(Dimension)",-1)
			call error_stop()
		end if
		Dimen%DimData(1:Dimen%lenDimData)=DimData
		return
	end subroutine
	subroutine resetDimension2(dimen,DimData)
		class(Dimension),intent(inout)::dimen
		class(Dimension),intent(in)::DimData
		if(dimen%size().ne.DimData%size()) then
			call writemess("ERROR in resetdim2",-1)
			call writemess(dimen%size()+','+DimData%size(),-1)
			call writemess(dimen,-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		dimen=DimData
		return
	end subroutine
	!***********************************************
	!********* assignment **************************
	!***********************************************

	subroutine DimInitialization(Dimen,DimData)
		class(Dimension),intent(inout) ::Dimen
		integer,intent(in) :: DimData(:)
		integer::i
		Dimen%lenDimData=size(DimData)
		call allocateCheck(Dimen%DimData,Dimen%lenDimData)
		Dimen%DimData(1:Dimen%lenDimData)=DimData
		Dimen%sample_dimension_flag=.true.
		Dimen%nameflag=.false.
		return
	end subroutine
	subroutine DimInitialization2(Dimen,Dimen2)
		class(Dimension),intent(inout) ::Dimen
		class(Dimension),intent(in) ::Dimen2
		integer::lenDim,lenBoun,i
		lenDim=Dimen2%lenDimData
		if(lenDim.eq.0) then !no Data
			Dimen%Dimsize=0
			Dimen%boundarysize=0
			Dimen%lenDimData=0
			Dimen%nameflag=.false.
			Dimen%sample_dimension_flag=.true.
			return
		end if
		Dimen%lenDimData=lenDim
		call allocateCheck(Dimen%DimData , lenDim)
		Dimen%DimData(1:lenDim)=Dimen2%DimData(1:lenDim)
		Dimen%nameflag=Dimen2%nameflag
		Dimen%sample_dimension_flag=Dimen2%sample_dimension_flag
		if((.not.Dimen%nameflag).and.Dimen%sample_dimension_flag) return
		
		if(Dimen2%nameflag)then
			call allocateCheck(dimen%dimName,lenDim)
			do i=1,lenDim
				dimen%dimName(i)=dimen2%dimName(i)
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
	subroutine copyDimToVec(dimenVec,Dimen)
		integer,intent(inout) :: dimenVec(:)
		class(Dimension),intent(in) :: Dimen
		integer :: i
		if(Dimen%sample_dimension_flag)then
			if(size(dimenVec).lt.dimen%LenDimData)then
				call writemess("ERROR in assignment dimension to array ",-1)
				call writemess("array(:)=dimension(:),size(array)<size(dimension)",-1)
				call writemess('size(dimenVec)='+size(dimenVec)+',Dimen%LenDimData='+Dimen%LenDimData,-1)
				call error_stop()
			end if
			dimenVec(1:dimen%LenDimData)=dimen%DimData(1:dimen%LenDimData)
			return
		end if
		if(size(dimenVec).lt.Dimen%Dimsize)then
			call writemess("ERROR in assignment dimension to array ",-1)
			call writemess("array(:)=dimension(:),size(array)<size(dimension)",-1)
			call writemess('size(dimenVec)='+size(dimenVec)+',Dimen%Dimsize='+Dimen%Dimsize,-1)
			call error_stop()
		end if
		do i=1,Dimen%Dimsize
			dimenVec(i)=Dim_i(Dimen,i)
		end do
	end subroutine


	!***********************************************
	!********* memory function  ******************
	!***********************************************


	subroutine emptyDimension(Dimen)!make dimension to a a empty dimension.but do not deallocate
		class(Dimension),intent(inout)::	Dimen
		if(deallocate_memory_flag)then
			call cleanDimension(Dimen)
			return
		end if
		Dimen%boundarysize=0
		Dimen%Dimsize=0
		Dimen%lenDimData=0
		Dimen%nameflag=.false.
		Dimen%sample_dimension_flag=.true.
		return
	end subroutine

	subroutine cleanDimension(Dimen)
		class(Dimension),intent(inout)::	Dimen
		Dimen%boundarysize=0
		Dimen%Dimsize=0
		Dimen%lenDimData=0
		Dimen%nameflag=.false.
		if(allocated(Dimen%boundary)) then
			deallocate(Dimen%boundary)
		end if
		if(allocated(Dimen%DimData)) then
			deallocate(Dimen%DimData)
		end if
		if(allocated(Dimen%DimName)) then
			deallocate(Dimen%DimName)
		end if
		Dimen%sample_dimension_flag=.true.
		return
	end subroutine


	!***********************************************
	!********* get dimension data ******************
	!***********************************************

	logical function out_simple_dimension_flag(dimen)
   	class(Dimension),intent(in)::dimen
		out_simple_dimension_flag=dimen%sample_dimension_flag
		return
	end function
	logical function outNameFlag(dimen)
		class(Dimension),intent(in)::dimen
		outNameFlag=dimen%NameFlag
		return
	end function

	!
		!	return the inde  dimension	,outpout in a vector of the dimenison    
		! If do the fuse, onedimenison will have more than one value
		! [2,3,4,5]	-->fuse the 2,3 index -->[2,(3,4),5]!the dimension is 2,12,5
		! then getSubDim(Dimen,2,dimenVec)==>dimenVec=[3,4]

	subroutine getSubDim(Dimen,inde,dimenVec)
		class(Dimension),intent(in) :: Dimen
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

	function getSubDim2(Dimen,inde) 
		class(Dimension),allocatable::getSubDim2
		class(Dimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		integer::boundarysize,boundary(2),Dimsize
		integer,allocatable :: Dimdata(:)
		integer::i,j,Dlen
		allocate(getSubDim2,mold=Dimen)
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
		if(.not.getSubDim2%NameFlag) return
		call getSubDimName(Dimen,inde,getSubDim2%dimName)
		return
	end function
	
	function getSubDim3(Dimen,inde) 
		class(Dimension),allocatable::getSubDim3
		class(Dimension),intent(in) :: Dimen
		integer,intent(in) :: inde(2)
		integer::boundarysize,Dimsize
		integer::i,j,Dlen,ith1,ith2
		allocate(getSubDim3,mold=Dimen)
		if(Dimen%sample_dimension_flag)then
			Dlen=inde(2)-inde(1)+1
			allocate(getSubDim3%DimData(Dlen))
			getSubDim3%DimData=Dimen%DimData(inde(1):inde(2))
			getSubDim3%lenDimData=Dlen
			getSubDim3%sample_dimension_flag=.true.
			getSubDim3%NameFlag=Dimen%NameFlag
			if(.not.getSubDim3%NameFlag) return
			allocate(getSubDim3%DimName(Dlen))
			getSubDim3%DimName=Dimen%DimName(inde(1):inde(2))
			return
		else
			ith1=Dimen%boundary(inde(1))+1
			ith2=Dimen%boundary(inde(2)+1)
			Dlen=ith2-ith1+1
			if((inde(2)-inde(1)+1).eq.Dlen)then
				allocate(getSubDim3%DimData(Dlen))
				getSubDim3%DimData=Dimen%DimData(ith1:ith2)
				getSubDim3%lenDimData=Dlen
				getSubDim3%sample_dimension_flag=.true.
				getSubDim3%NameFlag=Dimen%NameFlag
			else
				getSubDim3%sample_dimension_flag=.false.
				ith1=Dimen%boundary(inde(1))+1
				ith2=Dimen%boundary(inde(2)+1)
				Dlen=ith2-ith1+1
				allocate(getSubDim3%DimData(Dlen))
				getSubDim3%lenDimData=Dlen
				getSubDim3%DimData=Dimen%DimData(ith1:ith2)
				Dlen=inde(2)-inde(1)+2
				allocate(getSubDim3%boundary(Dlen))
				getSubDim3%boundary=Dimen%boundary(inde(1):inde(2)+1)
				getSubDim3%boundary=getSubDim3%boundary-getSubDim3%boundary(1)
				getSubDim3%boundarysize=Dlen
				getSubDim3%Dimsize=getSubDim3%boundarysize-1
				getSubDim3%NameFlag=Dimen%NameFlag
			end if
			if(.not.getSubDim3%NameFlag) return
			Dlen=ith2-ith1+1
			allocate(getSubDim3%DimName(Dlen))
			getSubDim3%DimName=Dimen%DimName(ith1:ith2)
			return
		end if
	end function

	function  getSubDim2_name(Dimen,w)
		class(Dimension),allocatable::getSubDim2_name
		class(Dimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w
		integer::inde
		allocate(getSubDim2_name,mold=Dimen)
		inde=Dimen%NameOrder(w)
		getSubDim2_name=getSubDim2(Dimen,inde)
		return
	end function
	function  getAllDim(Dimen)
		class(Dimension),allocatable::getAllDim
		class(Dimension),intent(in) :: Dimen
		allocate(getAllDim,mold=Dimen)
		getAllDim=Dimen
		return
	end function

	

	!	return  all the  dimension	,outpout in a vector	

	subroutine copydimension(dimenVec,Dimen)
		integer,allocatable,intent(inout) :: dimenVec(:)
		class(Dimension),intent(in) :: Dimen
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
			dimenVec=dimen%DimData(1:dimen%LenDimData)
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
	end subroutine	

	integer function getSubDimlen(Dimen,inde)
		class(Dimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		if(Dimen%sample_dimension_flag)then
			getSubDimlen=1
			return
		end if
		getSubDimlen=Dimen%boundary(inde+1)-Dimen%boundary(inde)
	end function


	subroutine outDimData(dimen,DimData)
		class(Dimension),intent(in)::dimen
		integer,intent(inout)::DimData(:)
		DimData=dimen%DimData(1:dimen%lenDimData)
		return
	end subroutine


	integer function Dim_i(Dimen,inde)
		class(Dimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		integer :: i,D1
		if(Dimen%sample_dimension_flag)then
			if(inde.gt.dimen%LenDimData) then 
				write(*,*) "ERROR in getting the ith dimension"
				write(*,*)"stop"
				call error_stop()
			end if
			Dim_i=Dimen%DimData(inde)
			return
		end if
		if(inde.gt.dimen%Dimsize) then 
			write(*,*) "ERROR in getting the ith dimension"
			write(*,*)"stop"
			call error_stop()
		end if
		D1=1
		do i=Dimen%boundary(inde) + 1,Dimen%boundary(inde+1)
			D1=D1*Dimen%DimData(i)
		end do
		Dim_i=D1
		return
	end function 
	integer function DimName_i(Dimen,w)
		class(Dimension),intent(in) :: Dimen
		character(len=*),intent(in)::w
		integer :: inde
		inde=Dimen%Nameorder(w)
		if(inde.eq.0)then
			DimName_i=0
			return
		end if
		DimName_i=Dim_i(Dimen,inde)
	end function 

	function AllDim(Dimen)
		integer,allocatable::AllDim(:)
		class(Dimension),intent(in) :: Dimen
		allocate(AllDim(Dimen%getRank()))
		call copyDimToVec(AllDim,Dimen)
		return
	end function 

	integer	function getDimSize(Dimen)
		class(Dimension),intent(in) :: Dimen
		if(Dimen%sample_dimension_flag)then
			getDimSize=Dimen%LenDimData
		else
			getDimSize=Dimen%dimSize
		end if
		return
	end function
	integer	function outlenDimData(Dimen)
		class(Dimension),intent(in) :: Dimen
		outlenDimData=Dimen%lenDimData
		return
	end function
	integer function productdim(Dimen)
		class(Dimension),intent(in) :: Dimen
		productdim=product(Dimen%DimData(1:dimen%lenDimData))
		return
	end function	

	!***********************************************
	!********* get dimension name ******************
	!***********************************************

	subroutine emptyDimensionName(Dimen)
		class(Dimension),intent(inout)::	Dimen
		Dimen%nameflag=.false.
		return
	end subroutine

	subroutine cleanDimensionName(Dimen)
		class(Dimension),intent(inout)::	Dimen
		if(allocated(Dimen%DimName))then
			deallocate(Dimen%DimName)
		end if
		Dimen%nameflag=.false.
		return
	end subroutine

	subroutine getSubDimName(Dimen,inde,dimenName)
		class(Dimension),intent(in) :: Dimen
		character(len=len_of_Name),allocatable,intent(inout) :: dimenName(:)
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


	!Nameorder2
		!find the index,whose name is w	,output the order of it in the dimension
		!If can not find , output 0

	function Nameorder2(dimen,w)
		integer::Nameorder2
		class(Dimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w
		integer::i
		character(len=len_of_Name)::nam
		if(.not.Dimen%nameflag)then
			call writemess("There is no character name in the dimension,Nameorder2",-1)
			call writemess('nameflag='+Dimen%nameflag,-1)
			call error_stop()
		end if
		if(.not.if_long_Name(w))then
			call writemess('input error, one should input name written as A'+indexsymbol+'B')
			call error_stop
		end if
		if(.not.if_original_dim(dimen))then
			call writemess("Nameorder2 is use in original dimension",-1)
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

	!Nameorder2_check
		!find the index,whose name is w	,output the order of it in the dimension
		!If can not find , error message

	function Nameorder2_check(dimen,w)
		integer::Nameorder2_check
		class(Dimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w
		integer::i
		character(len=len_of_Name)::nam
		if(.not.Dimen%nameflag)then
			call writemess("There is no character name in the dimension,Nameorder2",-1)
			call writemess('nameflag='+Dimen%nameflag,-1)
			call error_stop()
		end if
		if(.not.if_long_Name(w))then
			call writemess('input error, one should input name written as A'+indexsymbol+'B')
			call error_stop
		end if
		if(.not.if_original_dim(dimen))then
			call writemess("Nameorder2 is use in original dimension",-1)
			call error_stop()
		end if
		nam=w
		do i=1,dimen%lenDimData
			if(dimen%dimname(i).equ.nam)then
				Nameorder2_check=i
				return
			end if
		end do
		call writemess('Can Not Find the name:'+w)
		call dimen%print
		call error_stop
		return
	end function

	function NameorderArray(dimen,w)
		class(Dimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w(:)
		integer::NameorderArray(size(w))
		integer::i,lenn
		if(.not.Dimen%nameflag)then
			call writemess("There is no char name in the dimension,NameorderArray",-1)
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			call writemess("NameorderArray is use in original dimension",-1)
			call error_stop()
		end if
		lenn=size(w)
		if(lenn.gt.dimen%lenDimData)then
			call writemess("ERROR i NameorderArray",-1)
			call writemess('length='+lenn+', dimen%Dimsize='+dimen%Dimsize,-1)
			call error_stop()
		end if
		do i=1,lenn
			NameorderArray(i)=Nameorder2(dimen,w(i))
		end do
		return
	end function


	function NameorderArray_check(dimen,w)result(NameorderArray)
		class(Dimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w(:)
		integer::NameorderArray(size(w))
		integer::i,lenn
		if(.not.Dimen%nameflag)then
			call writemess("There is no char name in the dimension,NameorderArray",-1)
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			call writemess("NameorderArray is use in original dimension",-1)
			call error_stop()
		end if
		lenn=size(w)
		if(lenn.gt.dimen%lenDimData)then
			call writemess("ERROR i NameorderArray",-1)
			call writemess('length='+lenn+', dimen%Dimsize='+dimen%Dimsize,-1)
			call error_stop()
		end if
		do i=1,lenn
			NameorderArray(i)=Nameorder2_check(dimen,w(i))
		end do
		return
	end function

	logical function ifName(dimen,w,ch_)
		class(Dimension),intent(in) :: Dimen
		character(len=*),intent(in)::w
		character(len=*),intent(in),optional::ch_
		character(len=10)::ch
		integer::i
		if(.not.Dimen%nameflag)then
			if(dimen%lenDimData.eq.0)then
				call writemess("There is no data in the dimension",-1)
			end if
			call writemess("There is no CHARACTER name in the dimension",-1)
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			call writemess("NameID is use in original dimension",-1)
			call error_stop()
		end if
		if(present(ch_))then
			ch=ch_
		else
			ch='dimension'
		end if
		ifName=.true.
		if(if_long_Name(w))then
			do i=1,dimen%getrank()
				if(w.equ.dimen%DimName(i))return
			end do
		else
			if(ch.equ.'Tensor') then
				do i=1,dimen%getrank()
					if(w.equ.CharacterOutTensorName(Dimen%DimName(i)))return
				end do
			else
				do i=1,dimen%getrank()
					if(w.equ.CharacterOutDimensionName(Dimen%DimName(i)))return
				end do
			end if
		end if
		ifName=.false.
		return
	end function

	!***********************************************
	!********* print the data in dimension *********
	!***********************************************

	subroutine Dprint_File(Dimen,uni)
		class(Dimension),intent(in) ::Dimen
		integer,intent(in)::uni
		integer,allocatable :: dimenVec(:)
		CHARACTER*5000,allocatable::w(:)
		integer::i,maxdim
		CHARACTER*100::wlen,w2
		CHARACTER*5000::words
		write(uni,*) "***   Dimension Data    ***"
		if(Dimen%LenDimData.le.0)then
			write(uni,*) "***   There is no data in Dimension    ***"
			return
		end if
		if(Dimen%sample_dimension_flag)then
			maxdim=maxval(Dimen%DimData(1:Dimen%LenDimData))
			call outputform(w2,maxdim)
			wlen='('+(Dimen%LenDimData)+w2+')'
			write(uni,wlen)Dimen%DimData(1:Dimen%LenDimData)
		else
			call copydimension(dimenVec,Dimen)
			maxdim=maxval(dimenVec)
			call outputform(w2,maxdim)
			wlen='('+size(dimenVec)+w2+')'
			write(uni,wlen)	dimenVec
			write(uni,*)"It is not original dimension "
		end if
		write(uni,*) "***   Dimension END   ***"
		if(Dimen%nameflag)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName(1:Dimen%lenDimData)
			write(uni,*)"index Name are"
			write(uni,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		write(uni,*) " "
		return
	end subroutine

	subroutine Dprint(Dimen)
		class(Dimension),intent(in) ::Dimen
		integer,allocatable :: dimenVec(:)
		CHARACTER*5000,allocatable::w(:)
		integer::i,maxdim
		CHARACTER*100::wlen,w2
		CHARACTER*5000::words
		call writemess("***   Dimension Data    ***",-1)
		if(Dimen%LenDimData.le.0)then
			call writemess( "***   There is no data in Dimension    ***",-1)
			return
		end if
		if(Dimen%sample_dimension_flag)then
			words=' '+Dimen%DimData(1)
			do i=2,Dimen%LenDimData
				words=words+(' ,'+(' '+Dimen%DimData(i)))
			end do
			call writemess(words,-1)
		else
			call copydimension(dimenVec,Dimen)
			words=' '+dimenVec(1)
			do i=2,size(dimenVec)
				words=words+(' ,'+(' '+dimenVec(i)))
			end do
			call writemess(words,-1)
			call writemess("It is not original dimension ",-1)
		end if
		call writemess( "***   Dimension END   ***",-1)
		if(Dimen%nameflag)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName(1:Dimen%lenDimData)
			words=' '+w(1)
			do i=2,Dimen%LenDimData
				words=words+(' ,'+(' '+w(i)))
			end do
			call writemess("index Name are",-1)
			call writemess(words,-1)
		end if
		return
	end subroutine
	subroutine outputform(w,maxdim)
		integer,intent(in)::maxdim
		CHARACTER*100,intent(inout)::w
		integer::num
		num=log(real(maxdim))/log(10.)
		w='I'+(num+3)
		return
	end subroutine

	subroutine Dprint2(Dimen)
		class(Dimension),intent(in) ::Dimen
		integer,allocatable :: dimenVec(:)
		CHARACTER*500,allocatable::w(:)
		integer::i
		write(*,*) "***   Dimension Data    ***"
		write(*,*)"sample_flag ",Dimen%sample_dimension_flag
		write(*,*)"leng ",Dimen%LenDimData
		if(Dimen%LenDimData.le.0)then
			write(*,*) "***   There is no data in Dimension    ***"
			return
		end if
		write(*,*)Dimen%DimData(1:Dimen%LenDimData)
		if(.not.Dimen%sample_dimension_flag)then
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
		write(*,*)"nameflag ",Dimen%nameflag
		if(Dimen%nameflag)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName(1:Dimen%lenDimData)
			write(*,*)"index Name are"
			write(*,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		write(*,*) "***   END   ***"
	end subroutine
	subroutine Dprint2_file(Dimen,uni)
		class(Dimension),intent(in) ::Dimen
		integer,intent(in)::uni
		integer,allocatable :: dimenVec(:)
		CHARACTER*500,allocatable::w(:)
		integer::i
		write(uni,*) "***   Dimension Data    ***"
		write(uni,*)"sample_flag ",Dimen%sample_dimension_flag
		write(uni,*)"leng ",Dimen%LenDimData
		if(Dimen%LenDimData.le.0)then
			write(uni,*) "***   There is no data in Dimension    ***"
			return
		end if
		write(uni,*)Dimen%DimData(1:Dimen%LenDimData)
		if(.not.Dimen%sample_dimension_flag)then
			write(uni,*) "Dimsize"
			write(uni,*) Dimen%Dimsize
			write(uni,*) "boundarysize"
			write(uni,*) Dimen%boundarysize
			write(uni,*) "boundary"
			write(uni,*) Dimen%boundary(1:Dimen%boundarysize)
			call copydimension(dimenVec,Dimen)
			write(uni,*) "		Dimension   	"
			write(uni,*)	dimenVec
		end if
		write(uni,*)"nameflag ",Dimen%nameflag
		if(Dimen%nameflag)then
			allocate(w(Dimen%lenDimData))
			w=Dimen%DimName(1:Dimen%lenDimData)
			write(uni,*)"index Name are"
			write(uni,*)(trim(adjustl(w(i)))//"  ",i=1,Dimen%lenDimData)
		end if
		write(uni,*) "***   END   ***"
	end subroutine
	subroutine readdimension(dimen,uni)
		class(dimension),target,intent(inout)::dimen
		integer,intent(in)::uni
		logical::flag
		character*50::notused
		integer::lendimData,i,boundarysize,nameflag,Dimsize
		integer,allocatable::DimenData(:),boundary(:),dimnameint(:)
		character*50,allocatable::dimname(:)
		
		read(uni,*) notused
		read(uni,*)notused,flag
		read(uni,*)notused,lendimData
		if(lendimData.le.0)then
				read(uni,*) notused
				return
			end if
		allocate(DimenData(lendimData))
		read(uni,*)(DimenData(i),i=1,lendimData)
		if(.not.flag)then
			read(uni,*) notused
			read(uni,*)Dimsize
			read(uni,*) notused
			read(uni,*) boundarysize
			read(uni,*) notused
			allocate(boundary(boundarysize))
			read(uni,*)(boundary(i),i=1,boundarysize)
			read(uni,*) notused
			read(uni,*)	notused
			dimen=Dimensioniniti(boundarysize,Dimsize,boundary,DimenData)
		else
			dimen=DimenData
		end if
		read(uni,*)notused,nameflag
		if(nameflag.eq.1)then
			allocate(dimname(lendimData))
			read(uni,*)notused
			read(uni,*)(dimname(i),i=1,lendimData)
			do i=1,lendimData
				call dimen%setName(i,dimname(i),.true.)
			end do
		end if
		if(nameflag.eq.2)then
				call writemess("cannot read int name",-1)
				call error_stop()
		end if
		read(uni,*) notused
		return
	end subroutine

	!***********************************************
	!********* function for dimension name *********
	!***********************************************

	subroutine check_same_name_in_dimension(dimen)
		class(Dimension),intent(in) :: Dimen
		integer::i,j,rank
		character(len=len_of_Name)::na
		rank=dimen%getRank()
		do i=1,rank
			na=dimen%dimname(i)
			do j=i+1,rank
				if(na.equ.dimen%dimname(j))then
					call writemess('There are two legs with a same name',-1)
					call writemess('The name in the dimension can not be the same',-1)
					call writemess('The names are',-1)
					call writemess(na)
					call writemess(dimen%dimname(j))
					call dimen%print()
					open(unit=1234,file='_ERROR_DIMENSION'+output_ProID+'.err',status='replace')
					call dimen%print_file(1234)
					write(1234,*)na
					close(unit=1234)
					call error_stop
				end if
			end do
		end do
	end subroutine

	function CharacterOutTensorName(w)result(Res)
		character(len=len_of_Name)::Res
		character(len=len_of_Name),intent(in)::w
		Res=trim(adjustl(w)).subl.indexsymbol
		return
	end function
	function CharacterOutDimensionName(w)result(Res)
		character(len=len_of_Name)::Res
		character(len=len_of_Name),intent(in)::w
		Res=trim(adjustl(w)).subr.indexsymbol
		return
	end function
	function changeTensorName(NewTensorName,oldFullName)result(Res)
		character(len=len_of_Name)::Res
		character(len=len_of_Name),intent(in)::NewTensorName,oldFullName
		Res=NewTensorName+indexsymbol+CharacterOutDimensionName(oldFullName)
		return
	end function


	subroutine setDimName0(dimen,w)
		class(Dimension),intent(inout)::dimen
		character(len=*),intent(in)::w
		character(len=len_of_Name)::TensorName,DimenName
		integer::lenDim,i
		logical::flag
		flag=if_long_Name(w)
		if(flag)then
			call writemess("ERROR in Set the name to the Tensor",-1)
			call writemess("input:"+w,-1)
			call writemess('could not contains the indexsymbol:"'+indexsymbol+'"',-1)
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			write(*,*)"The dimension should be in its original dimension"
			call error_stop()
		end if
		lenDim=dimen%LenDimData
		call allocateCheck(dimen%dimname,lenDim)
		if(dimen%nameflag)then!There are name in the dimension,rename
			do i=1,lenDim
				dimen%dimname(i)=changeTensorName(w,dimen%dimname(i))
			end do
		else
			do i=1,lenDim
				dimen%DimName(i)=w+indexsymbol+i
			end do
			dimen%nameflag=.true.
		end if
		if(check_same_name_Flag)call check_same_name_in_dimension(dimen)
		return
	end subroutine

	subroutine setDimName1(dimen,oldname,newname)
		class(Dimension),intent(inout) :: dimen
		CHARACTER(len=*),intent(in)::oldname,newname
		integer::i
		CHARACTER(len=len_of_Name)::OldTname
		logical::oldfullname,newfullname
		if(.not.Dimen%nameflag)then
			call writemess("There is no char name in the dimension,resetindexname",-1)
			call error_stop()
		end if
		oldfullname=if_long_Name(oldname)
		newfullname=if_long_Name(newname)
		if(oldfullname.neqv.newfullname)then
			call writemess("ERROR in setName",-1)
			call error_stop()
		end if
		if(newfullname)then
			do i=1,Dimen%lenDimData
				if(dimen%dimname(i).equ.oldname)then
					dimen%dimname(i)=newname
				end if
			end do
		else
			do i=1,Dimen%lenDimData
				OldTname=CharacterOutTensorName(dimen%dimname(i))
				if(OldTname.equ.newname)then
					dimen%dimname(i)=changeTensorName(newname,dimen%dimname(i))
				end if
			end do
		end if
		if(check_same_name_Flag)call check_same_name_in_dimension(dimen)
		return
	end subroutine

	subroutine setDimName2(dimen,ith,newname,no_check)
		class(Dimension),intent(inout) :: dimen
		CHARACTER(len=*),intent(in)::newname
		integer,intent(in)::ith
		logical,optional,intent(in)::no_check
		CHARACTER(len=len_of_Name)::newTensorName
		logical::fullname
		integer::i,lenName
		if(ith.gt.dimen%LenDimData)then
			write(*,*)"The index is larger than the size of the name"
			write(*,*)ith,dimen%LenDimData
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			if(.not.present(no_check))then
				write(*,*)"setName should be use in original dimension(try call T%split() and setName)"
				call error_stop()
			end if
		end if
		fullname=if_long_Name(newname)
		if(Dimen%nameflag)then!There is no name in dimension
			if(fullname)then
				dimen%dimname(ith)=newname
			else
				dimen%dimname(ith)=changeTensorName(newname,dimen%dimname(ith))
			end if
			if(check_same_name_Flag)call check_same_name_in_dimension(dimen)
			return
		end if
		lenName=dimen%LenDimData
		if(fullname)then
			newTensorName=CharacterOutTensorName(newname)+indexsymbol
		else
			newTensorName=newname+indexsymbol
		end if
		call allocateCheck(dimen%dimname,lenName)
		do i=1,lenName
			dimen%dimname(i)=newTensorName+i
		end do
		if(fullname)then
			dimen%dimname(ith)=newname
		end if
		Dimen%nameflag=.true.
		if(check_same_name_Flag)call check_same_name_in_dimension(dimen)
		return
	end subroutine

	function getName1(dimen,ith)result(outName)
		CHARACTER(len=len_of_Name)::outName
		class(Dimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(.not.Dimen%nameflag)then
			call writemess("There is no name in the dimension",-1)
			call error_stop()
		end if
		if(ith.gt.size(Dimen%DimName))then
			call writemess("The index is larger than the size of the name",-1)
			call writemess("ith="+ith+", size(Dimen%DimName)="+size(Dimen%DimName),-1)
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			call writemess("indexID is use in original dimension",-1)
			call error_stop()
		end if
		outName=Dimen%DimName(ith)
		return
	end function
	function outNameAll(Dimen)
		CHARACTER(len=characterlen),allocatable::outNameAll(:)
		class(Dimension),intent(in) :: Dimen
		integer::i,length
		length=Dimen%Getrank()
		allocate(outNameAll(length))
		do i=1,length
			outNameAll(i)=getName1(dimen,i)
		end do
		return
	end function

	function outNameTenChar(dimen,w,ch)
		CHARACTER(len=len_of_Name)::outNameTenChar
		class(Dimension),intent(in) :: Dimen
		character(len=*),intent(in)::w
		character(len=*),intent(in),optional::ch
		CHARACTER(len=len_of_Name)::dimenname
		integer::i,check
		if(.not.Dimen%nameflag)then
			if(dimen%lenDimData.eq.0)then
				call writemess("There is no data in the dimension",-1)
			end if
			call writemess("There is no CHARACTER name in the dimension",-1)
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			call writemess("NameID is use in original dimension",-1)
			call error_stop()
		end if
		check=0
		do i=1,dimen%getrank()
			dimenname=CharacterOutDimensionName(Dimen%DimName(i))
			if(w.equ.dimenname)then
				outNameTenChar=Dimen%DimName(i)
				check=check+1
			end if
		end do
		if(check.eq.0)then
			call writemess('Can Not find the name in the dimension',-1)
			call writemess('name='+w)
			call dimen%print()
			call error_stop()
		end if
		if(check.gt.1)then
			call writemess('there more than 2 legs have the same dim name in the dimension',-1)
			call writemess('name='+w)
			call dimen%print()
			call error_stop()
		end if
		if(present(ch))then
			if(ch.equ.'Tensor')outNameTenChar=outNameTenChar.subl.indexsymbol
		end if
		return
	end function

	logical function outAllNameTenChar(dimen,outchar,w)
		class(Dimension),intent(in) :: Dimen
		character(len=*),intent(in)::w
		character(len=max_len_of_char_in_TData),allocatable,intent(inout)::outchar(:)
		character(len=max_len_of_char_in_TData),allocatable::temp(:)
		integer::i,k
		if(.not.Dimen%nameflag)then
			if(dimen%lenDimData.eq.0)then
				call writemess("There is no data in the dimension",-1)
			end if
			call writemess("There is no CHARACTER name in the dimension",-1)
			call error_stop()
		end if
		if(.not.if_original_dim(dimen))then
			call writemess("NameID is use in original dimension",-1)
			call error_stop()
		end if
		allocate(temp(dimen%getRank()))
		k=0
		do i=1,dimen%getRank()
			if(w.equ.(dimen%DimName(i).subr.indexsymbol))then
				k=k+1
				temp(k)=Dimen%DimName(i)
			end if
		end do
		if(k.eq.0)then
			outAllNameTenChar=.false.
			return
		end if
		outAllNameTenChar=.true.
		allocate(outchar(k))
		outchar=temp(1:k)
		return
	end function


	!***********************************************
	!********* trim dimension              *********
	!***********************************************

	!RNDimFun
		!input dimension [1,1,2,1,3,1,1,4,1]
		!output dimenison [2,3,4]
		! if input [1,1,1,1,1]
		!  ouput [1] without name

	function RNDimFun(dimen) Result(RNDim)
		class(Dimension),allocatable::RNDim
		class(Dimension),intent(in) :: Dimen
		integer::i,lenDimen,lenNewD
		integer,allocatable::Dimindex(:)
		allocate(RNDim,mold=Dimen)
		if(.not.if_original_dim(dimen)) then
			write(*,*)"ERROR IN RNDim,in Diemnsion.f90"
			call error_stop()
		end if
		lenDimen=dimen%LenDimData
		allocate(Dimindex(lenDimen))
		lenNewD=0
		do i=1,lenDimen
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
		RNDim%LenDimData=lenNewD
		RNDim%sample_dimension_flag=Dimen%sample_dimension_flag
		if(.not.Dimen%sample_dimension_flag)then
			RNDim%boundarysize=RNDim%Dimsize+1
			allocate(RNDim%boundary(RNDim%boundarysize))
			do i=1,RNDim%boundarysize
				RNDim%boundary(i)=i-1
			end do
		end if
		RNDim%NameFlag=dimen%NameFlag
		if(dimen%NameFlag)then
			allocate(RNDim%DimName(lenNewD))
			do i=1,lenNewD
				RNDim%DimName(i)=dimen%DimName(Dimindex(i))
			end do
		end if
		return
	end function

	subroutine RNDimRoutine(Dimen)
		class(Dimension),intent(inout) :: Dimen
		Dimen=RNDimFun(dimen)
		return
	end subroutine

	function RNDimFunint(dimen,notkillleg,killFlag) Result(RNDim)
		class(Dimension),allocatable::RNDim
		class(Dimension),intent(in) :: Dimen
		integer,intent(in)::notkillleg
		character(len=*),intent(in),optional::killFlag
		integer::i,lenDimen,lenNewD
		integer,allocatable::Dimindex(:)
		allocate(RNDim,mold=Dimen)
		if(.not.if_original_dim(dimen)) then
			write(*,*)"ERROR IN RNDim,in Diemnsion.f90"
			call error_stop()
		end if
		lenDimen=dimen%LenDimData
		allocate(Dimindex(lenDimen))
		lenNewD=0
		if(present(killFlag))then
			if(killFlag.equ.'kill')then
				do i=1,lenDimen
					if(i.ne.notkillleg)then
						lenNewD=lenNewD+1
						Dimindex(lenNewD)=i
					end if
				end do
			else
				do i=1,lenDimen
					if((Dimen%Dimdata(i).ne.1).or.(i.eq.notkillleg))then
						lenNewD=lenNewD+1
						Dimindex(lenNewD)=i
					end if
				end do
			end if
		else
			do i=1,lenDimen
				if((Dimen%Dimdata(i).ne.1).or.(i.eq.notkillleg))then
					lenNewD=lenNewD+1
					Dimindex(lenNewD)=i
				end if
			end do
		end if
		if(lenNewD.eq.0)then
			RNDim=(/1/)
			return
		end if
			
		allocate(RNDim%DimData(lenNewD))
		do i=1,lenNewD
			RNDim%DimData(i)=dimen%DimData(Dimindex(i))
		end do
		RNDim%Dimsize=lenNewD
		RNDim%LenDimData=lenNewD
		RNDim%sample_dimension_flag=Dimen%sample_dimension_flag
		RNDim%NameFlag=dimen%NameFlag
		if(dimen%NameFlag)then
			allocate(RNDim%DimName(lenNewD))
			do i=1,lenNewD
				RNDim%DimName(i)=dimen%DimName(Dimindex(i))
			end do
		end if
		return
	end function
	function RNDimFunchar(dimen,notkillleg,killFlag) Result(RNDim)
		class(Dimension),allocatable::RNDim
		class(Dimension),intent(in) :: Dimen
		character(len=*),intent(in)::notkillleg
		character(len=*),intent(in),optional::killFlag
		integer::i,lenDimen,lenNewD
		integer,allocatable::Dimindex(:)
		allocate(RNDim,mold=Dimen)
		if(.not.if_original_dim(dimen)) then
			write(*,*)"ERROR IN RNDim,in Diemnsion.f90"
			call error_stop()
		end if
		lenDimen=dimen%LenDimData
		allocate(Dimindex(lenDimen))
		lenNewD=0
		if(present(killFlag))then
			if(killFlag.equ.'kill')then
				do i=1,lenDimen
					if(Dimen%getName(i).nequ.notkillleg)then
						lenNewD=lenNewD+1
						Dimindex(lenNewD)=i
					end if
				end do
			else
				do i=1,lenDimen
					if((Dimen%Dimdata(i).ne.1).or.(Dimen%getName(i).equ.notkillleg))then
						lenNewD=lenNewD+1
						Dimindex(lenNewD)=i
					end if
				end do
			end if
		else
			do i=1,lenDimen
				if((Dimen%Dimdata(i).ne.1).or.(Dimen%getName(i).equ.notkillleg))then
					lenNewD=lenNewD+1
					Dimindex(lenNewD)=i
				end if
			end do
		end if
		if(lenNewD.eq.0)then
			RNDim=(/1/)
			return
		end if
			
		allocate(RNDim%DimData(lenNewD))
		do i=1,lenNewD
			RNDim%DimData(i)=dimen%DimData(Dimindex(i))
		end do
		RNDim%Dimsize=lenNewD
		RNDim%LenDimData=lenNewD
		RNDim%sample_dimension_flag=Dimen%sample_dimension_flag
		RNDim%NameFlag=dimen%NameFlag
		if(dimen%NameFlag)then
			allocate(RNDim%DimName(lenNewD))
			do i=1,lenNewD
				RNDim%DimName(i)=dimen%DimName(Dimindex(i))
			end do
		end if
		return
	end function

	subroutine RNDimRoutineint(Dimen,notkillleg,killFlag)
		class(Dimension),intent(inout) :: Dimen
		integer,intent(in)::notkillleg
		character(len=*),intent(in),optional::killFlag
		Dimen=RNDimFunint(dimen,notkillleg,killFlag)
		return
	end subroutine
	subroutine RNDimRoutinechar(Dimen,notkillleg,killFlag)
		class(Dimension),intent(inout) :: Dimen
		character(len=*),intent(in)::notkillleg
		character(len=*),intent(in),optional::killFlag
		Dimen=RNDimFunchar(dimen,notkillleg,killFlag)
		return
	end subroutine


	!***********************************************
	!********* fuse and split dimension    *********
	!***********************************************

	subroutine outDimension_boundry(dimen1,boundary,boundarysize,Dimsize)
	   	class(Dimension),intent(in)::dimen1
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

	function fuseIndexFunc(dimen,inde,num)result(Res)
		class(Dimension),allocatable::Res
		class(Dimension)::Dimen
		integer,intent(in):: inde,num
		integer :: i,boundarysize,Dimsize
		integer,allocatable::boundary(:)
		allocate(Res,mold=dimen)
		if(Dimen%size().le.0)then
			call Res%empty()
			return
		end if
		if(dimen%sample_dimension_flag)then
			if(inde.gt.dimen%LenDimData) then 
				call writemess("ERROR in fusedimension",-1)
				call error_stop()
			end if
			if((num.le.0).or.(inde.eq.dimen%LenDimData)) then
				Res=dimen
				return
			end if
		else
			if(inde.gt.dimen%Dimsize) then 
				call writemess("ERROR in fusedimension",-1)
				call error_stop()
			end if
			if((num.le.0).or.(inde.eq.dimen%Dimsize)) then
				Res=dimen
				return
			end if
		end if
		call outDimension_boundry(dimen,boundary,boundarysize,Dimsize)
		Res%sample_dimension_flag=.false.
		if(inde+num.ge.boundarysize) then
			Res%boundarysize=inde+1
			Res%Dimsize=inde
			allocate(Res%boundary(Res%boundarysize))
			Res%boundary(:inde)=boundary(:inde)
			Res%boundary(inde+1)=boundary(boundarysize)
			Res%LenDimData=dimen%LenDimData
			allocate(Res%DimData(Res%LenDimData))
			Res%DimData=dimen%DimData(1:dimen%LenDimData)
			Res%Nameflag=dimen%Nameflag
			if(.not.dimen%Nameflag) return
			allocate(Res%DimName(Dimen%lenDimData))
			Res%DimName=dimen%DimName(1:Dimen%lenDimData)
			return
		end if
		Res%boundarysize=boundarysize-num
		Res%Dimsize=Dimsize-num
		allocate(Res%boundary(Res%boundarysize))
		Res%LenDimData=dimen%LenDimData
		allocate(Res%DimData(dimen%LenDimData))
		Res%boundary(:inde)=boundary(:inde)
		do i=inde+1,Res%boundarysize
			Res%boundary(i)=boundary(i+num)
		end do
		Res%DimData=dimen%DimData(1:dimen%LenDimData)
		Res%Nameflag=dimen%Nameflag
		if(dimen%Nameflag) return
		allocate(Res%DimName(Dimen%lenDimData))
		Res%DimName=dimen%DimName(1:Dimen%lenDimData)
		return
	end function


	function fusedimensionFunc(dimen,index1,index2)result(Res)
		class(Dimension),allocatable::Res
		class(Dimension)::Dimen
		integer,intent(in)::index1
		integer,intent(in),optional::index2
		integer:: inde,num
		integer :: i,boundarysize,Dimsize
		integer,allocatable::boundary(:)
		allocate(Res,mold=dimen)
		if(present(index2))then
			num=index2-index1
		else
			num=1
		end if
		if(Dimen%size().le.0)then
			call Res%empty()
			return
		end if
		if(dimen%sample_dimension_flag)then
			if(inde.gt.dimen%LenDimData) then 
				call writemess("ERROR in fusedimension",-1)
				call error_stop()
			end if
			if((num.le.0).or.(inde.eq.dimen%LenDimData)) then
				Res=dimen
				return
			end if
		else
			if(inde.gt.dimen%Dimsize) then 
				call writemess("ERROR in fusedimension",-1)
				call error_stop()
			end if
			if((num.le.0).or.(inde.eq.dimen%Dimsize)) then
				Res=dimen
				return
			end if
		end if
		call outDimension_boundry(dimen,boundary,boundarysize,Dimsize)
		Res%sample_dimension_flag=.false.
		if(inde+num.ge.boundarysize) then
			Res%boundarysize=inde+1
			Res%Dimsize=inde
			allocate(Res%boundary(Res%boundarysize))
			Res%boundary(:inde)=boundary(:inde)
			Res%boundary(inde+1)=boundary(boundarysize)
			Res%LenDimData=dimen%LenDimData
			allocate(Res%DimData(Res%LenDimData))
			Res%DimData=dimen%DimData(1:dimen%LenDimData)
			Res%Nameflag=dimen%Nameflag
			if(.not.dimen%Nameflag) return
			allocate(Res%DimName(Dimen%lenDimData))
			Res%DimName=dimen%DimName(1:Dimen%lenDimData)
			return
		end if
		Res%boundarysize=boundarysize-num
		Res%Dimsize=Dimsize-num
		allocate(Res%boundary(Res%boundarysize))
		Res%LenDimData=dimen%LenDimData
		allocate(Res%DimData(dimen%LenDimData))
		Res%boundary(:inde)=boundary(:inde)
		do i=inde+1,Res%boundarysize
			Res%boundary(i)=boundary(i+num)
		end do
		Res%DimData=dimen%DimData(1:dimen%LenDimData)
		Res%Nameflag=dimen%Nameflag
		if(dimen%Nameflag) return
		allocate(Res%DimName(Dimen%lenDimData))
		Res%DimName=dimen%DimName(1:Dimen%lenDimData)
		return
	end function

	!
		!combine two index of the Tensor,which is index1 index1+1,..,index2-1,index2
		!if index2 larger than rank,the function will contract the index of index1 -> rank


	function fuseDimension_val(Dimen,con_index)result(fuseDim)
		class(Dimension),allocatable::fuseDim
		class(Dimension),intent(in) :: Dimen
		integer,intent(in) :: con_index
		allocate(fuseDim,mold=Dimen)
		fuseDim=fusedimensionFunc(dimen,con_index,1)
		return
	end function
	function fuseDimension_vec(dimen,vector)result(fuseDim)
		class(Dimension),allocatable::fuseDim
		integer,intent(in) ::vector(2)
		class(Dimension),intent(in) :: dimen
		integer ::num
		num=vector(2)-vector(1)
		allocate(fuseDim,mold=Dimen)
		fuseDim=fusedimensionFunc(dimen,vector(1),num)
		return
	end function	
	subroutine fuseDimension(dimen,index1,index2)
		class(Dimension),intent(inout) :: dimen
		integer,intent(in) :: index1
		integer,optional,intent(in)::index2
		integer ::num
		if(present(index2))then
			num=index2-index1
		else
			num=1
		end if
		dimen=fusedimensionFunc(dimen,index1,num)
		return
	end subroutine	

	!
		!(inde,inde+1,..midindde,midindde+1,...rank)-->(inde,inde+1,..midindde),(midindde+1,...rank)	     		
		!if (inde,inde+1 ..midindde),	midindde larger then next elmemt,do nothing
		!for example:	
		!D=[2,2,(3,4,5),2,(3,4)],	boundary=[0,1,2,5,6,8]
		!SplitDimensionFunc(D,3,2)=[2,2,(3,4),5,2,(3,4)],boundary=[0,1,2,4,5,6,8]
		!SplitDimensionFunc(D,3,4) will do nothing

	function SplitDimensionFunc(dimen,inde,midindde)
		class(Dimension),allocatable::SplitDimensionFunc
		class(Dimension),intent(in) :: Dimen
		integer,intent(in):: inde,midindde
		integer :: i
		allocate(SplitDimensionFunc,mold=Dimen)
		if(Dimen%outlenDimdata().le.0)then
			call SplitDimensionFunc%empty()
			return
		end if
		if(dimen%sample_dimension_flag)then
			call writemess("ERROR in SplitDimensionFunc,1",-1)
			call error_stop()
		end if
		if(inde.gt.dimen%DimSize) then
			call writemess("ERROR in SplitDimensionFunc",-1)
			call error_stop()
		end if
		if(dimen%boundary(inde) +midindde .ge.dimen%boundary(inde+1)) then
			SplitDimensionFunc=dimen
			return
		end if
		SplitDimensionFunc%LenDimData=dimen%LenDimData
		allocate(SplitDimensionFunc%DimData(dimen%LenDimData))
		SplitDimensionFunc%DimData=dimen%DimData(1:dimen%LenDimData)
		SplitDimensionFunc%Dimsize=dimen%Dimsize+1
		if(SplitDimensionFunc%LenDimData.eq.SplitDimensionFunc%Dimsize) then !It is sample dimension
			SplitDimensionFunc%sample_dimension_flag=.true.
			SplitDimensionFunc%Dimsize=0
			SplitDimensionFunc%boundarysize=0
		else
			SplitDimensionFunc%sample_dimension_flag=.false.
			SplitDimensionFunc%boundarysize=dimen%boundarysize+1
			allocate(SplitDimensionFunc%boundary(SplitDimensionFunc%boundarysize))
			SplitDimensionFunc%boundary(:inde)=dimen%boundary(:inde)
			SplitDimensionFunc%boundary(inde+1)=dimen%boundary(inde) + midindde
			do i=inde+1,SplitDimensionFunc%boundarysize-1
				SplitDimensionFunc%boundary(i+1)=dimen%boundary(i)
			end do
		end if
		SplitDimensionFunc%Nameflag=dimen%Nameflag
		if(.not.dimen%Nameflag) return
		allocate(SplitDimensionFunc%DimName(Dimen%lenDimData))
		SplitDimensionFunc%DimName=dimen%DimName(1:Dimen%lenDimData)
		return
	end function

	!
		!Will not store boundary as it is sample dimension	
		!If dimen is already the simple type, do nothing	

	function SplitDimensionFuncAll(dimen)
		class(Dimension),allocatable::SplitDimensionFuncAll
		class(Dimension),intent(in) :: Dimen
		integer :: i
		allocate(SplitDimensionFuncAll,mold=dimen)
		if(Dimen%outlenDimdata().le.0)then
			call SplitDimensionFuncAll%empty()
			return
		end if
		if(dimen%sample_dimension_flag)then
			SplitDimensionFuncAll=dimen
			return
		end if
     	SplitDimensionFuncAll%sample_dimension_flag=.true.
     	SplitDimensionFuncAll%LenDimData=dimen%LenDimData
     	allocate(SplitDimensionFuncAll%DimData(dimen%LenDimData))
		SplitDimensionFuncAll%DimData=dimen%DimData(1:dimen%LenDimData)
		SplitDimensionFuncAll%boundarysize=0
		SplitDimensionFuncAll%Dimsize=0
		SplitDimensionFuncAll%Nameflag=dimen%Nameflag
		if(.not.dimen%Nameflag) return
		allocate(SplitDimensionFuncAll%DimName(Dimen%lenDimData))
		SplitDimensionFuncAll%DimName=dimen%DimName(1:Dimen%lenDimData)
		return
	end function

	!
		! decompose the de_index index of the Tensor into n(1),n(2)
		!		for example the de_index index is (1,2,3,4,..inde,inde+1,...rank)
		!		(1,2,3,4,..inde,inde+1,...rank)-->(1,2,3,4,..inde),(inde+1,...rank)		
		!		if inde larger than rank ,the function will return no change	

	subroutine splitDimension(Dimen,de_index,inde)
		class(Dimension),intent(inout) :: Dimen
		integer,optional,intent(in) :: de_index
		integer,optional,intent(in)::inde
		if(present(de_index))then
			if(present(inde))then
				Dimen=SplitDimensionFunc(dimen,de_index,inde)
			else
				Dimen=SplitDimensionFunc(dimen,de_index,1)
			end if
		else
			Dimen=SplitDimensionFuncAll(Dimen)
		end if
		return
	end subroutine	
	function splitDimension2(dimen,vector)result(splitDimension)
		class(Dimension),allocatable::splitDimension
		class(Dimension),intent(in) :: dimen
		integer,intent(in) ::vector(2)
		integer:: de_index,inde
		allocate(splitDimension,mold=dimen)
		de_index=vector(1)
		inde=vector(2)
		splitDimension=SplitDimensionFunc(dimen,de_index,inde)
		return
	end function
	function splitDimension3(dimen,de_index)result(splitDimension)
		class(Dimension),allocatable::splitDimension
		class(Dimension),intent(in) :: dimen
		integer,intent(in) :: de_index
		allocate(splitDimension,mold=dimen)
		splitDimension=SplitDimensionFunc(dimen,de_index,1)
		return
	end function
	function splitDimensionAll(Dimen)result(splitDimension)
		class(Dimension),allocatable::splitDimension
		class(Dimension),intent(in) :: Dimen
		allocate(splitDimension,mold=dimen)
		splitDimension=SplitDimensionFuncAll(Dimen)
		return
	end function	

	!***********************************************
	!********* permuation                  *********
	!***********************************************


	function Dimpermute(dimen,v)
		class(Dimension),allocatable::Dimpermute
		class(Dimension),intent(in) :: Dimen
		integer,intent(in):: v(Dimen%dimSize)
		character(len=len_of_Name),allocatable::dimenName(:)
		integer::i,datalen,subDlen,k
		integer,allocatable :: dimenVec(:)
		allocate(Dimpermute,mold=dimen)
		datalen=dimen%LenDimData
		if(datalen.eq.0) then
			write(*,*)"There is no data in Dimension"
			call error_stop()
		end if
		Dimpermute%LenDimData=datalen
		allocate(Dimpermute%DimData(datalen))
		if(dimen%NameFlag)then
			allocate(Dimpermute%DimName(datalen))
		end if
		Dimpermute%Nameflag=dimen%NameFlag
		Dimpermute%sample_dimension_flag=dimen%sample_dimension_flag
		if(dimen%sample_dimension_flag)then
			do i=1,datalen
				if(dimen%NameFlag)then
					Dimpermute%DimName(i)=dimen%DimName(v(i))
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
			if(dimen%NameFlag)then
				call getSubDimName(Dimen,v(i),dimenName)
				Dimpermute%DimName(k:subDlen+k-1)=dimenName
			end if
			k=k+subDlen
		end do
	end function

	! Dimpermute_forwards: output order is [ith,1,2,3...]

	subroutine  Dimpermute_forwards(outdim,dimen,ith)
		class(Dimension),intent(inout)::outdim
		class(Dimension),intent(in) :: Dimen
		integer,intent(in):: ith
		integer,pointer::order(:)
		integer::length,i
		call WorkingMemory%Check()
		length=dimen%getRank()
		!allocate(order(length))
		call WorkingMemory%get_memory(order,length)
		order(1)=ith
		do i=2,ith
			order(i)=i-1
		end do
		do i=ith+1,length
			order(i)=i
		end do
		outdim=Dimpermute(dimen,order)
		call WorkingMemory%free()
		return
	end subroutine

	! Dimpermute_backwards: output order is [1,2,3...,n,ith]

	subroutine Dimpermute_backwards(outdim,dimen,ith)
		class(Dimension),intent(inout)::outdim
		class(Dimension),intent(in) :: Dimen
		integer,intent(in):: ith
		integer,pointer::order(:)
		integer::length,i
		call WorkingMemory%Check()
		length=dimen%getRank()
		!allocate(order(length))
		call WorkingMemory%get_memory(order,length)
		order(length)=ith
		do i=1,ith-1
			order(i)=i
		end do
		do i=ith,length-1
			order(i)=i+1
		end do
		outdim=Dimpermute(dimen,order)
		call WorkingMemory%free()
		return
	end subroutine

	!Dimpermute_forwards_index: output oeder is [2,3,4,...,ith,1,ith+1,...]

	subroutine Dimpermute_forwards_index(outdim,dimen,ith)
		class(Dimension),intent(inout)::outdim
		class(Dimension),intent(in) :: Dimen
		integer,intent(in):: ith
		integer,pointer::order(:)
		integer::length,i
		call WorkingMemory%Check()
		length=dimen%getRank()
		call WorkingMemory%get_memory(order,length)
		order(ith)=1
		do i=1,ith-1
			order(i)=i+1
		end do
		do i=ith+1,length
			order(i)=i
		end do
		outdim=Dimpermute(dimen,order)
		call WorkingMemory%free()
		return
	end subroutine

	! Dimpermute_backwards_index: output order is [1,2,3,4,...,ith,n,ith+1,...,n-1]

	subroutine Dimpermute_backwards_index(outdim,dimen,ith)
		class(Dimension),intent(inout)::outdim
		class(Dimension),intent(in) :: Dimen
		integer,intent(in):: ith
		integer,pointer::order(:)
		integer::length,i
		call WorkingMemory%Check()
		length=dimen%getRank()
		call WorkingMemory%get_memory(order,length)
		order(ith)=length
		do i=1,ith-1
			order(i)=i
		end do
		do i=ith+1,length
			order(i)=i-1
		end do
		outdim=Dimpermute(dimen,order)
		call WorkingMemory%free()
		return
	end subroutine


	!***********************************************
	!*********  logical function           *********
	!***********************************************

	function  if_long_Name(w_) result(res)
		logical::Res
		character(len=*),intent(in)::w_
		character(len=len(trim(adjustl(w_))))::w
		integer::lenw
		w=trim(adjustl(w_))
		lenw=index(w,indexsymbol)
		if(lenw.eq.0)then!input asd
			res=.false.
			return
		end if
		res=.true.
		return
	end function

	logical function equal_of_dim(dim1,dim2)
		class(Dimension),intent(in) :: dim1
		class(Dimension),intent(in) :: dim2
		integer,allocatable :: dimenVec1(:)
		integer,allocatable :: dimenVec2(:)
		call copydimension(dimenVec1,dim1)
		call copydimension(dimenVec2,dim2)
		equal_of_dim=dimenVec1.equ.dimenVec2
		return
	end function

	logical function equal_of_array_old(a,b)result(equal_of_array)
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

	logical function if_original_dim(dimen)
		class(Dimension),intent(in)::dimen
		if_original_dim=.true.
		if(dimen%sample_dimension_flag) return
		if(dimen%LenDimData.eq.dimen%Dimsize) return	
		if_original_dim=.false.
		return
	end function

	logical function Same_Name_Order_forwards(dimen,charorder)result(Res)
		class(dimension),intent(in)::dimen
		character(len=*),intent(in)::charorder(:)
		integer::i,lenchar
		lenchar=size(charorder)
		Res=.true.
		do i=1,lenchar
			if(dimen%getName(i).nequ.charorder(i))then
				Res=.false.
				return
			end if
		end do
		return
	end function
	logical function Same_Name_Order_Backwards(dimen,charorder)result(Res)
		class(dimension),intent(in)::dimen
		character(len=*),intent(in)::charorder(:)
		integer::i,lenchar,rank,j
		lenchar=size(charorder)
		rank=dimen%getRank()
		Res=.true.
		j=rank
		do i=lenchar,1,-1
			if(dimen%getName(j).nequ.charorder(i))then
				Res=.false.
				return
			end if
			j=j-1
		end do
		return
	end function

	!***********************************************
	!******** conbine two dimesion function   ******
	!***********************************************

	!
		! The type of input are both dimension	
		!If dimen have a name but Dimen2 do not or Dimen2 have but not dimen
		!then the name for the one with no name will be ['0',0,i]

	function  AddFunc(Dimen,Dimen2)
		class(Dimension),allocatable::AddFunc
		class(Dimension),intent(in) :: Dimen
		class(Dimension),intent(in) :: Dimen2
		integer,allocatable::boundary1(:),boundary2(:)
		integer::boundarysize1,boundarysize2,Dimsize1,Dimsize2
		integer::i,j,l1,l2,templen
		allocate(AddFunc,mold=Dimen)
		l1=dimen%LenDimData
		l2=dimen2%LenDimData
		if(l1.eq.0)then
			call writemess(' There is no Data in the first dimension when dim1 + dim2')
			call error_stop
		end if
		if(l2.eq.0)then
			call writemess(' There is no Data in the second dimension when dim1 + dim2')
			call error_stop
		end if
		AddFunc%lenDimData=l1+l2
		allocate(AddFunc%Dimdata(l1+l2))
		AddFunc%Dimdata(:l1)=Dimen%DimData(1:Dimen%lenDimData)
		AddFunc%Dimdata(l1+1:)=Dimen2%DimData(1:Dimen2%lenDimData)
		if( (.not.Dimen%nameFlag) .and. (.not.Dimen2%nameFlag) ) then
			AddFunc%nameFlag=.false.
			if(Dimen%sample_dimension_flag.and.Dimen2%sample_dimension_flag)then
				AddFunc%sample_dimension_flag=.true.
				return
			end if
		end if
		if((.not.Dimen%nameFlag) .and. (Dimen2%nameFlag))then
			allocate(AddFunc%DimName(l1+l2))
			do i=1,l1
				AddFunc%DimName(i)='0'+indexsymbol+i
			end do
			AddFunc%DimName(l1+1:)=Dimen2%DimName(1:l2)
			AddFunc%nameFlag=.true.
		end if
		if((Dimen%nameFlag) .and. (Dimen2%nameFlag))then
			allocate(AddFunc%DimName(l1+l2))
			AddFunc%DimName(1:l1)=Dimen%DimName(1:l1)
			AddFunc%DimName(l1+1:)=Dimen2%DimName(1:l2)
			AddFunc%nameFlag=.true.
		end if
		if((Dimen%nameFlag) .and. (.not.Dimen2%nameFlag))then
			allocate(AddFunc%DimName(l1+l2))
			AddFunc%DimName(1:l1)=Dimen%DimName(1:l1)
			do i=l1+1,l1+l2
				AddFunc%DimName(i)='0'+indexsymbol+i
			end do
			AddFunc%nameFlag=.true.
		end if
		
		if(check_same_name_Flag)then
			if(AddFunc%nameFlag)call check_same_name_in_dimension(AddFunc)
		end if
		
		if(Dimen%sample_dimension_flag.and.Dimen2%sample_dimension_flag)then
			AddFunc%sample_dimension_flag=.true.
			return
		end if
		AddFunc%sample_dimension_flag=.false.
		call outDimension_boundry(dimen,boundary1,boundarysize1,Dimsize1)
		call outDimension_boundry(dimen2,boundary2,boundarysize2,Dimsize2)
		AddFunc%boundarysize=boundarysize1+boundarysize2-1
		allocate(AddFunc%boundary(AddFunc%boundarysize))
		do i=1,boundarysize1
			AddFunc%boundary(i)=boundary1(i)
		end do
		templen=boundary1(boundarysize1)
		do i=1,boundarysize2-1
			AddFunc%boundary(i+boundarysize1)=boundary2(i+1)+templen
		end do
		AddFunc%Dimsize=Dimsize1+Dimsize2
		return
	end function

	! The type of input are one dimension and the other vector	

	function  AddFunc2(Dimen,Dimenvec)
		class(Dimension),allocatable::AddFunc2
		class(Dimension),intent(in) :: Dimen
		integer,intent(in) :: Dimenvec(:)
		integer,allocatable::boundary(:)
		integer::boundarysize,Dimsize
		integer::i,j,l1,l2,tempLen
		allocate(AddFunc2,mold=Dimen)
		l1=dimen%LenDimData
		l2=size(Dimenvec)
		if(l1.eq.0)then
			call writemess(' There is no Data in the first dimension when dim1 + vec(:)')
			call error_stop
		end if
		if(l2.eq.0)then
			call writemess(' There is no Data in the array when dim1 + vec(:)')
			call error_stop
		end if
		AddFunc2%lenDimData=l1+l2
		allocate(AddFunc2%Dimdata(l1+l2))
		AddFunc2%Dimdata(:l1)=Dimen%DimData(1:dimen%lenDimData)
		AddFunc2%Dimdata(l1+1:)=Dimenvec
		AddFunc2%nameFlag=Dimen%nameFlag
		AddFunc2%sample_dimension_flag=Dimen%sample_dimension_flag
		if((.not.AddFunc2%nameFlag).and.(AddFunc2%sample_dimension_flag)) return
				
			
		if(Dimen%nameFlag)then
			allocate(AddFunc2%DimName(l1+l2))
			AddFunc2%DimName(1:l1)=Dimen%DimName(1:l1)
			do i=l1+1,l1+l2
				AddFunc2%DimName(i)='0'+indexsymbol+i
			end do
		end if
		
		if(check_same_name_Flag)then
			if(AddFunc2%nameFlag)call check_same_name_in_dimension(AddFunc2)
		end if
		if(AddFunc2%sample_dimension_flag) return
		
		call outDimension_boundry(dimen,boundary,boundarysize,Dimsize)
		
		AddFunc2%boundarysize=boundarysize+size(Dimenvec)
		allocate(AddFunc2%boundary(AddFunc2%boundarysize))
		do i=1,boundarysize
			AddFunc2%boundary(i)=boundary(i)
		end do
		tempLen=boundary(boundarysize)
		do i=1,size(Dimenvec)
			AddFunc2%boundary(i+boundarysize)=i+tempLen
		end do
		AddFunc2%Dimsize=Dimsize+size(Dimenvec)
		return
	end function
	
	function  AddFunc3(Dimenvec,Dimen)
		class(Dimension),allocatable::AddFunc3
		class(Dimension),intent(in) :: Dimen
		integer,intent(in) :: Dimenvec(:)
		class(Dimension),allocatable::tempdimension
		allocate(AddFunc3,mold=Dimen)
		allocate(tempdimension,mold=Dimen)
		tempdimension=Dimenvec
		AddFunc3=AddFunc(tempdimension,Dimen)
	end function

	!**********************************************************************
	!**********************************************************************
	!	the code below is for MPI
	!**********************************************************************

	
	subroutine send_Dimension(Dimen1,Dimen2,ID1,ID2,ierr,MPIcommon)
		class(Dimension),intent(in)::Dimen1
		class(Dimension),intent(inout)::Dimen2
		integer,intent(in)::ID1,ID2
		integer,optional::ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,lenDimData,istatus(MPI_STATUS_SIZE),i,mpi_comm
		logical::nameflag
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		lenDimData=0
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		
		if(present(MPIcommon))then
			if((ID1.ge.proNum).or.(ID2.ge.proNum))return
		end if
		
		
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
			call mpi_send(Dimen1%NameFlag,1,MPI_logical,ID2,tag,mpi_comm,ierr)
			nameflag=Dimen1%NameFlag
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%NameFlag,1,MPI_logical,ID1,tag,mpi_comm,istatus,ierr)
			nameflag=Dimen2%NameFlag
			if(Dimen2%NameFlag) then
				call allocateCheck(Dimen2%DimName,dimen2%lenDimData)
			end if
		end if
		!******************   dimensionName  ***********************************************************
		if(nameflag)then
			if(proID.eq.ID1) then
				call mpi_send(Dimen1%DimName(1:dimen1%lenDimData),len_of_Name*dimen1%lenDimData,MPI_CHARACTER,ID2,tag,mpi_comm,ierr)
			end if
			if(proID.eq.ID2) then
				call mpi_recv(Dimen2%DimName(1:dimen2%lenDimData),len_of_Name*dimen2%lenDimData,MPI_CHARACTER,ID1,tag,mpi_comm,istatus,ierr)
			end if
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
		class(Dimension),intent(inout)::Dimen1
		integer,intent(in)::ID
		integer,optional::ierr
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
		
		if(present(MPIcommon))then
			if(ID.ge.proNum)return
		end if
		
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
		if(Dimen1%NameFlag)then
			if(proID.ne.ID) then
				call allocateCheck(Dimen1%DimName,Dimen1%LenDimData)
			end if
			call MPI_BCAST(Dimen1%DimName(1:dimen1%lenDimData),len_of_Name*dimen1%lenDimData,MPI_CHARACTER,ID,mpi_comm,ierr)	
		end if
		!******************   sample_dimension_flag  *********************************
		call MPI_BCAST(Dimen1%sample_dimension_flag,1,MPI_logical,ID,mpi_comm,ierr)
		if(Dimen1%sample_dimension_flag) return
		!*******************   boundarysize   ***************************************		
		call MPI_BCAST(Dimen1%boundarysize,1,MPI_integer,ID,mpi_comm,ierr)	
		!*********************   boundary   *****************************************
		if(proID.ne.ID) then
			call allocateCheck(Dimen1%boundary,Dimen1%boundarysize)
		end if
		call MPI_BCAST(Dimen1%boundary(1:Dimen1%boundarysize),Dimen1%boundarysize,MPI_integer,ID,mpi_comm,ierr)		
		!********************   Dimsize   ******************************************
		call MPI_BCAST(Dimen1%Dimsize,1,MPI_integer,ID,mpi_comm,ierr)	
		return
	end subroutine

end module