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
!
!************************************************************
!************* START OF Tensor **********************
!************************************************************
!*****************     worning    *******************
!			1.when link1=link2,and want link1 unchange while operate on link2,use
!		copylink but not (=)
!			2.The code is fortran90 version, gfortran4.8.4.
!			3.to compile the code one should link the files to lapack and blas.
!			4.Send email to sj.dong@outlook.com to report any bugs.
!
!
!
!
!1.  allocateTensor
!2.  dallocateTensor
!3.  assignment to Tensor
!4.  assignment to array
!5 . Type Overloading
!6.  print Tensor 
!7.  get  dimension   data
!8.  TensorName
!9.  element
!10.  generate data
!11. modify element in Tensor
!12. + - * /
!13. overwrite int real dble cmplx dcmplx aimag char
!14. max or min element in Tensor
!15. operaction on dimension
!16. permutation
!17. contract
!18. useful function
!19. lapack function (SVD, linear equation inverse and so on)
!20. Tensorlink and Tensornode
!21. MPI function

module Tensor_type
	use TData_module
	use Dimension_typede
	use Contract_real8_Tools
	use mpi
	use Tools
	implicit none
	private
	logical::SVD_S_matrix_flag=.false.
	character(len=1)::array_character_divider='|'
	
	public::Tensor
	type,extends (Dimension) :: Tensor
		private
		type(TData) :: Data
		integer:: rank=0
	contains
		procedure::ipointer,spointer,dpointer,cpointer,zpointer,lpointer,apointer
		procedure::ipointer_,spointer_,dpointer_,cpointer_,zpointer_,lpointer_,apointer_
		procedure::ipointer2,spointer2,dpointer2,cpointer2,zpointer2,lpointer2,apointer2
		procedure::ipointer2_,spointer2_,dpointer2_,cpointer2_,zpointer2_,lpointer2_,apointer2_
		procedure::ipointer3,spointer3,dpointer3,cpointer3,zpointer3,lpointer3,apointer3
		procedure::ipointer3_,spointer3_,dpointer3_,cpointer3_,zpointer3_,lpointer3_,apointer3_
		procedure::ipointer4,spointer4,dpointer4,cpointer4,zpointer4,lpointer4,apointer4
		procedure::ipointer4_,spointer4_,dpointer4_,cpointer4_,zpointer4_,lpointer4_,apointer4_
		generic,public::pointer=>ipointer,spointer,dpointer,cpointer,zpointer,lpointer,apointer,&
								ipointer_,spointer_,dpointer_,cpointer_,zpointer_,lpointer_,apointer_,&
								ipointer2,spointer2,dpointer2,cpointer2,zpointer2,lpointer2,apointer2,&
								ipointer2_,spointer2_,dpointer2_,cpointer2_,zpointer2_,lpointer2_,apointer2_,&
								ipointer3,spointer3,dpointer3,cpointer3,zpointer3,lpointer3,apointer3,&
								ipointer3_,spointer3_,dpointer3_,cpointer3_,zpointer3_,lpointer3_,apointer3_,&
								ipointer4,spointer4,dpointer4,cpointer4,zpointer4,lpointer4,apointer4,&
								ipointer4_,spointer4_,dpointer4_,cpointer4_,zpointer4_,lpointer4_,apointer4_

		! function on T%Data
		procedure,public::getFlag
		procedure::setclassType,setclassType2
		generic,public::setType =>setclassType,setclassType2
		procedure,public::getclassType
		procedure,public::getType
		procedure,public::Dynamic =>setDynamic
		procedure,public::Static=>setStatic
		procedure,public::ifDynamic
		! over write diemsnion
		procedure:: emptyDimension!empty
		procedure:: cleanDimension!deallocate
		procedure::getDimSize!getRank
		procedure::getTotalData
!		procedure::TMprint2,TMprint3,TMprint4
		procedure::Dprint,TMprint2,TMprint4
 		generic,public::print =>TMprint2,TMprint4
 		procedure::Dprint2,Tprint2,Tprint4
 		generic,public::info =>Tprint2,Tprint4
 		!procedure::readdimension
 		!procedure,public::readData =>Tread_data
! 		procedure,public::read=>readdimension
! 		procedure,public::fuse=>fuseDimension
! 		procedure,public::split =>splitDimension
! 		procedure::resetDimension
! 		generic,public::resetDim =>resetDimension

! 		procedure::RNDimRoutine,RNDimRoutineint,RNDimRoutinechar
! 		generic,public::killLeg=>RNDimRoutine,RNDimRoutineint,RNDimRoutinechar
! 		procedure::RNDimFun,RNDimFunInt,RNDimFunChar
! 		generic,public::getkillLeg=>RNDimFun,RNDimFunInt,RNDimFunChar
! 		procedure,public::getSubDimlen
! 		procedure,public::if_original_dim
!       	procedure,public::ifName
!       	procedure::send_Dimension,BCAST_Dimension
!       	procedure,public::MPI_send_Dimension=>send_Dimension
!       	generic,public::Send=>send_Dimension
!       	procedure,public::MPI_BCAST_Dimension=>BCAST_Dimension
!       	generic,public::BCAST=>BCAST_Dimension
 		!********* assignment **************************
 		procedure::DimInitialization
 		procedure::DimInitialization2
 		procedure::assignmentTenarray_int2
 		procedure::assignmentTenarray_int3
 		procedure::assignmentTenarray_int4
 		procedure::assignmentTenarray_real4_1
 		procedure::assignmentTenarray_real4_2
 		procedure::assignmentTenarray_real4_3
 		procedure::assignmentTenarray_real4_4
 		procedure::assignmentTenarray_real8_1
 		procedure::assignmentTenarray_real8_2
 		procedure::assignmentTenarray_real8_3
 		procedure::assignmentTenarray_real8_4
 		procedure::assignmentTenarray_com4_1
 		procedure::assignmentTenarray_com4_2
 		procedure::assignmentTenarray_com4_3
 		procedure::assignmentTenarray_com4_4
 		procedure::assignmentTenarray_com8_1
 		procedure::assignmentTenarray_com8_2
 		procedure::assignmentTenarray_com8_3
 		procedure::assignmentTenarray_com8_4
 		procedure::assignmentTenarray_logi1
 		procedure::assignmentTenarray_logi2
 		procedure::assignmentTenarray_logi3
 		procedure::assignmentTenarray_logi4
 		procedure::assignmentTenarray_char1
 		procedure::assignmentTenarray_char2
 		procedure::assignmentTenarray_char3
 		procedure::assignmentTenarray_char4
 		procedure::assignmentTenNum_int
 		procedure::assignmentTenNum_real4
 		procedure::assignmentTenNum_real8
 		procedure::assignmentTenNum_com4
 		procedure::assignmentTenNum_com8
 		procedure::assignmentTenNum_logi
 		procedure::assignmentTenNum_char
 		procedure,pass(T)::assignmentNumTen_int
 		procedure,pass(T)::assignmentNumTen_real4
 		procedure,pass(T)::assignmentNumTen_real8
 		procedure,pass(T)::assignmentNumTen_com4
 		procedure,pass(T)::assignmentNumTen_com8
 		procedure,pass(T)::assignmentNumTen_logi
 		procedure,pass(T)::assignmentNumTen_char
 		!procedure,pass(Dimen)::copyDimToVec
 		generic,public :: assignment(=) =>  assignmentTenarray_int2,&
									 		assignmentTenarray_int3,&
									 		assignmentTenarray_int4,&
									 		assignmentTenarray_real4_1,&
									 		assignmentTenarray_real4_2,&
									 		assignmentTenarray_real4_3,&
									 		assignmentTenarray_real4_4,&
									 		assignmentTenarray_real8_1,&
									 		assignmentTenarray_real8_2,&
									 		assignmentTenarray_real8_3,&
									 		assignmentTenarray_real8_4,&
									 		assignmentTenarray_com4_1,&
									 		assignmentTenarray_com4_2,&
									 		assignmentTenarray_com4_3,&
									 		assignmentTenarray_com4_4,&
									 		assignmentTenarray_com8_1,&
									 		assignmentTenarray_com8_2,&
									 		assignmentTenarray_com8_3,&
									 		assignmentTenarray_com8_4,&
									 		assignmentTenarray_logi1,&
									 		assignmentTenarray_logi2,&
									 		assignmentTenarray_logi3,&
									 		assignmentTenarray_logi4,&
									 		assignmentTenarray_char1,&
									 		assignmentTenarray_char2,&
									 		assignmentTenarray_char3,&
									 		assignmentTenarray_char4,&
									 		assignmentTenNum_int,&
									 		assignmentTenNum_real4,&
									 		assignmentTenNum_real8,&
									 		assignmentTenNum_com4,&
									 		assignmentTenNum_com8,&
									 		assignmentTenNum_logi,&
									 		assignmentTenNum_char,&
									 		assignmentNumTen_int,&
									 		assignmentNumTen_real4,&
									 		assignmentNumTen_real8,&
									 		assignmentNumTen_com4,&
									 		assignmentNumTen_com8,&
									 		assignmentNumTen_logi,&
									 		assignmentNumTen_char

 ! 		!********* operators **************************
! 		procedure::fuseDimension_val,fuseDimension_vec
! 		generic,public :: operator(.fuse.)=>fuseDimension_val,fuseDimension_vec
! 		procedure::splitDimension2,splitDimension3,splitDimensionAll
! 		generic,public :: operator(.split.)=>splitDimension2,splitDimension3,splitDimensionAll
! 		procedure::AddFunc,AddFunc2
! 		procedure,pass(Dimen)::AddFunc3
! 		generic,public :: operator(+) =>AddFunc,AddFunc2,AddFunc3
! 		procedure::getSubDim2,getSubDim3,getSubDim2_name
! 		generic,public :: operator(.subdim.)=>getSubDim2,getSubDim3,getSubDim2_name
! 		generic,public :: operator(.dim.)=>DimName_i,Dim_i,AllDim
! 		procedure::equal_of_dim
! 		generic,public :: operator(.equ.)=>equal_of_dim
	end type

contains
	!***********************************************
	!********* assignment **************************
	!***********************************************

	subroutine DimInitialization2(Dimen,Dimen2)
		class(Tensor),intent(inout) ::Dimen
		class(Dimension),intent(in) ::Dimen2
		select type(Dimen2)
		class is (Tensor)
			if(.not.Dimen2%getFlag())then
				call Dimen%empty()
				return
			end if
			call Dimen%allocate(dimen2)
			call assignmentTData_routine(Dimen%Data,Dimen2%Data)
		class default
			call writemess('ERROR in assignment for tensor',-1)
			call error_stop
		end select
		return
	end subroutine
	subroutine DimInitialization(Dimen,DimData)
		class(Tensor),intent(inout) ::Dimen
		integer,intent(in) :: DimData(:)
		integer::length
		length=size(DimData)
		call allocatedTensor(Dimen,[length],1)
		call assignmentTData_int(Dimen%Data,DimData,length)
		return
	end subroutine
	subroutine assignmentTenarray_int2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		integer,intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,1)
		call assignmentTData_int(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_int3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		integer,intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,1)
		call assignmentTData_int(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_int4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		integer,intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,1)
		call assignmentTData_int(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real4_1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call allocatedTensor(T,(/length/),2)
		call assignmentTData_real4(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real4_2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,2)
		call assignmentTData_real4(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real4_3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,2)
		call assignmentTData_real4(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real4_4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,2)
		call assignmentTData_real4(T%Data,Tensor_data,length)
		return
	end subroutine
	
	subroutine assignmentTenarray_real8_1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call allocatedTensor(T,(/length/),3)
		call assignmentTData_real8(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real8_2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,3)
		call assignmentTData_real8(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real8_3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,3)
		call assignmentTData_real8(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real8_4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,3)
		call assignmentTData_real8(T%Data,Tensor_data,length)
		return
	end subroutine
	

	subroutine assignmentTenarray_com4_1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call allocatedTensor(T,(/length/),4)
		call assignmentTData_com4(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com4_2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,4)
		call assignmentTData_com4(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com4_3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,4)
		call assignmentTData_com4(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com4_4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,4)
		call assignmentTData_com4(T%Data,Tensor_data,length)
		return
	end subroutine
	
	subroutine assignmentTenarray_com8_1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call allocatedTensor(T,(/length/),5)
		call assignmentTData_com8(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com8_2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,5)
		call assignmentTData_com8(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com8_3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,5)
		call assignmentTData_com8(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com8_4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,5)
		call assignmentTData_com8(T%Data,Tensor_data,length)
		return
	end subroutine	

	

	subroutine assignmentTenarray_logi1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		logical,intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call allocatedTensor(T,(/length/),6)
		call assignmentTData_logi(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_logi2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		logical,intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,6)
		call assignmentTData_logi(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_logi3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		logical,intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,6)
		call assignmentTData_logi(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_logi4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		logical,intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,6)
		call assignmentTData_logi(T%Data,Tensor_data,length)
		return
	end subroutine
	
	subroutine assignmentTenarray_char1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call allocatedTensor(T,(/length/),7)
		call assignmentTData_char(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_char2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,7)
		call assignmentTData_char(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_char3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,7)
		call assignmentTData_char(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_char4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call allocatedTensor(T,dimen,7)
		call assignmentTData_char(T%Data,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenNum_int(T,num)
		class(Tensor),intent(inout)::T
		integer,intent(in)::num
		call DimInitialization(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_real4(T,num)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::num
		call assignmentTenarray_real4_1(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_real8(T,num)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::num
		call assignmentTenarray_real8_1(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_com4(T,num)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::num
		call assignmentTenarray_com4_1(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_com8(T,num)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::num
		call assignmentTenarray_com8_1(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_logi(T,num)
		class(Tensor),intent(inout)::T
		logical,intent(in)::num
		call assignmentTenarray_logi_1(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_char(T,num)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::num
		call assignmentTenarray_char_1(T,[num])
		return
	end subroutine
	subroutine assignmentNumTen_int(val,T)
		integer,intent(inout)::val
		class(Tensor),intent(in)::T
		if(getTotalData(T).eq.1) then
			call assignment_int_Tdata_value(val,T%Data)
		else 
			write(*,*)"ERROR in assignment for Tensor to integer"
			call error_stop()
		end if
		return
	end subroutine
	subroutine assignmentNumTen_real4(val,T)
		real(kind=4),intent(inout)::val
		class(Tensor),intent(in)::T
		if(getTotalData(T).eq.1) then
			call assignment_real4_Tdata_value(val,T%Data)
		else if(.not.T%getFlag())then
			call writemess("ERROR in assignment for Tensor to real",-1)
			call writemess("The Tensor is a empty Tensor",-1)
			call error_stop()
		else
			call writemess("ERROR in assignment for Tensor to real",-1)
			call T%print()
			call error_stop()
		end if
		return
	end 	subroutine 
	subroutine assignmentNumTen_real8(val,T)
		real(kind=8),intent(inout)::val
		class(Tensor),intent(in)::T
		if(getTotalData(T).eq.1) then
			call assignment_real8_Tdata_value(val,T%Data)
		else if(.not.T%getFlag())then
			call writemess("ERROR in assignment for Tensor to real",-1)
			call writemess("The Tensor is a empty Tensor",-1)
			call error_stop()
		else
			call writemess("ERROR in assignment for Tensor to real",-1)
			call T%print()
			call error_stop()
		end if
		return
	end subroutine
	subroutine assignmentNumTen_com4(val,T)
		complex(kind=4),intent(inout)::val
		class(Tensor),intent(in)::T
		if(getTotalData(T).eq.1) then
			call assignment_com4_Tdata_value(val,T%Data)
		else
			call writemess("ERROR in assignment for Tensor to complex",-1)
			call error_stop()
		end if
		return
	end 	subroutine
	subroutine  assignmentNumTen_com8(val,T)
		complex(kind=8),intent(inout)::val
		class(Tensor),intent(in)::T
		if(getTotalData(T).eq.1) then
			call assignment_com8_Tdata_value(val,T%Data)
		else
			write(*,*)"ERROR in assignment for Tensor to complex"
			call error_stop()
		end if
		return
	end subroutine
	subroutine  assignmentNumTen_logi(val,T)
		logical,intent(inout)::val
		class(Tensor),intent(in)::T
		if(getTotalData(T).eq.1) then
			call assignment_logi_Tdata_value(val,T%Data)
		else
			write(*,*)"ERROR in assignment for Tensor to complex"
			call error_stop()
		end if
		return
	end subroutine
	subroutine  assignmentNumTen_char(val,T)
		character(len=*),intent(inout)::val
		class(Tensor),intent(in)::T
		if(getTotalData(T).eq.1) then
			call assignment_char_Tdata_value(val,T%Data)
		else
			write(*,*)"ERROR in assignment for Tensor to complex"
			call error_stop()
		end if
		return
	end subroutine

	!*********************************************************
	!   get The info of the tensor
	!*********************************************************

	logical function getFlag(T)
		class(Tensor),intent(in) :: T
		getFlag=T%Data%getFlag()
	end function
	function getclassType(T)
		character(len=20)::getclassType
		class(Tensor),intent(in) :: T
		getclassType=T%Data%getclassType()
		return
	end function
	integer function getType(T)
		class(Tensor),intent(in) :: T
		getType=T%Data%getType()
	end function	
	subroutine setclassType(T,classType_)
		class(Tensor),intent(inout) :: T
		character(len=*),intent(in)::classType_
		call T%Data%setType(classType_)
		return
	end subroutine
	subroutine setclassType2(T,classType)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::classType
		call T%Data%setType(classType)
		return
	end subroutine
	subroutine setDynamic(T)
		class(Tensor),intent(inout) :: T
		call T%Data%Dynamic()
		return
	end subroutine
	subroutine setStatic(T)
		class(Tensor),intent(inout) :: T
		call T%Data%Static()
		return
	end subroutine
	logical function ifDynamic(T)
		class(Tensor),intent(in) :: T
		ifDynamic=T%Data%ifDynamic()
		return
	end function

	!*********************************************************
	!   over write the function of dimension
	!*********************************************************

	subroutine emptyDimension(Dimen)!make dimension to a a empty dimension.but do not deallocate
		class(Tensor),intent(inout)::	Dimen
		call Dimen%Dimension%empty()
		call Dimen%Data%empty()
		Dimen%rank=0
		return
	end subroutine
	subroutine cleanDimension(Dimen)
		class(Tensor),intent(inout)::	Dimen
		call Dimen%Dimension%deallocate()
		call Dimen%Data%deallocate()
		Dimen%rank=0
		return
	end subroutine
	integer	function getDimSize(Dimen)
		class(Tensor),intent(in) :: Dimen
		getDimSize=Dimen%rank
		return
	end function
	integer function getTotalData(Dimen)
		class(Tensor),intent(in) :: Dimen
		getTotalData=Dimen%Data%getTotalData()
		return
	end function	
	subroutine Dprint(Dimen,uni)
		class(Tensor),intent(in) ::Dimen
		integer,optional,intent(in)::uni
		if(present(uni))then
			call writemess(' This function is no longer use',-1)
			call error_stop
		else
			call TMprint1(Dimen)
		end if
	end subroutine

	subroutine Dprint2(Dimen,uni)
		class(Tensor),intent(in) ::Dimen
		integer,optional,intent(in)::uni
		if(present(uni))then
			call writemess(' This function is no longer use',-1)
			call error_stop
		else
			call Tprint1(Dimen)
		end if
	end subroutine

	!*********************************************************
	!  print the data of Tensor
	!*********************************************************

	subroutine TMprint1(T)
		class(Tensor),intent(in) :: T
		CHARACTER(len=20)::classTypeChar
		integer,allocatable::dimen(:)
		write(*,*)"=================="
		write(*,*)"------------------"
		classTypeChar=getclassType(T)
		if(ifDynamic(T))then
			write(*,*)'Dynamic,',classTypeChar
		else
			write(*,*)'Static,',classTypeChar
		end if
		write(*,*) "*** START ***"
		if(.not.T%getFlag())then
			write(*,*)"There is no data"
			write(*,*) "*** END ***"
			return
		end if
		
		select case(T%getRank())
			case(1)
				call Tprintdata(T%Data,0)
				write(*,*) "*** END ***"
			case(2)
				allocate(dimen(2))
				dimen=.dim.T
				call Tprint_as_matrix(T%Data,0,T%rank,dimen)
				write(*,*) "*** END ***"
			case(3)
				allocate(dimen(3))
				dimen=.dim.T
				call Tprint_as_matrix(T%Data,0,T%rank,dimen)
				write(*,*) "*** END ***"
			case(4)
				allocate(dimen(4))
				dimen=.dim.T
				call Tprint_as_matrix(T%Data,0,T%rank,dimen)
				write(*,*) "*** END ***"
			case default
				write(*,*) "rank of the Tensor is large than 4"
				classTypeChar=getclassType(T)
				if(ifDynamic(T))then
					write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
				else
					write(*,*)'static class Tensor,data type is,',classTypeChar
				end if
				write(*,*) "The data of the Tensor is"
				call Tprintdata(T%Data,0)
				write(*,*) "***end***"
				write(*,*) "The dimension of the Tensor is"
				call T%Dimension%print()
				write(*,*) "The rank,total data are"
				write(*,*) T%getRank(),T%getTotalData()
		end select
		return
	end subroutine

	subroutine TMprint2(T,words,printType)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		integer,allocatable::dimen(:)
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)trim(words)
		classTypeChar=T%getclassType()
		if(ifDynamic(T))then
			write(*,*)'Dynamic,',classTypeChar
		else
			write(*,*)'Static,',classTypeChar
		end if
		write(*,*) "*** START ***"
		if(.not.T%getFlag())then
			write(*,*)"There is no data"
			write(*,*) "*** END ***"
			return
		end if
		
		select case(T%getRank())
			case(1)
				call Tprintdata(T%Data,0,printType)
				write(*,*) "*** END ***"
			case(2)
				allocate(dimen(2))
				dimen=.dim.T
				call Tprint_as_matrix(T%Data,0,T%rank,dimen,printType)
				write(*,*) "*** END ***"
			case(3)
				allocate(dimen(3))
				dimen=.dim.T
				call Tprint_as_matrix(T%Data,0,T%rank,dimen,printType)
				write(*,*) "*** END ***"
			case(4)
				allocate(dimen(4))
				dimen=.dim.T
				call Tprint_as_matrix(T%Data,0,T%rank,dimen,printType)
				write(*,*) "*** END ***"
			case default
				write(*,*) "rank of the Tensor is large than 4"
				classTypeChar=getclassType(T)
				if(ifDynamic(T))then
					write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
				else
					write(*,*)'static class Tensor,data type is,',classTypeChar
				end if
				write(*,*) "The data of the Tensor is"
				call Tprintdata(T%Data,0,printType)
				write(*,*) "***end***"
				write(*,*) "The dimension of the Tensor is"
				call T%Dimension%print()
				write(*,*) "The rank,total data are"
				write(*,*) T%getRank(),T%getTotalData()
		end select
		return
	end subroutine
	subroutine TMprint4(T,words,realpart,printType)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		integer,intent(in)::realpart
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		integer,allocatable::dimen(:)
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)trim(words)
		classTypeChar=getclassType(T)
		if(ifDynamic(T))then
			write(*,*)'Dynamic,',classTypeChar
		else
			write(*,*)'Static,',classTypeChar
		end if
		write(*,*) "*** START ***"
		if(.not.T%getFlag())then
			write(*,*)"There is no data"
			write(*,*) "*** END ***"
			return
		end if
		
		select case(T%getRank())
			case(1)
				call Tprintdata(T%Data,realpart,printType)
				write(*,*) "*** END ***"
			case(2)
				allocate(dimen(2))
				dimen=.dim.T
				call Tprint_as_matrix(T%Data,realpart,T%rank,dimen,printType)
				write(*,*) "*** END ***"
			case(3)
				allocate(dimen(3))
				dimen=.dim.T
				call Tprint_as_matrix(T%Data,realpart,T%rank,dimen,printType)
				write(*,*) "*** END ***"
			case(4)
				allocate(dimen(4))
				dimen=.dim.T
				call Tprint_as_matrix(T%Data,realpart,T%rank,dimen,printType)
				write(*,*) "*** END ***"
			case default
				write(*,*) "rank of the Tensor is large than 4"
				classTypeChar=getclassType(T)
				if(ifDynamic(T))then
					write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
				else
					write(*,*)'static class Tensor,data type is,',classTypeChar
				end if
				write(*,*) "The data of the Tensor is"
				call Tprintdata(T,realpart,printType)
				write(*,*) "***end***"
				write(*,*) "The dimension of the Tensor is"
				call T%Dimension%print()
				write(*,*) "The rank,total data are"
				write(*,*) T%getRank(),T%getTotalData()
		end select
		return
	end subroutine

	subroutine Tprint1(T)
		class(Tensor),intent(in) :: T
		CHARACTER(len=20)::classTypeChar
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)""
		if(T%getflag()) then
			classTypeChar=getclassType(T)
			if(ifDynamic(T))then
				write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(*,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(*,*) "The rank of the Tensor is"
			write(*,*) T%getRank()
			write(*,*) "The number of  data of the Tensor is"
			write(*,*) getTotalData(T)
			write(*,*) "The data of the Tensor is"
			call Tprintdata(T%Data,0)
			call T%Dimension%print()
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
	subroutine Tprint2(T,words,printType)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)trim(words)
		if(getflag(T)) then
			classTypeChar=getclassType(T)
			if(ifDynamic(T))then
				write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(*,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(*,*) "The rank of the Tensor is"
			write(*,*) T%getRank()
			write(*,*) "The number of  data of the Tensor is"
			write(*,*) getTotalData(T)
			write(*,*) "The data of the Tensor is"
			call Tprintdata(T%Data,0,printType)
			call T%Dimension%print()
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
	subroutine Tprint4(T,words,realpart,printType)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		CHARACTER(len=*),optional,intent(in)::printType
		integer,intent(in)::realpart
		CHARACTER(len=20)::classTypeChar
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)trim(words)
		if(getflag(T)) then
			classTypeChar=getclassType(T)
			if(ifDynamic(T))then
				write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(*,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(*,*) "The rank of the Tensor is"
			write(*,*) T%getRank()
			write(*,*) "The number of  data of the Tensor is"
			write(*,*) getTotalData(T)
			write(*,*) "The data of the Tensor is"
			call Tprintdata(T%Data,realpart,printType)
			call T%Dimension%print()
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine

	!********************* print dimension *********************   
	
	subroutine TDprint1(T,words,uni)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		integer,intent(in)::uni
		CHARACTER(len=20)::classTypeChar
		write(uni,*)"=================="
		write(uni,*)"------------------"
		write(uni,*)trim(words)
		if(getflag(T)) then!if1
			write(uni,*) "*** START ***"
			classTypeChar=getclassType(T)
			if(ifDynamic(T))then
				write(uni,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(uni,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(uni,*) "The rank of the Tensor is"
			write(uni,*) T%rank
			write(uni,*) "The number of  data of the Tensor is"
			write(uni,*) getTotalData(T)
			call T%Dimension%print(uni)
			write(uni,*) "***end***"
			write(uni,*) ""
		else!if1
			write(uni,*) "There is no data"
		end if!if1
		return
	end subroutine
	subroutine TDprint2(T,words)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),optional,intent(in)::words
		CHARACTER(len=20)::classTypeChar
		call writemess("==================",-1)
		call writemess("------------------",-1)
		call writemess(words)
		if(getflag(T)) then!if1
			call writemess("*** START ***",-1)
			classTypeChar=getclassType(T)
			if(ifDynamic(T))then
				call writemess('Dynamic class Tensor,data type is '+(' '+classTypeChar),-1)
			else
				call writemess('static class Tensor,data type is '+(' '+classTypeChar),-1)
			end if
			call writemess("The rank of the Tensor is",-1)
			call writemess(''+T%rank)
			call writemess("The number of  data of the Tensor is",-1)
			call writemess(''+getTotalData(T),-1)
			call T%Dimension%print()
			call writemess( "***end***",-1)
			call writemess("",-1)
		else!if1
			call writemess("There is no data",-1)
		end if!if1
		return
	end subroutine
	subroutine TDprint3(T,uni)
		class(Tensor),intent(in) :: T
		integer,intent(in)::uni
		CHARACTER(len=20)::classTypeChar
		write(uni,*)"=================="
		write(uni,*)"------------------"
		if(getflag(T)) then!if1
			write(uni,*) "*** START ***"
			classTypeChar=getclassType(T)
			if(ifDynamic(T))then
				write(uni,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(uni,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(uni,*) "The rank of the Tensor is"
			write(uni,*) T%rank
			write(uni,*) "The number of  data of the Tensor is"
			write(uni,*) getTotalData(T)
			call T%Dimension%print(uni)
			write(uni,*) "***end***"
			write(uni,*) ""
		else!if1
			write(uni,*) "There is no data"
		end if!if1
		return
	end subroutine


	subroutine readdimension2(dimenaas,uni)
		class(Tensor),intent(inout)::dimenaas
		integer,intent(in)::uni
		CHARACTER(len=20)::classTypeChar
		CHARACTER(len=50)::notused
		integer::rank,TotalData,i
		integer,allocatable::idata(:)
		real*4,allocatable::sdata(:),sdata2(:)
		real*8,allocatable::ddata(:),ddata2(:)
		character(len=max_len_of_char_in_TData),allocatable::adata(:)
		logical,allocatable::ldata(:)
		read(uni,*)notused
		read(uni,*)notused
		if(notused.ne.'readable')then
			call writemess("error in reading",-1)
			call error_stop()
		end if
		read(uni,*)notused
		read(uni,*)notused
		read(uni,*)classTypeChar
		call dimenaas%empty()
		if(classTypeChar.equ.'END')return
		call dimenaas%setType(classTypeChar)
		read(uni,*) notused
		read(uni,*) rank
		read(uni,*) notused
		read(uni,*) TotalData
		select case(dimenaas%getType())
			case(1)
				allocate(idata(TotalData))
				read(uni,*) notused
				read(uni,*)(idata(i),i=1,TotalData)
				dimenaas=idata
			case(2)
				allocate(sdata(TotalData))
				read(uni,*) notused
				read(uni,*)(sdata(i),i=1,TotalData)
				dimenaas=sdata
			case(3)
				allocate(ddata(TotalData))
				read(uni,*) notused
				read(uni,*)(ddata(i),i=1,TotalData)
				dimenaas=ddata
			case(4)
				allocate(sdata(TotalData))
				allocate(sdata2(TotalData))
				read(uni,*) notused
				read(uni,*)(sdata(i),i=1,TotalData)
				read(uni,*) notused
				read(uni,*)(sdata2(i),i=1,TotalData)
				dimenaas=cmplx(sdata,sdata2,kind=4)
			case(5)
				allocate(ddata(TotalData))
				allocate(ddata2(TotalData))
				read(uni,*) notused
				read(uni,*)(ddata(i),i=1,TotalData)
				read(uni,*) notused
				read(uni,*)(ddata2(i),i=1,TotalData)
				dimenaas=dcmplx(ddata,ddata2)
			case(6)
				allocate(ldata(TotalData))
				read(uni,*) notused
				read(uni,*)(ldata(i),i=1,TotalData)
				dimenaas=ldata
			case(7)
				allocate(adata(TotalData))
				read(uni,*) notused
				read(uni,*)(adata(i),i=1,TotalData)
				dimenaas=adata
		end select
		call dimenaas%Dimension%read(uni)
		read(uni,*) notused
		return
	end subroutine
	
! 	subroutine Tread_data(T,uni)
! 		class(Tensor),intent(inout) :: T
! 		integer,intent(in)::uni
! 		CHARACTER(len=50)::notused
! 		integer::rank,TotalData,i,dim1,dim2,j
! 		integer,pointer::idata(:,:),iidata(:)
! 		real*4,pointer::sdata(:,:),ssdata(:)
! 		real*8,pointer::ddata(:,:),dddata(:)
! 		character(len=max_len_of_char_in_TData),pointer::adata(:,:),aadata(:)
! 		logical,pointer::ldata(:,:),lldata(:)
! 		if(T%getRank().gt.2)then
! 			call writemess('ERROR in reading data, only allow for rank<=2 Tensor',-1)
! 			call error_stop
! 		end if
! 		if(T%getRank().eq.1)then
! 			TotalData=T%getTotalData()
! 			select case(T%getType())
! 				case(1)
! 					call T%pointer(iidata)
! 					read(uni,*)(iidata(i),i=1,TotalData)
! 					nullify(iidata)
! 				case(2)
! 					call T%pointer(ssdata)
! 					read(uni,*)(ssdata(i),i=1,TotalData)
! 					nullify(ssdata)
! 				case(3)
! 					call T%pointer(dddata)
! 					read(uni,*)(dddata(i),i=1,TotalData)
! 					nullify(dddata)
! 				case(4)
! 					call writemess('ERROR in reading data, Tensor',-1)
! 					call writemess('Do not finished for complex data yet',-1)
! 					call writemess('You can real two real data and combine them into a complex one',-1)
! 					call writemess('for example:A and B are real*4 Tensor.',-1)
! 					call writemess(' call A%readData(unit1).',-1)
! 					call writemess(' call B%readData(unit2).',-1)
! 					call writemess(' C=cmplex(A,B).',-1)
! 					call error_stop
! 				case(5)
! 					call writemess('ERROR in reading data, Tensor',-1)
! 					call writemess('Do not finished for complex data yet',-1)
! 					call writemess('You can real two real data and combine them into a complex one',-1)
! 					call writemess('for example:A and B are real*8 Tensor.',-1)
! 					call writemess(' call A%readData(unit1).',-1)
! 					call writemess(' call B%readData(unit2).',-1)
! 					call writemess(' C=dcmplex(A,B).',-1)
! 					call error_stop
! 				case(6)
! 					call T%pointer(lldata)
! 					read(uni,*)(lldata(i),i=1,TotalData)
! 					nullify(lldata)
! 				case(7)
! 					call T%pointer(aadata)
! 					read(uni,*)(aadata(i),i=1,TotalData)
! 					nullify(aadata)
! 			end select
! 			return
! 		endif
! 		dim1=T%dim(1)
! 		dim2=T%dim(2)
! 		TotalData=T%getTotalData()
! 		select case(T%getType())
! 			case(1)
! 				call T%pointer(idata)
! 				do i=1,dim1
! 					read(uni,*)(idata(i,j),j=1,dim2)
! 				end do
! 				nullify(idata)
! 			case(2)
! 				call T%pointer(sdata)
! 				do i=1,dim1
! 					read(uni,*)(sdata(i,j),j=1,dim2)
! 				end do
! 				nullify(sdata)
! 			case(3)
! 				call T%pointer(ddata)
! 				do i=1,dim1
! 					read(uni,*)(ddata(i,j),j=1,dim2)
! 				end do
! 				nullify(ddata)
! 			case(4)
! 				call writemess('ERROR in reading data, Tensor',-1)
! 				call writemess('Do not finished for complex data yet',-1)
! 				call writemess('You can real two real data and combine them into a complex one',-1)
! 				call writemess('for example:A and B are real*4 Tensor.',-1)
! 				call writemess(' call A%readData(unit1).',-1)
! 				call writemess(' call B%readData(unit2).',-1)
! 				call writemess(' C=cmplex(A,B).',-1)
! 				call error_stop
! 			case(5)
! 				call writemess('ERROR in reading data, Tensor',-1)
! 				call writemess('Do not finished for complex data yet',-1)
! 				call writemess('You can real two real data and combine them into a complex one',-1)
! 				call writemess('for example:A and B are real*8 Tensor.',-1)
! 				call writemess(' call A%readData(unit1).',-1)
! 				call writemess(' call B%readData(unit2).',-1)
! 				call writemess(' C=dcmplex(A,B).',-1)
! 				call error_stop
! 			case(6)
! 				call T%pointer(ldata)
! 					do i=1,dim1
! 						read(uni,*)(ldata(i,j),j=1,dim2)
! 					end do
! 					nullify(ldata)
! 			case(7)
! 				call T%pointer(adata)
! 				do i=1,dim1
! 					read(uni,*)(adata(i,j),j=1,dim2)
! 				end do
! 				nullify(adata)
! 		end select
! 		return
! 	end subroutine

	!*********************************************************
	!  pointer function of Tensor
	!*********************************************************

	subroutine ipointer(T,p)
		class(Tensor),intent(in)::T
		integer,pointer,intent(inout)::p(:)
		call T%Data%pointer(p)
		return
	end subroutine
	subroutine ipointer_(T,p,i1i2)
		class(Tensor),intent(in)::T
		integer,pointer,intent(inout)::p(:)
		integer,intent(in)::i1i2(2)
		call T%Data%pointer(p,i1i2)
		return
	end subroutine
	subroutine spointer(T,p)
		class(Tensor),intent(in)::T
		real*4,pointer,intent(inout)::p(:)
		call T%Data%pointer(p)
		return
	end subroutine
	subroutine spointer_(T,p,i1i2)
		class(Tensor),intent(in)::T
		real*4,pointer,intent(inout)::p(:)
		integer,intent(in)::i1i2(2)
		call T%Data%pointer(p,i1i2)
		return
	end subroutine
	subroutine dpointer(T,p)
		class(Tensor),intent(in)::T
		real*8,pointer,intent(inout)::p(:)
		call T%Data%pointer(p)
		return
	end subroutine
	subroutine dpointer_(T,p,i1i2)
		class(Tensor),intent(in)::T
		real*8,pointer,intent(inout)::p(:)
		integer,intent(in)::i1i2(2)
		call T%Data%pointer(p,i1i2)
		return
	end subroutine
	subroutine cpointer(T,p)
		class(Tensor),intent(in)::T
		complex(kind=4),pointer,intent(inout)::p(:)
		call T%Data%pointer(p)
		return
	end subroutine
	subroutine cpointer_(T,p,i1i2)
		class(Tensor),intent(in)::T
		complex(kind=4),pointer,intent(inout)::p(:)
		integer,intent(in)::i1i2(2)
		call T%Data%pointer(p,i1i2)
		return
	end subroutine
	subroutine zpointer(T,p)
		class(Tensor),intent(in)::T
		complex(kind=8),pointer,intent(inout)::p(:)
		call T%Data%pointer(p)
		return
	end subroutine
	subroutine zpointer_(T,p,i1i2)
		class(Tensor),intent(in)::T
		complex(kind=8),pointer,intent(inout)::p(:)
		integer,intent(in)::i1i2(2)
		call T%Data%pointer(p,i1i2)
		return
	end subroutine
	subroutine lpointer(T,p)
		class(Tensor),intent(in)::T
		logical,pointer,intent(inout)::p(:)
		call T%Data%pointer(p)
		return
	end subroutine
	subroutine lpointer_(T,p,i1i2)
		class(Tensor),intent(in)::T
		logical,pointer,intent(inout)::p(:)
		integer,intent(in)::i1i2(2)
		call T%Data%pointer(p,i1i2)
		return
	end subroutine
	subroutine apointer(T,p)
		class(Tensor),intent(in)::T
		character(len=max_len_of_char_in_TData),pointer,intent(inout)::p(:)
		call T%Data%pointer(p)
		return
	end subroutine
	subroutine apointer_(T,p,i1i2)
		class(Tensor),intent(in)::T
		character(len=max_len_of_char_in_TData),pointer,intent(inout)::p(:)
		integer,intent(in)::i1i2(2)
		call T%Data%pointer(p,i1i2)
		return
	end subroutine

	subroutine ipointer2(T,p,i,j)
		class(Tensor),intent(in)::T
		integer,pointer,intent(inout)::p(:,:)
		integer,intent(in)::i(2),j(2)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.1)then
			call writemess('The type of Tensor is not  integer',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,d)
		return
	end subroutine
	subroutine ipointer2_(T,p)
		class(Tensor),intent(in)::T
		integer,pointer,intent(inout)::p(:,:)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.1)then
			call writemess('The type of Tensor is not  integer',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],d)
		return
	end subroutine
	subroutine spointer2(T,p,i,j)
		class(Tensor),intent(in)::T
		real*4,pointer,intent(inout)::p(:,:)
		integer,intent(in)::i(2),j(2)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.2)then
			call writemess('The Tensor is of real*4, one should use the real*4 pointer',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,d)
		return
	end subroutine
	subroutine spointer2_(T,p)
		class(Tensor),intent(in)::T
		real*4,pointer,intent(inout)::p(:,:)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.2)then
			call writemess('The type of Tensor is not  real*4',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],d)
		return
	end subroutine
	subroutine dpointer2(T,p,i,j)
		class(Tensor),intent(in)::T
		real*8,pointer,intent(inout)::p(:,:)
		integer,intent(in)::i(2),j(2)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.3)then
			call writemess('The type of Tensor is not  real*8',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,d)
		return
	end subroutine
	subroutine dpointer2_(T,p)
		class(Tensor),intent(in)::T
		real*8,pointer,intent(inout)::p(:,:)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.3)then
			call writemess('The type of Tensor is not  real*8',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],d)
		return
	end subroutine
	subroutine cpointer2(T,p,i,j)
		class(Tensor),intent(in)::T
		complex(kind=4),pointer,intent(inout)::p(:,:)
		integer,intent(in)::i(2),j(2)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.4)then
			call writemess('The type of Tensor is not  complex(kind=4)',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,d)
		return
	end subroutine
	subroutine cpointer2_(T,p)
		class(Tensor),intent(in)::T
		complex(kind=4),pointer,intent(inout)::p(:,:)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.4)then
			call writemess('The type of Tensor is not  complex(kind=4)',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],d)
		return
	end subroutine
	subroutine zpointer2(T,p,i,j)
		class(Tensor),intent(in)::T
		complex(kind=8),pointer,intent(inout)::p(:,:)
		integer,intent(in)::i(2),j(2)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.5)then
			call writemess('The type of Tensor is not  complex(kind=8)',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,d)
		return
	end subroutine
	subroutine zpointer2_(T,p)
		class(Tensor),intent(in)::T
		complex(kind=8),pointer,intent(inout)::p(:,:)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.5)then
			call writemess('The type of Tensor is not  complex(kind=8)',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],d)
		return
	end subroutine
	subroutine lpointer2(T,p,i,j)
		class(Tensor),intent(in)::T
		logical,pointer,intent(inout)::p(:,:)
		integer,intent(in)::i(2),j(2)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.6)then
			call writemess('The type of Tensor is not  logical',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,d)
		return
	end subroutine
	subroutine lpointer2_(T,p)
		class(Tensor),intent(in)::T
		logical,pointer,intent(inout)::p(:,:)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.6)then
			call writemess('The type of Tensor is not  logical',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],d)
		return
	end subroutine
	subroutine apointer2(T,p,i,j)
		class(Tensor),intent(in)::T
		character(len=max_len_of_char_in_TData),pointer,intent(inout)::p(:,:)
		integer,intent(in)::i(2),j(2)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.7)then
			call writemess('The type of Tensor is not  character(len=max_len_of_char_in_TData)',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,d)
		return
	end subroutine
	subroutine apointer2_(T,p)
		class(Tensor),intent(in)::T
		character(len=max_len_of_char_in_TData),pointer,intent(inout)::p(:,:)
		integer::d(2)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.7)then
			call writemess('The type of Tensor is not  character(len=max_len_of_char_in_TData)',-1)
			call error_stop
		end if
		if(T%getRank().ne.2)then
			call writemess('The rank of the Tensor is not 2, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],d)
		return
	end subroutine
	
	
	subroutine ipointer3(T,p,i,j,k)
		class(Tensor),intent(in)::T
		integer,pointer,intent(inout)::p(:,:,:)
		integer,intent(in)::i(2),j(2),k(2)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.1)then
			call writemess('The type of Tensor is not  integer',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,d)
		return
	end subroutine
	subroutine ipointer3_(T,p)
		class(Tensor),intent(in)::T
		integer,pointer,intent(inout)::p(:,:,:)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.1)then
			call writemess('The type of Tensor is not  integer',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],d)
		return
	end subroutine
	subroutine spointer3(T,p,i,j,k)
		class(Tensor),intent(in)::T
		real*4,pointer,intent(inout)::p(:,:,:)
		integer,intent(in)::i(2),j(2),k(2)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.2)then
			call writemess('The type of Tensor is not  real*4',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,d)
		return
	end subroutine
	subroutine spointer3_(T,p)
		class(Tensor),intent(in)::T
		real*4,pointer,intent(inout)::p(:,:,:)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.2)then
			call writemess('The type of Tensor is not  real*4',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],d)
		return
	end subroutine
	subroutine dpointer3(T,p,i,j,k)
		class(Tensor),intent(in)::T
		real*8,pointer,intent(inout)::p(:,:,:)
		integer,intent(in)::i(2),j(2),k(2)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.3)then
			call writemess('The type of Tensor is not  real*8',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,d)
		return
	end subroutine
	subroutine dpointer3_(T,p)
		class(Tensor),intent(in)::T
		real*8,pointer,intent(inout)::p(:,:,:)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.3)then
			call writemess('The type of Tensor is not  real*8',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],d)
		return
	end subroutine
	subroutine cpointer3(T,p,i,j,k)
		class(Tensor),intent(in)::T
		complex(kind=4),pointer,intent(inout)::p(:,:,:)
		integer,intent(in)::i(2),j(2),k(2)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.4)then
			call writemess('The type of Tensor is not  complex(kind=4)',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,d)
		return
	end subroutine
	subroutine cpointer3_(T,p)
		class(Tensor),intent(in)::T
		complex(kind=4),pointer,intent(inout)::p(:,:,:)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.4)then
			call writemess('The type of Tensor is not  complex(kind=4)',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],d)
		return
	end subroutine
	subroutine zpointer3(T,p,i,j,k)
		class(Tensor),intent(in)::T
		complex(kind=8),pointer,intent(inout)::p(:,:,:)
		integer,intent(in)::i(2),j(2),k(2)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.5)then
			call writemess('The type of Tensor is not  complex(kind=8)',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,d)
		return
	end subroutine
	subroutine zpointer3_(T,p)
		class(Tensor),intent(in)::T
		complex(kind=8),pointer,intent(inout)::p(:,:,:)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.5)then
			call writemess('The type of Tensor is not  complex(kind=8)',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],d)
		return
	end subroutine
	subroutine lpointer3(T,p,i,j,k)
		class(Tensor),intent(in)::T
		logical,pointer,intent(inout)::p(:,:,:)
		integer,intent(in)::i(2),j(2),k(2)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.6)then
			call writemess('The type of Tensor is not  logical',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,d)
		return
	end subroutine
	subroutine lpointer3_(T,p)
		class(Tensor),intent(in)::T
		logical,pointer,intent(inout)::p(:,:,:)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.6)then
			call writemess('The type of Tensor is not  logical',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],d)
		return
	end subroutine
	subroutine apointer3(T,p,i,j,k)
		class(Tensor),intent(in)::T
		character(len=max_len_of_char_in_TData),pointer,intent(inout)::p(:,:,:)
		integer,intent(in)::i(2),j(2),k(2)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.7)then
			call writemess('The type of Tensor is not  character(len=max_len_of_char_in_TData)',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,d)
		return
	end subroutine
	subroutine apointer3_(T,p)
		class(Tensor),intent(in)::T
		character(len=max_len_of_char_in_TData),pointer,intent(inout)::p(:,:,:)
		integer::d(3)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.7)then
			call writemess('The type of Tensor is not  character(len=max_len_of_char_in_TData)',-1)
			call error_stop
		end if
		if(T%getRank().ne.3)then
			call writemess('The rank of the Tensor is not 3, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],d)
		return
	end subroutine

	subroutine ipointer4(T,p,i,j,k,l)
		class(Tensor),intent(in)::T
		integer,pointer,intent(inout)::p(:,:,:,:)
		integer,intent(in)::i(2),j(2),k(2),l(2)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.1)then
			call writemess('The type of Tensor is not  integer',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,l,d)
		return
	end subroutine
	subroutine ipointer4_(T,p)
		class(Tensor),intent(in)::T
		integer,pointer,intent(inout)::p(:,:,:,:)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.1)then
			call writemess('The type of Tensor is not  integer',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],[1,d(4)],d)
		return
	end subroutine
	subroutine spointer4(T,p,i,j,k,l)
		class(Tensor),intent(in)::T
		real*4,pointer,intent(inout)::p(:,:,:,:)
		integer,intent(in)::i(2),j(2),k(2),l(2)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.2)then
			call writemess('The type of Tensor is not  real*4',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,l,d)
		return
	end subroutine
	subroutine spointer4_(T,p)
		class(Tensor),intent(in)::T
		real*4,pointer,intent(inout)::p(:,:,:,:)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.2)then
			call writemess('The type of Tensor is not  real*4',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],[1,d(4)],d)
		return
	end subroutine
	subroutine dpointer4(T,p,i,j,k,l)
		class(Tensor),intent(in)::T
		real*8,pointer,intent(inout)::p(:,:,:,:)
		integer,intent(in)::i(2),j(2),k(2),l(2)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.3)then
			call writemess('The type of Tensor is not  real*8',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,l,d)
		return
	end subroutine
	subroutine dpointer4_(T,p)
		class(Tensor),intent(in)::T
		real*8,pointer,intent(inout)::p(:,:,:,:)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.3)then
			call writemess('The type of Tensor is not  real*8',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],[1,d(4)],d)
		return
	end subroutine
	subroutine cpointer4(T,p,i,j,k,l)
		class(Tensor),intent(in)::T
		complex*8,pointer,intent(inout)::p(:,:,:,:)
		integer,intent(in)::i(2),j(2),k(2),l(2)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.4)then
			call writemess('The type of Tensor is not  complex(kind=4)',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,l,d)
		return
	end subroutine
	subroutine cpointer4_(T,p)
		class(Tensor),intent(in)::T
		complex*8,pointer,intent(inout)::p(:,:,:,:)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.4)then
			call writemess('The type of Tensor is not  complex(kind=4)',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],[1,d(4)],d)
		return
	end subroutine
	subroutine zpointer4(T,p,i,j,k,l)
		class(Tensor),intent(in)::T
		complex*16,pointer,intent(inout)::p(:,:,:,:)
		integer,intent(in)::i(2),j(2),k(2),l(2)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.5)then
			call writemess('The type of Tensor is not  complex(kind=8)',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,l,d)
		return
	end subroutine
	subroutine zpointer4_(T,p)
		class(Tensor),intent(in)::T
		complex*16,pointer,intent(inout)::p(:,:,:,:)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.5)then
			call writemess('The type of Tensor is not  complex(kind=8)',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],[1,d(4)],d)
		return
	end subroutine
	subroutine lpointer4(T,p,i,j,k,l)
		class(Tensor),intent(in)::T
		logical,pointer,intent(inout)::p(:,:,:,:)
		integer,intent(in)::i(2),j(2),k(2),l(2)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.6)then
			call writemess('The type of Tensor is not  logical',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,l,d)
		return
	end subroutine
	subroutine lpointer4_(T,p)
		class(Tensor),intent(in)::T
		logical,pointer,intent(inout)::p(:,:,:,:)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.6)then
			call writemess('The type of Tensor is not  logical',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],[1,d(4)],d)
		return
	end subroutine
	subroutine apointer4(T,p,i,j,k,l)
		class(Tensor),intent(in)::T
		character(len=characterLen),pointer,intent(inout)::p(:,:,:,:)
		integer,intent(in)::i(2),j(2),k(2),l(2)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.7)then
			call writemess('The type of Tensor is not  character(len=characterLen)',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,i,j,k,l,d)
		return
	end subroutine
	subroutine apointer4_(T,p)
		class(Tensor),intent(in)::T
		character(len=characterLen),pointer,intent(inout)::p(:,:,:,:)
		integer::d(4)
		if(T%gettotalData().eq.0)then
			call writemess('There is no data in the Tensor',-1)
			call error_stop
		end if
		if(T%getType().ne.7)then
			call writemess('The type of Tensor is not  character(len=characterLen)',-1)
			call error_stop
		end if
		if(T%getRank().ne.4)then
			call writemess('The rank of the Tensor is not 4, input error dimension pointer',-1)
			call error_stop
		end if
		d=T%dim()
		call T%Data%pointerDim(p,[1,d(1)],[1,d(2)],[1,d(3)],[1,d(4)],d)
		return
	end subroutine
end module