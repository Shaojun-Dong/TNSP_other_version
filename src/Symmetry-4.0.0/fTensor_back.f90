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

module fermiTensor
	use SymTensor_type
	use SymDimension_typede
	use Tensor_type
	use Dimension_typede
	use Tools
	use mpi
	use U1Tool
	use ParityTool
	use QuantumNumber_Type
	implicit none
	private

#if NoPointer

#if U1Flag
	interface ifParity
		module procedure U1ifParity4
	end interface

	interface QaunNumParity
		module procedure U1Parity1
	end interface

	interface fuseOrder
		module procedure U1fuseOrder
	end interface
	interface NonZeroElement
		module procedure U1NonZeroElement
	end interface

	interface Rule
		module procedure U1rule4
	end interface

	interface checkSymmetryRule
		module procedure checkU1Rule
	end interface

	interface reverseSymmetryRule
		module procedure reverseU1Rule
	end interface

#else

	interface ifParity
		module procedure ifParityIndex4
	end interface

	interface QaunNumParity
		module procedure pParity1
	end interface

	interface fuseOrder
		module procedure ParityfuseOrder
	end interface
	interface NonZeroElement
		module procedure ParityNonZeroElement
	end interface

	interface Rule
		module procedure PRule4
	end interface

	interface checkSymmetryRule
		module procedure checkParityRule
	end interface

	interface reverseSymmetryRule
		module procedure reverseParityRule
	end interface

#endif

#else
	procedure( ),public,pointer::checkSymmetryRule=>defaultcheckSymmetryRule
	procedure( ),public,pointer::reverseSymmetryRule=>defaultreverseSymmetryRule
	procedure(type(Tensor)),public,pointer::NonZeroElement=>defaultNonZeroElement
	procedure(logical),public,pointer::Rule=>defaultRule
	procedure(type(QuanNum)),public,pointer::fuseOrder=>defaultfuseOrder
	procedure(logical),public,pointer::ifParity=>defaultifParity
	procedure(integer),public,pointer::QaunNumParity=>defaultQaunNumParity
#endif

	logical::mointer_order_flag=.true.
	integer::fermi_factor=-1
	!if mointer_order_flag=.false. Do not moniter the rule=0
	!
	
	public::fTensor
	type,extends (SymTensor) :: fTensor
		
	contains
		generic,public::allocatefBlock=>allocatefBlock1,allocatefBlock2,allocatefBlock3,allocatefBlockAll
		generic,public::randomfTensor=>set_random_fTensor_val,set_random_fTensor_vec,set_random_fTensor_all	
		generic,public::setfValue=>set_Symelement_vec,set_Symelement_val!there is symmetry check in the input
		generic,public::SVDfTensor=>SVD,SVDName,SVDNum,SVDNameNum,SVDNameNumMaxValue
		generic,public::SVDfRoutine=>SVDRoutine,SVDRoutineNum,SVDNameNumRoutine,SVDNameNumMaxValueRoutine
		procedure::SVDRoutineu1f,SVDRoutineNumu1f,SVDNameNumMaxValueRoutine
		procedure::SVDNameNum,SVDNum,SVDNameNumRoutine,SVDNameNumMaxValue
		procedure,public::SymmetryCheck
		procedure,public::ifSymmetry
		generic,public::FermiSplit=>FermiSplitfTensor2,FermiSplitfTensor
		generic,public::FermiFuse=>FermiFusefTensor,FermiFusefTensor2,FermiFusefTensor3,FermiFusefTensor4&
									,FermiFusefTensor_,FermiFusefTensor2_,FermiFusefTensor4_,&
									FermiFusefTensorRule,FermiFusefTensorRule_
		generic,public::QRfTensor=>QRfTensor_noName,QRfTensor_name,QRfTensor_name2!It use QuanFuse not fermiFuse
		generic,public::QRfRoutine=>QRdecompositionfTensor,QRfTensorRoutine_name2
		generic,public::LQfTensor=>LQfTensor_noName,LQfTensor_name,LQfTensor_name2!It use QuanFuse not fermiFuse
		procedure,public::LQfRoutine=>LQdecompositionfTensor
		procedure,public::reorder=>DimOrder
		! reorder the Tensor, order the leg as the order in allTensorName(:)
		procedure,public::ReverseFermiArrow=>Reverse_Fermi_Rule_specify3
		procedure::contract_name_routine,contract_name_routine1,contract_name_routine2
		procedure::contract_name_routine4,contract_name_routine5,contract_name_routine6
		

		procedure::QRdecompositionfTensor,QRfTensorRoutine_name2
		procedure::permutefo_routine
		procedure::permutefo_name_routine
		procedure::permutefo_vec_routine
		procedure::permutefo_vec_name_routine
		procedure::permutefo_Tensor_routine
		procedure::permutation_routine
		procedure::permutation_name_routine
		procedure::permutation_Tensor_routine
		procedure::permuteback_routine
		procedure::permuteback_name_routine
		procedure::permuteback_vec_routine
		procedure::permuteback_vec_name_routine
		procedure::permuteback_Tensor_routine
		procedure::allocatefBlock1
		procedure::allocatefBlock2
		procedure::allocatefBlock3
		procedure::allocatefBlockAll
		procedure::set_random_fTensor_val
		procedure::set_random_fTensor_vec
		procedure::set_random_fTensor_all	
		procedure::set_Symelement_vec
		procedure::set_Symelement_val!As the dimension of the block is fix by the degeneracy of the symDimension
		                                           !This subroutine will ignore the input Tensor dimension
		procedure::SVD
		procedure::SVDName
		procedure::FermiSplitfTensor
		procedure::FermiSplitfTensor2
		procedure::FermiFusefTensor
		procedure::FermiFusefTensorRule
		procedure::FermiFusefTensor2
		procedure::FermiFusefTensor3
		procedure::FermiFusefTensor4
		procedure::FermiFusefTensor_
		procedure::FermiFusefTensorRule_
		procedure::FermiFusefTensor2_
		procedure::FermiFusefTensor4_
		procedure::QRfTensor_noName
		procedure::QRfTensor_name
		procedure::QRfTensor_name2
		procedure::LQfTensor_noName
		procedure::LQfTensor_name2
		procedure::LQfTensor_name
	end type fTensor
	
	public::Reverse_Fermi_Rule
	interface Reverse_Fermi_Rule
		module procedure Reverse_Fermi_Rule_specify
		module procedure Reverse_Fermi_Rule_specify2
		module procedure Reverse_Fermi_Rule_specify3
		module procedure Reverse_Fermi_Rule1
	end interface

	public::ReverseFermiArrow
	interface ReverseFermiArrow
		module procedure Reverse_Fermi_Rule_specify
		module procedure Reverse_Fermi_Rule_specify2
		module procedure Reverse_Fermi_Rule_specify3
		module procedure Reverse_Fermi_Rule1
	end interface
	
	public::un_set_mointer_order_flag
	public::set_mointer_order_flag
	
	
	
	public::expm
	interface expm
		module procedure  expmfTensor
	end interface
	
	public::operator(.pf.)
	interface operator(.pf.)
		module procedure FermiPermuteForWard
		module procedure FermiPermuteForWard_vec
		module procedure FermiPermuteForWard_name
		module procedure FermiPermuteForWard_name_vec
		module procedure FermiPermuteForWard_Tensor
	end interface
	public::operator(.pb.)
	interface operator(.pb.)
		module procedure FermiPermuteBackWard
		module procedure FermiPermuteBackWard_vec
		module procedure FermiPermuteBackWard_name
		module procedure FermiPermuteBackWard_name_vec
		module procedure FermiPermuteBackWard_Tensor
	end interface
	public::operator(.p.)
	interface operator(.p.)
		module procedure FermiPermutation
		module procedure FermiPermutation_name
		module procedure FermiPermutation_Tensor
	end interface
	
	
	public::contract
	interface contract
		module procedure contract_Same_name
		module procedure contract_noName
		module procedure contract_noName2
		module procedure contract_Name!In put dimension as character
		module procedure contract_Name2
		module procedure contract_name_int
		module procedure contract_name_int2
		module procedure contract_int_name
		module procedure contract_int_name2
	end interface
	public::operator(*)
	interface operator(*)
		module procedure ProductTensor
		module procedure multiply_number_real8
		module procedure multiply_number_com4
		module procedure multiply_number_real4
		
	end interface
	public::operator(/)
	interface operator(/)
		module procedure divide_real8
		module procedure divide_real4
		module procedure divide_Tensor
	end interface
	
	public::operator(+)
	interface operator(+)
		module procedure add_fTensor
	end interface
	
	public::operator(-)
	interface operator(-)
		module procedure minu_fTensor
	end interface
	
	public::operator(.H.)
	interface operator(.H.)
		module procedure Htranspose! Htranspose, and the symmetry rule will go inverse
	end interface
	public::operator(.Hn.)
	interface operator(.Hn.)
		module procedure Htranspose2! Htranspose, and the symmetry rule will go inverse,and Rename all TensorName
	end interface
	
	public::operator(.subdim.)!overwrite the function in type(Dimension)
	interface operator(.subdim.)
		module procedure getSubDim2
		module procedure getSubDim3
		module procedure getSubDim4
		module procedure getSubDim2_name
	end interface
	
	public::operator(.con.)
	interface operator(.con.)
		module procedure conjugate! conjugate
	end interface
	
	
	
	public::operator(.kron.)!direct Product,support any rank tensor and keep their TensorName,see more in help/operator
	interface operator(.kron.)
		module procedure directProductTensor
	end interface
	
	public::operator(.x.)
	interface operator(.x.)! dot product conjugating the first vector,The Tensor will be regard as a vector
		module procedure TdotTensor
	end interface
	
	public::assignment(=)
	interface assignment(=)
		module procedure assignmentTen
		module procedure assignmentreal4
		module procedure assignmentreal8
		module procedure assignmentcom4
		module procedure assignmentcom8
		module procedure fTensorToSymTensor
		module procedure fTensorToTensor
		module procedure SymTensorTofTensor!Do not check the Symmetry rule
		module procedure TensorTofTensor!Do not check the Symmetry rule
	end interface
	
	public::MPI_BCAST_fTensor,MPI_send_fTensor
	public::MPI_SUM_fTensor!Do not check the input fTensor
	public::set_fermi_symmetry
	interface MPI_SUM_fTensor
		module procedure MPI_SUM_fTensor1
	end interface

	logical,private::contract_firt_print=.true.

#if old_gcc
	logical,private::old_verion_gcc_flag=.true.
#else
	logical,private::old_verion_gcc_flag=.false.
#endif

contains

	subroutine print_contract_message()
		if(contract_firt_print)then
			call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
			call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
			call writemess('@@@                 NOTE                          @@@  ')
			call writemess('@@@ you have set old_verion_gcc_flag=1 in make.in @@@  ')
			call writemess('@@@ call T%contract(A,nameA,B,nameB) may cost     @@@  ')
			call writemess('@@@ more time !                                   @@@  ')
			call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
			call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
			contract_firt_print=.false.
		end if
	end subroutine


	subroutine defaultcheckSymmetryRule(Rule1,Rule2,lenname1,lenname2)
		integer,intent(in)::Rule1,Rule2
		character(len=*),intent(in)::lenname1,lenname2
	end subroutine

	subroutine defaultreverseSymmetryRule(T)
		Type(SymTensor),intent(inout)::T
		call error_mess()
	end subroutine

	type(Tensor) function defaultNonZeroElement(dimen)
		type(SymDimension),intent(in)::dimen
		call error_mess()
	end function
	logical function defaultRule(dimen,indices,length)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::length
		integer,intent(in)::indices(length)
		call error_mess()
	end function

	type(QuanNum) function defaultfuseOrder(order,Q1,Q2,newRule_)
		type(QuanNum),intent(in)::Q1,Q2
		type(Tensor),intent(inout)::order
		integer,intent(in)::newRule_
		call error_mess()
	end function
	logical function defaultifParity(dimen,vec,ith,jth,rank)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::rank
		integer,intent(in)::vec(rank)
		integer,intent(in)::ith,jth
		call error_mess()
	end function

	integer function defaultQaunNumParity(dimen,ith,jth)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::ith,jth
		call error_mess()
	end function

	subroutine error_mess()
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call writemess('@@@                 ERROR                   @@@  ')
		call writemess('@@@ You have not set the symmetry group yet @@@  ')
		call writemess('@@@                                         @@@  ')
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call error_stop
	end subroutine
#if NoPointer
	subroutine set_fermi_symmetry(group)
		character(len=*),intent(in)::group
		call writemess(' Can not set the symmetry group')
		return
	end subroutine
#else
	subroutine set_fermi_symmetry(group)
		character(len=*),intent(in)::group
		call writemess(' Set the symmetry group as '+(' '+group))
		if(group.equ.'U(1)')then
			checkSymmetryRule=>checkU1Rule
			reverseSymmetryRule=>reverseU1Rule
			NonZeroElement=>U1NonZeroElement
			Rule=>U1rule4
			fuseOrder=>U1fuseOrder
			ifParity=>U1ifParity4
			QaunNumParity=>U1Parity1
		else if(group.equ.'Parity')then
			checkSymmetryRule=>checkParityRule
			reverseSymmetryRule=>reverseParityRule
			NonZeroElement=>ParityNonZeroElement
			Rule=>PRule4
			fuseOrder=>ParityfuseOrder
			ifParity=>ifParityIndex4
			QaunNumParity=>pParity1
		else
			call writemess(' Only support the symmetry group of U(1) or Parity')
			call error_stop
		end if
		return
	end subroutine
#endif

	subroutine set_fermi_factor(factor)
		integer,intent(in)::factor
		fermi_factor=factor
		call writemess('Set fermi_factor as '+(' '+factor))
		return
	end subroutine
	logical function ifSymmetry(T)
		class(fTensor),intent(in)::T
		integer::i,rank
		integer,allocatable::Tdim(:),Maxdim(:)
		rank=T%GetRank()
		allocate(TDim(rank))
		allocate(Maxdim(rank))
		Maxdim=T%dim()
		ifSymmetry=.true.
		do i=1,T%getTotalData()
			if(T%getFlag(i))then
				call IndesToaddress(Maxdim,Tdim,i)	
				ifSymmetry=rule(T%SymDimension,TDim,rank)
				if(.not.ifSymmetry)return
			end if
		end do
	end function

	subroutine SymmetryCheck(T)
		class(fTensor),intent(in)::T
		integer::i,rank
		logical::checkFlag
		integer,allocatable::Tdim(:),Maxdim(:)
		rank=T%GetRank()
		allocate(TDim(rank))
		allocate(Maxdim(rank))
		Maxdim=T%dim()
		checkFlag=.true.
		do i=1,T%getTotalData()
			if(T%getFlag(i))then
				call IndesToaddress(Maxdim,Tdim,i)	
				checkFlag=rule(T%SymDimension,TDim,rank)
				if(.not.checkFlag)then
					call writemess(' The Tensor is not a Symmetry Tensor')
					call error_stop
				end if
			end if
		end do
	end subroutine
	
	subroutine assignmentTen(T,T2)
		type(fTensor),intent(inout) ::T
		type(fTensor),intent(in) :: T2
		T%SymTensor=T2%SymTensor
		return
	end subroutine
	subroutine assignmentreal4(r,T2)
		real*4,intent(inout)::r
		type(fTensor),intent(in) :: T2
		r=T2%SymTensor
		return
	end subroutine
	subroutine assignmentreal8(r,T2)
		real*8,intent(inout)::r
		type(fTensor),intent(in) :: T2
		r=T2%SymTensor
		return
	end subroutine
	subroutine assignmentcom4(r,T2)
		complex(kind=4),intent(inout)::r
		type(fTensor),intent(in) :: T2
		r=T2%SymTensor
		return
	end subroutine
	subroutine assignmentcom8(r,T2)
		complex(kind=8),intent(inout)::r
		type(fTensor),intent(in) :: T2
		r=T2%SymTensor
		return
	end subroutine
	subroutine fTensorToSymTensor(T,T2)
		type(SymTensor),intent(inout) ::T
		type(fTensor),intent(in) :: T2
		T=T2%SymTensor
		return
	end subroutine
	subroutine fTensorToTensor(T,T2)
		type(Tensor),intent(inout) ::T
		type(fTensor),intent(in) :: T2
		T=T2%SymTensor
		return
	end subroutine
	subroutine SymTensorTofTensor(T,T2)
		type(fTensor),intent(inout) ::T
		type(SymTensor),intent(in) :: T2
		T%SymTensor=T2
		return
	end subroutine
	subroutine TensorTofTensor(T,T2)
		type(fTensor),intent(inout) ::T
		type(Tensor),intent(in) :: T2
		if(.not.T%getFlag())then
			call writemess('Allocate fTensor before setting value')
			call error_stop()
		end if
		T%SymTensor=T2
		return
	end subroutine
	
	subroutine allocatefBlockAll(T)
		class(fTensor),intent(inout)::T
		type(Tensor)::elementindex
		elementindex=NonZeroElement(T%symDimension)
		if(.not.elementindex%getFlag())then
			call T%diminfo()
			call writemess('There is no non-zero block')
			call error_stop
		end if
		call T%SymTensor%allocateBlock(elementindex)
		return
	end subroutine
	subroutine allocatefBlock1(T,vec)
		class(fTensor),intent(inout)::T
		integer,intent(in)::vec(:)
		logical::goon
		integer::i
		character*500,allocatable::w
		real*4,allocatable::qnum(:)
		goon=rule(T%SymDimension,vec,T%GetRank())
		if(goon)then
			call T%SymTensor%allocateBlock(vec)
			return
		end if
		call writemess('Can not allocate the element of a fTensor')
		allocate(w)
		w='The indices are: ('
		do i=1,size(vec)-1
			w=w+vec(i)+','
		end do
		w=w+vec(i)+')'
		call writemess(w)
		allocate(qnum(size(vec)))
		qnum=T%QNDim(vec)
		w='The quantum number are: ('
		do i=1,size(vec)-1
			w=w+qnum(i)+','
		end do
		w=w+qnum(i)+')'
		call writemess(w)
		call error_stop()
	end 	subroutine
	
	subroutine allocatefBlock2(T,element)
		class(fTensor),intent(inout)::T
		type(Tensor),intent(in)::element
		integer,allocatable::Adim(:)
		if(element.gt.T%getTotalData())Then
			write(*,*)"Index is larger than the len of SymTensor"
			call error_stop()
		end if
		allocate(Adim(T%getRank()))
		call IndesToaddress(T%dim(),Adim,element%ii(1))
		call allocatefBlock1(T,Adim)
		return
	end subroutine
	subroutine allocatefBlock3(T,element)
		class(fTensor),intent(inout)::T
		integer,intent(in)::element
		integer,allocatable::Adim(:)
		if(element.gt.T%getTotalData())Then
			write(*,*)"Index is larger than the len of SymTensor"
			call error_stop()
		end if
		allocate(Adim(T%getRank()))
		call IndesToaddress(T%dim(),Adim,element)
		call allocatefBlock1(T,Adim)
		return
	end subroutine
	
	subroutine set_Symelement_vec(T,vec,element) 
		class(fTensor),intent(inout)::T
		integer,intent(in)::vec(:)
		type(Tensor),intent(in)::element
		logical::goon
		integer::i
		character*500,allocatable::w
		real*4,allocatable::qnum(:)
		goon=rule(T%SymDimension,vec,T%getRank())
		if(goon)then
			call T%SymTensor%setValue(vec,element)
			return
		end if
		call writemess('Can not set the value to the element of a fTensor')
		allocate(w)
		w='The indices are: ('
		do i=1,size(vec)-1
			w=w+vec(i)+','
		end do
		w=w+vec(i)+')'
		call writemess(w)
		allocate(qnum(size(vec)))
		qnum=T%QNDim(vec)
		w='The quantum number are: ('
		do i=1,size(vec)-1
			w=w+qnum(i)+','
		end do
		w=w+qnum(i)+')'
		call writemess(w)
		call error_stop()
	end 	subroutine
	subroutine set_Symelement_val(T,vec,element) 
		class(fTensor),intent(inout)::T
		integer,intent(in)::vec
		type(Tensor),intent(in)::element
		integer,allocatable::Adim(:)
		if(vec.gt.T%getTotalData())Then
			write(*,*)"Index is larger than the len of SymTensor"
			call error_stop()
		end if
		allocate(Adim(T%getRank()))
		call IndesToaddress(T%dim(),Adim,vec)
		call set_Symelement_vec(T,Adim,element) 
		return
	end subroutine
	
	subroutine set_random_fTensor_all(T,region)
		class(fTensor),intent(inout)::T
		real*8,optional,intent(in)::region(2)
		integer::I
		type(Tensor)::elementindex
		elementindex=NonZeroElement(T%symDimension)
		if(.not.elementindex%getFlag())then
			call T%diminfo()
			call writemess('There is no non-zero block')
			call error_stop
		end if
		call T%SymTensor%allocateBlock(elementindex)
		call T%SymTensor%random(region)
		return
	end subroutine
	subroutine set_random_fTensor_vec(T,vec,region)
		class(fTensor),intent(inout)::T
		integer,intent(in)::vec(:)
		real*8,optional,intent(in)::region(2)
		logical::goon
		integer::i
		character*500,allocatable::w
		real*4,allocatable::qnum(:)
		goon=rule(T%SymDimension,vec,T%getRank())
		if(goon)then
			call T%SymTensor%random(vec,region)
			return
		end if
		call writemess('Can not set random value to the element of a fTensor')
		allocate(w)
		w='The indices are: ('
		do i=1,size(vec)-1
			w=w+vec(i)+','
		end do
		w=w+vec(i)+')'
		call writemess(w)
		allocate(qnum(size(vec)))
		qnum=T%QNDim(vec)
		w='The quantum number are: ('
		do i=1,size(vec)-1
			w=w+qnum(i)+','
		end do
		w=w+qnum(i)+')'
		call writemess(w)
		call error_stop()
	end subroutine
	subroutine set_random_fTensor_val(T,ith,region)
		class(fTensor),intent(inout)::T
		integer,intent(in)::ith
		real*8,optional,intent(in)::region(2)
		integer,allocatable::Adim(:)
		if(ith.gt.T%getTotalData())Then
			write(*,*)"Index is larger than the len of SymTensor"
			call error_stop()
		end if
		allocate(Adim(T%getRank()))
		call IndesToaddress(T%dim(),Adim,ith)
		call set_random_fTensor_vec(T,Adim,region)
		return
	end subroutine
	
	
	type(Symdimension) function  getSubDim2(T,inde)
		type(fTensor),intent(in) :: T
		integer,intent(in)::inde
		getSubDim2=T%SymDimension.subdim.inde
		return
	end function
	type(Symdimension) function  getSubDim3(T,inde)
		type(fTensor),intent(in) :: T
		integer,intent(in)::inde(2)
		getSubDim3=T%SymDimension.subdim.inde
		return
	end function
	type(Symdimension) function  getSubDim4(T)
		type(fTensor),intent(in) :: T
		getSubDim4=T%SymDimension
		return
	end function
	type(Symdimension) function  getSubDim2_name(T,w)
		type(fTensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::w
		integer::inde
		inde=T%FindOrder(w)
		getSubDim2_name=T%SymDimension.subdim.inde
		return
	end function
	
	
	
	
	
	
	
	
	
	
	
	type(fTensor) function FermiPermutation(T,newOrder)
		type(fTensor),intent(in) :: T
		integer,intent(in)::newOrder(:)
		integer,allocatable ::inde(:)
		integer::lenOrder,i,j
		lenorder=size(newOrder)-1
		allocate(inde(lenorder))
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		FermiPermutation=T
		do i=lenorder,1,-1
			call permutefo_routine(FermiPermutation,inde(i))
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		return
	end function
	type(fTensor) function FermiPermutation_name(T,newOrderchar)result(Res)
		type(fTensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::newOrderchar(:)
		integer,allocatable::newOrder(:)
		integer::i
		allocate(newOrder(size(newOrderchar)))
		newOrder=T%FindOrder(newOrderchar)
		Res= FermiPermutation(T,newOrder)
		return
	end function
	type(fTensor) function FermiPermutation_Tensor(T,Order)result(Res)
		type(fTensor),intent(in) :: T
		type(Tensor),intent(in)::Order
		select case(Order%getType())
			case (1)
				Res=FermiPermutation(T,Order%ii())
			case (7)
				Res=FermiPermutation_name(T,Order%ai())
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()
		end select
		return
	end function
	subroutine permutation_Tensor_routine(T,Tenorder)
		class(fTensor),intent(inout)::T
		type(Tensor),intent(in)::Tenorder
		character(len=max_len_of_char_in_TData),allocatable::indechar(:)
		integer,allocatable::indeint(:)
		select case (Tenorder%getType())
			case (1)
				allocate(indeint(Tenorder%getTotalData()))
				indeint=Tenorder
				call permutation_routine(T,indeint)
			case (7)
				allocate(indechar(Tenorder%getTotalData()))
				indechar=Tenorder
				call permutation_name_routine(T,indechar)
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()	
		end select
		return
	end subroutine

	subroutine permutefo_routine(T,inde)
		class(fTensor),intent(inout)::T
		integer,intent(in)::inde
		call T%SymTensor%forward(inde)
		call T%external(externalForWard,Tensor(inde))
		return
	end subroutine
	
	subroutine permutefo_name_routine(T,indechar)
		class(fTensor),intent(inout)::T
		character(len=*),intent(in)::indechar
		integer::ith,sizeindechar,i,lenintindex
		character(len=50)::indexname
		integer,allocatable::charToInt(:)
		if(long_Name_logi(indechar))then
			ith=T%FindOrder(indechar)
			call T%SymTensor%forward(indechar)
			call T%external(externalForWard,Tensor(ith))
		else
			sizeindechar=T%getRank()
			allocate(charToInt(sizeindechar))
			lenintindex=0
			do i=1,sizeindechar
				indexname=T%outTensorName(i)
				if(indechar.equ.indexname)then
					lenintindex=lenintindex+1
					charToInt(lenintindex)=i
				end if
			end do
			call permutefo_vec_routine(T,charToInt(1:lenintindex))
		end if
		return
	end subroutine
	
	subroutine permutefo_vec_routine(T,vec_)
		class(fTensor),intent(inout)::T
		integer,intent(in)::vec_(:)
		integer::lenVec,i,j
		integer::vec(size(vec_))
		lenVec=size(vec_)
		vec=vec_
		do i=lenVec,1,-1
			call permutefo_routine(T,vec(i))
			do j=1,i-1
				if(vec(j).lt.vec(i))then
					vec(j)=vec(j)+1
				end if
			end do
		end do
		return
	end subroutine
	
	subroutine permutefo_vec_name_routine(T,indechar)
		class(fTensor),intent(inout)::T
		character(len=*),intent(in)::indechar(:)
		integer::i
		do i=size(indechar),1,-1
			call permutefo_name_routine(T,indechar(i))
		end do
		return
	end subroutine

	Type(fTensor) function FermiPermuteForWard(T,ith)Result(Res)
		Type(fTensor),intent(in)::T
		integer,intent(in)::ith
		Res%SymTensor=T%SymTensor.pf.ith
		call Res%external(externalForWard,Tensor(ith))
		return
	end function
	
	Type(fTensor) function FermiPermuteForWard_vec(T,vec)Result(Res)
		Type(fTensor),intent(in)::T
		integer,intent(in)::vec(:)
		integer::i
		Res%SymTensor=T%SymTensor
		call permutefo_vec_routine(Res,vec)
		return
	end function
	Type(fTensor) function FermiPermuteForWard_name(T,a)Result(Res)
		Type(fTensor),intent(in)::T
		character(len=*),intent(in)::a
		Res%SymTensor=T%SymTensor
		call permutefo_name_routine(Res,a)
		return
	end function
	Type(fTensor) function FermiPermuteForWard_name_vec(T,a)Result(Res)
		Type(fTensor),intent(in)::T
		character(len=*),intent(in)::a(:)
		integer::i,ith
		Res%SymTensor=T%SymTensor
		call permutefo_vec_name_routine(Res,a)
		return
	end function
	
	type(fTensor) function FermiPermuteForWard_Tensor(T,Tenorder)result(permutefo)
		type(fTensor),intent(in)::T
		type(Tensor),intent(in)::Tenorder
		character(len=max_len_of_char_in_TData),allocatable::indechar(:)
		integer,allocatable::indeint(:)
		select case (Tenorder%getType())
			case (1)
				allocate(indeint(Tenorder%getTotalData()))
				indeint=Tenorder
				permutefo=FermiPermuteForWard_vec(T,indeint)
			case (7)
				allocate(indechar(Tenorder%getTotalData()))
				indechar=Tenorder
				permutefo=FermiPermuteForWard_name_vec(T,indechar)
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()	
		end select
		return
	end function
	subroutine permutefo_Tensor_routine(T,Tenorder)
		class(fTensor),intent(inout)::T
		type(Tensor),intent(in)::Tenorder
		character(len=max_len_of_char_in_TData),allocatable::indechar(:)
		integer,allocatable::indeint(:)
		select case (Tenorder%getType())
			case (1)
				allocate(indeint(Tenorder%getTotalData()))
				indeint=Tenorder
				call permutefo_vec_routine(T,indeint)
			case (7)
				allocate(indechar(Tenorder%getTotalData()))
				indechar=Tenorder
				call permutefo_vec_name_routine(T,indechar)
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()	
		end select
		return
	end subroutine
	
	subroutine externalForWard(block,lenofblock,dimen,info)
		integer::lenofblock
		type(Tensor)::block(lenofblock)
		Type(Tensor)::info
		Type(SymDimension)::dimen
		integer,allocatable::indices(:),maxinde(:),mininde(:)
		logical::goon
		integer::i,indexEnd,rank
		indexEnd=info%ii(1)
		if(indexEnd.eq.1)return
		rank=dimen%getRank()
		allocate(indices(rank))
		allocate(mininde(rank))
		allocate(maxinde(rank))
		maxinde=dimen%dim()
		indices=1
		mininde=1
		goon=.true.
		i=1
		do while (goon)
			if(block(i)%getFlag())then
				if( QaunNumParity(dimen,1,indices(1) ).eq.(-1))then
					if(ifParity(dimen,indices,1,indexEnd,rank)) then
						block(i)=fermi_factor*block(i)
					end if
				end if
			end if
			i=i+1
			goon=inde_counter(indices,mininde,maxinde,1)
		end do
		return
	end subroutine
	
	
	
	subroutine permuteback_routine(T,inde)
		class(fTensor),intent(inout)::T
		integer,intent(in)::inde
		call T%SymTensor%Backward(inde)
		call T%external(externalBackWard,Tensor(inde))
		return
	end subroutine
	
	subroutine permuteback_name_routine(T,indechar)
		class(fTensor),intent(inout)::T
		character(len=*),intent(in)::indechar
		integer::inde,i,sizeinde,lenintindex,j
		character*50::indexname
		integer,allocatable::indes(:)
		if(long_Name_logi(indechar))then
			inde=T%FindOrder(indechar)
			call T%SymTensor%Backward(inde)
			call T%external(externalBackWard,Tensor(inde))
		else
			sizeinde=T%getRank()
			allocate(indes(sizeinde))
			lenintindex=0
			do i=1,sizeinde
				indexname=T%outTensorName(i)
				if(indechar.equ.indexname)then
					lenintindex=lenintindex+1
					indes(lenintindex)=i
				end if
			end do
			call permuteback_vec_routine(T,indes(1:lenintindex))
		end if
		return
	end subroutine
	
	subroutine permuteback_vec_routine(T,vec_)
		class(fTensor),intent(inout)::T
		integer,intent(in)::vec_(:)
		integer::vec(size(vec_)),lenVec,i,j
		lenVec=size(vec_)
		vec=vec_
		do i=1,lenVec
			call permuteback_routine(T,vec(i))
			do j=lenVec,i+1,-1
				if(vec(j).gt.vec(i))then
					vec(j)=vec(j)-1
				end if
			end do
		end do
		return
	end subroutine
	
	subroutine permuteback_vec_name_routine(T,indechar)
		class(fTensor),intent(inout)::T
		character(len=*),intent(in)::indechar(:)
		integer::vec(size(indechar)),rank
		vec=T%FindOrder(indechar)
		call permuteback_vec_routine(T,vec)
		return
	end subroutine
	
	
	Type(fTensor) function FermiPermuteBackWard(T,ith)Result(Res)
		class(fTensor),intent(in)::T
		integer,intent(in)::ith
		Res%SymTensor=T%SymTensor.pb.ith
		call Res%external(externalBackWard,Tensor(ith))
		return
	end function
	Type(fTensor) function FermiPermuteBackWard_name(T,a)Result(Res)
		class(fTensor),intent(in)::T
		character(len=*),intent(in)::a
		integer::ith
		Res%SymTensor=T%SymTensor
		call permuteback_name_routine(Res,a)
		return
	end function
	Type(fTensor) function FermiPermuteBackWard_name_vec(T,a)Result(Res)
		class(fTensor),intent(in)::T
		character(len=*),intent(in)::a(:)
		integer::i,ith
		Res%SymTensor=T%SymTensor
		call permuteback_vec_name_routine(Res,a)
		return
	end function
	Type(fTensor) function FermiPermuteBackWard_vec(T,vec)Result(Res)
		class(fTensor),intent(in)::T
		integer,intent(in)::vec(:)
		integer::i
		Res%SymTensor=T%SymTensor
		call permuteback_vec_routine(Res,vec)
		return
	end function
	
	type(fTensor) function FermiPermuteBackWard_Tensor(T,Tenorder)result(permutefo)
		type(fTensor),intent(in)::T
		type(Tensor),intent(in)::Tenorder
		character(len=max_len_of_char_in_TData),allocatable::indechar(:)
		integer,allocatable::indeint(:)
		select case (Tenorder%getType())
			case (1)
				allocate(indeint(Tenorder%getTotalData()))
				indeint=Tenorder
				permutefo=FermiPermuteBackWard_vec(T,indeint)
			case (7)
				allocate(indechar(Tenorder%getTotalData()))
				indechar=Tenorder
				permutefo=FermiPermuteBackWard_name_vec(T,indechar)
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()	
		end select
		return
	end function
	subroutine permuteback_Tensor_routine(T,Tenorder)
		class(fTensor),intent(inout)::T
		type(Tensor),intent(in)::Tenorder
		character(len=max_len_of_char_in_TData),allocatable::indechar(:)
		integer,allocatable::indeint(:)
		select case (Tenorder%getType())
			case (1)
				allocate(indeint(Tenorder%getTotalData()))
				indeint=Tenorder
				call permuteback_vec_routine(T,indeint)
			case (7)
				allocate(indechar(Tenorder%getTotalData()))
				indechar=Tenorder
				call permuteback_vec_name_routine(T,indechar)
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()	
		end select
		return
	end subroutine
	subroutine externalBackWard(block,lenofblock,dimen,info)
		integer::lenofblock
		type(Tensor)::block(lenofblock),info
		Type(SymDimension)::dimen
		integer,allocatable::indices(:),maxinde(:),mininde(:)
		logical::goon
		integer::i,rank,indexStart
		indexStart=info%ii(1)
		rank=dimen%getRank()
		if(indexStart.eq.rank)return
		allocate(indices(rank))
		allocate(mininde(rank))
		allocate(maxinde(rank))
		maxinde=dimen%dim()
		indices=1
		mininde=1
		goon=.true.
		i=1
		do while (goon)
			if(block(i)%getFlag())then
				if( QaunNumParity(dimen,rank,indices(rank)).eq.(-1))then
					if(ifParity(dimen,indices,indexStart,rank,rank))then
						block(i)=fermi_factor*block(i)
					end if
				end if
			end if
			i=i+1
			goon=inde_counter(indices,mininde,maxinde,1)
		end do
		return
	end subroutine
	
	subroutine permutation_routine(T,newOrder)
		class(fTensor),intent(inout) :: T
		integer,intent(in)::newOrder(:)
		integer,allocatable ::inde(:)
		integer::lenOrder,i,j
		lenorder=size(newOrder)-1
		allocate(inde(lenorder))
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		do i=lenorder,1,-1
			 call permutefo_routine(T,inde(i))
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		return
	end subroutine
	subroutine permutation_name_routine(T,newOrderchar)
		class(fTensor),intent(inout) :: T
		CHARACTER(len=*),intent(in)::newOrderchar(:)
		integer,allocatable::newOrder(:)
		integer,allocatable ::inde(:)
		integer::lenOrder,i,j
		allocate(newOrder(size(newOrderchar)))
		newOrder=T%FindOrder(newOrderchar)
		lenorder=size(newOrder)-1
		allocate(inde(lenorder))
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		do i=lenorder,1,-1
			call permutefo_routine(T,inde(i))
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		return
	end subroutine
	
	
!************************************************************************
	type(fTensor) function ProductTensor (T1,T2)
		type(fTensor),intent(in) :: T1,T2
		call ProductTensor%ProductTensorRoutine(T1,T2)
		return
	end function
	
	type(fTensor) function divide_real8(T1,num) result(Res)
		type(fTensor),intent(in) :: T1
		real(kind=8),intent(in) ::   num
		Res%SymTensor=T1%SymTensor/num
		return
	end function
	type(fTensor) function divide_real4(T1,num) result(Res)
		type(fTensor),intent(in) :: T1
		real(kind=4),intent(in) ::   num
		Res%SymTensor=T1%SymTensor/num
		return
	end function
	type(fTensor) function divide_Tensor(T1,num) result(Res)
		type(fTensor),intent(in) :: T1
		type(Tensor),intent(in) ::   num
		Res%SymTensor=T1%SymTensor/num
		return
	end function
	
	type(fTensor) function add_fTensor(T1,T2) result(Res)
		type(fTensor),intent(in) :: T1,T2
		Res%SymTensor=T1%SymTensor+T2%SymTensor
		return
	end function
	
	type(fTensor) function minu_fTensor(T1,T2) result(Res)
		type(fTensor),intent(in) :: T1,T2
		Res%SymTensor=T1%SymTensor-T2%SymTensor
		return
	end function
	
	type(fTensor) function multiply_number_real8(T1,num) result(Res)
		type(fTensor),intent(in) :: T1
		real(kind=8),intent(in) ::   num
		Res%SymTensor=T1%SymTensor*num
		return
	end function
	type(fTensor) function multiply_number_real4(T1,num) result(Res)
		type(fTensor),intent(in) :: T1
		real(kind=4),intent(in) ::   num
		Res%SymTensor=T1%SymTensor*num
		return
	end function
	type(fTensor) function multiply_number_com4(T1,num) result(Res)
		type(fTensor),intent(in) :: T1
		complex(kind=4),intent(in) ::   num
		Res%SymTensor=T1%SymTensor*num
		return
	end function



!**************************************************************************************************************
!**************************************************************************************************************
!
!                                  contract
!
!**************************************************************************************************************
!**************************************************************************************************************	



! T(:,col)=P(col)*T(:,col) ,input Tensor(:,:)
!where P(col)=+-1 , the parity of col

	subroutine contract_sgin_order_subroutine(T,LD1,LD2,col)
		integer,intent(in)::LD1,LD2,col
		type(Tensor),intent(inout) :: T(LD1,LD2)
		integer::i
		do i=1,LD1
			if(T(i,col)%getFlag())T(i,col)=fermi_factor*T(i,col)
		end do
		return
	end subroutine

! T(:,ith,:)=P(ith)*T(:,ith,:) ,input Tensor(:,:,:)
!where P(col)=+-1 , the parity of col

	subroutine contract_sgin_order_3dim_subroutine(T,LD1,LD2,LD3,ith)
		integer,intent(in)::LD1,LD2,LD3,ith
		type(Tensor),intent(inout) :: T(LD1,LD2,LD3)
		integer::i,j
		do j=1,LD3
			do i=1,LD1
				if(T(i,ith,j)%getFlag())T(i,ith,j)=fermi_factor*T(i,ith,j)
			end do
		end do
		return
	end subroutine

! T(row,:)=P(row)*T(row,:) ,input Tensor(:,:)
!where P(row)=+-1 , the parity of row

	subroutine contract_sgin_order_row_subroutine(T,LD1,LD2,row)
		integer,intent(in)::LD1,LD2,row
		type(Tensor),intent(inout) :: T(LD1,LD2)
		integer::i
		do i=1,LD2
			if(T(row,i)%getFlag())T(row,i)=fermi_factor*T(row,i)
		end do
		return
	end subroutine
!
! T(:,col)=P(col)*T(:,col) ,input fTensor
!where P(col)=+-1 , the parity of col

	subroutine externalContractSginOrder(block,lenofblock,dimen)
		integer::lenofblock
		type(Tensor)::block(lenofblock)
		Type(SymDimension)::dimen
		integer::col,LD1,LD2,i,rank
		rank=dimen%getRank()
		LD2=dimen%dim(rank)
		LD1=1
		do i=1,rank-1
			LD1=LD1*dimen%dim(i)
		end do
		do col=1,LD2
			if(QaunNumParity(dimen,rank,col).eq.(-1))then
				call contract_sgin_order_subroutine(block,LD1,LD2,col)
			end if
		end do
		return
	end subroutine
! T(row,:)=P(row)*T(row,:) ,input fTensor
!where P(row)=+-1 , the parity of row

	subroutine externalContractSginOrderRow(block,lenofblock,dimen)
		integer::lenofblock
		type(Tensor)::block(lenofblock)
		Type(SymDimension)::dimen
		integer::row,LD1,LD2,i,rank
		rank=dimen%getRank()
		LD1=dimen%dim(1)
		LD2=1
		do i=2,rank
			LD2=LD2*dimen%dim(i)
		end do
		do row=1,LD1
			if(QaunNumParity(dimen,1,row).eq.(-1))then
				call contract_sgin_order_row_subroutine(block,LD1,LD2,row)
			end if
		end do
		return
	end subroutine
!
! T(:,ith,:)=P(ith)*T(:,ith,:) ,input fTensor
!where P(ith)=+-1 , the parity of ith
!routine(T%Block,T%getTotalData(),T%SymDimension,T2)

	subroutine externalContractSginOrderDim(block,lenofblock,dimen,Info)
		integer::lenofblock
		type(Tensor)::block(lenofblock)
		Type(SymDimension)::dimen
		type(Tensor)::Info
		integer::ith,LD1,LD2,LD3,i,rank,ithleg
		ithleg=Info%ii(1)
		rank=dimen%getRank()
		if(ithleg.eq.1)then
			call externalContractSginOrderRow(block,lenofblock,dimen)
			return
		end if
		if(ithleg.eq.rank)then
			call externalContractSginOrder(block,lenofblock,dimen)
			return
		end if
		LD1=dimen%dim(1)
		do i=2,ithleg-1
			LD1=LD1*dimen%dim(i)
		end do
		LD2=dimen%dim(ithleg)
		LD3=dimen%dim(ithleg+1)
		do i=ithleg+2,rank
			LD3=LD3*dimen%dim(i)
		end do

		do ith=1,dimen%dim(ithleg)
			if(QaunNumParity(dimen,ithleg,ith).eq.(-1))then
				call contract_sgin_order_3dim_subroutine(block,LD1,LD2,LD3,ith)
			end if
		end do
		return
	end subroutine	
!
! [T1]---->---[T2]
! change to
! [T1]----<---[T2]
!
!modify the rule of last leg of T1 and the first leg of T2
!
!	

	subroutine Reverse_Fermi_Rule1(T1,T2)
		type(fTensor),intent(inout)::T1,T2
		integer::rank1
		if(T1%getTotalBlock().gt.T2%getTotalBlock())then
			call T2%external(externalContractSginOrderRow)
		else
			call T1%external(externalContractSginOrder)
		end if
		rank1=T1%getRank()
		call T1%setFermiArrow(rank1,-1*T1%getFermiArrow(rank1))
		call T2%setFermiArrow(1,-1*T2%getFermiArrow(1))
		return
	end subroutine
	
	subroutine Reverse_Fermi_Rule_specify(T1,T2,row)
		character*1,intent(in)::row
		type(fTensor),intent(inout)::T1,T2
		integer::rank1
		if(row.eq.'r')then
			call T2%external(externalContractSginOrderRow)
		else if(row.eq.'c')then
			call T1%external(externalContractSginOrder)
		else
			call writemess('ERROR in Reverse_Fermi_Rule, input parameter,row='+row,-1)
			call writemess('row="r", change the row of T2',-1)
			call writemess('row="c", change the col of T1',-1)
			call error_stop()
		end if
		rank1=T1%getRank()
		call T1%setFermiArrow(rank1,-1*T1%getFermiArrow(rank1))
		call T2%setFermiArrow(1,-1*T2%getFermiArrow(1))
		return
	end subroutine
	
	subroutine Reverse_Fermi_Rule_specify2(T1,ith1,T2,ith2)
		integer,intent(in)::ith1,ith2
		type(fTensor),intent(inout)::T1,T2
		if(T1%getTotalBlock().gt.T2%getTotalBlock())then
			call T2%external(externalContractSginOrderDim,Tensor((/ith2/)))
		else
			call T1%external(externalContractSginOrderDim,Tensor((/ith1/)))
		end if
		call T1%setFermiArrow(ith1,-1*T1%getFermiArrow(ith1))
		call T2%setFermiArrow(ith2,-1*T2%getFermiArrow(ith2))
		return
	end subroutine
	subroutine Reverse_Fermi_Rule_specify3(T1,ith1)
		integer,intent(in)::ith1
		class(fTensor),intent(inout)::T1
		call T1%external(externalContractSginOrderDim,Tensor((/ith1/)))
		call T1%setFermiArrow(ith1,-1*T1%getFermiArrow(ith1))
		return
	end subroutine


	subroutine monitor_check(rule1,rule2)
		integer,intent(in)::rule1,rule2
		if(mointer_order_flag)then
			if(rule1*rule2.ge.0)then
				call writemess('ERROR in fermi-arrow',-1)
				write(*,*)rule1,rule2
				call error_stop
			end if
		else
			if(rule1*rule2.gt.0)then
				call writemess('ERROR in fermi-arrow',-1)
					write(*,*)rule1,rule2
				call error_stop
			end if
		end if
		return
	end subroutine
	
	subroutine un_set_mointer_order_flag()
		mointer_order_flag=.false.
	end subroutine
	subroutine set_mointer_order_flag()
		mointer_order_flag=.true.
	end subroutine

!fermi-arrow=1     mean the leg is <a| :  --->----[A]
!fermi-arrow=-1 	mean the leg is |b> :  ---<----[B]
!
!   the contraction should be    <a|b>
!   if not do the Correction
!      __         __
!     |  |       |  |
!     |A |--->---|B |
!     |  |--->---|  | 
!     |  |---<---|  |
!     |__|       |__|
!
!the fermi-arrow of A are   -1,  -1,  1
!the fermi-arrow of B are    1,   1, -1
!The first two leg will not cause the SignCorrect
! the last leg, fermi-arrow of A is 1 and fermi-arrow of B is -1, it cause SignCorrect
!
!
!If fermi-arrow=0 , will not monitor the order

	subroutine RuleSignCorrect(T1,T1_,T1Name,T2,T2_,T2Name,lenName)
		type(fTensor),intent(inout)::T1,T2
		type(fTensor),intent(in)::T1_,T2_
		character(len=*),intent(in)::T2Name(:),T1Name(:)
		integer,intent(in)::lenName
		integer::i,rank
		integer::rule1,rule2,SymRule1,SymRule2
		T2=T2_.pf.T2Name(1:lenName)
		rank=T1_%getRank()
		call T1%empty()
		do i=lenName,1,-1
			if(T1%getFlag())then
				call T1%backWard(T1Name(i))
			else
				T1=T1_.pb.T1Name(i)
			end if
			rule1=T1%getFermiArrow(rank)
			rule2=T2%getFermiArrow(i)
			call monitor_check(rule1,rule2)
			SymRule1=T1%getRule(rank)
			SymRule2=T2%getRule(i)
			call checkSymmetryRule(SymRule1,SymRule2,T1%getName(rank),T2%getName(i))
			if((rule2.lt.0).and.(rule1.gt.0))then
				call T1%external(externalContractSginOrder)
			end if
		end do
		call T1%SymTensor%backWard(T1Name(1:lenName))
		return
	end subroutine
	subroutine RuleSignCorrectInout1(T1,T1_,T1Name,T2,T2Name,lenName)!the subroutine will change the order of legs in T2, will not change T1_
		type(fTensor),intent(inout)::T1,T2
		type(fTensor),intent(in)::T1_
		character(len=*),intent(in)::T2Name(:),T1Name(:)
		integer,intent(in)::lenName
		integer::i,rank
		integer::rule1,rule2,SymRule1,SymRule2
		call T2%forward(T2Name(1:lenName))
		rank=T1_%getRank()
		call T1%empty()
		do i=lenName,1,-1
			if(T1%getFlag())then
				call T1%backWard(T1Name(i))
			else
				T1=T1_.pb.T1Name(i)
			end if
			rule1=T1%getFermiArrow(rank)
			rule2=T2%getFermiArrow(i)
			call monitor_check(rule1,rule2)
			SymRule1=T1%getRule(rank)
			SymRule2=T2%getRule(i)
			call checkSymmetryRule(SymRule1,SymRule2,T1%getName(rank),T2%getName(i))
			if((rule2.lt.0).and.(rule1.gt.0))then
				call T1%external(externalContractSginOrder)
			end if
		end do
		call T1%SymTensor%backWard(T1Name(1:lenName))
		return
	end subroutine
	subroutine RuleSignCorrectInout2(T1,T1Name,T2,T2Name,lenName)!the subroutine will change the order of legs in T2, and destroy T1
		type(fTensor),intent(inout)::T1,T2
		character(len=*),intent(in)::T2Name(:),T1Name(:)
		integer,intent(in)::lenName
		integer::i,rank
		integer::rule1,rule2,SymRule1,SymRule2
		call T2%forward(T2Name(1:lenName))
		rank=T1%getRank()
		do i=lenName,1,-1
			call T1%backWard(T1Name(i))
			rule1=T1%getFermiArrow(rank)
			rule2=T2%getFermiArrow(i)
			call monitor_check(rule1,rule2)
			SymRule1=T1%getRule(rank)
			SymRule2=T2%getRule(i)
			call checkSymmetryRule(SymRule1,SymRule2,T1%getName(rank),T2%getName(i))
			if((rule2.lt.0).and.(rule1.gt.0))then
				call T1%external(externalContractSginOrder)
			end if
		end do
		call T1%SymTensor%backWard(T1Name(1:lenName))
		return
	end subroutine
	subroutine RuleSignCorrectInout3(T1,T1Name,T2,T2Name,lenName)!the subroutine will change the order of legs in T1, and destroy T2
		type(fTensor),intent(inout)::T1,T2
		character(len=*),intent(in)::T2Name(:),T1Name(:)
		integer,intent(in)::lenName
		integer::i,rank,firsti
		integer::rule1,rule2,SymRule1,SymRule2
		call T1%backward(T1Name(1:lenName))
		rank=T1%getRank()
		firsti=rank-lenName
		do i=1,lenName
			call T2%forWard(T2Name(i))
			rule1=T1%getFermiArrow(firsti+i)
			rule2=T2%getFermiArrow(1)
			call monitor_check(rule1,rule2)
			SymRule1=T1%getRule(firsti+i)
			SymRule2=T2%getRule(1)
			call checkSymmetryRule(SymRule1,SymRule2,T1%getName(firsti+i),T2%getName(1))
			if((rule2.lt.0).and.(rule1.gt.0))then
				call T2%external(externalContractSginOrderRow)
			end if
		end do
		call T2%SymTensor%forWard(T2Name(1:lenName))
		return
	end subroutine
	subroutine RuleSignCorrect_int(T1,T1_,T1index_,T2,T2_,T2index,lenindex)
		type(fTensor),intent(inout)::T1,T2
		type(fTensor),intent(in)::T1_,T2_
		integer,intent(in)::T1index_(:),T2index(:)
		integer,intent(in)::lenindex
		integer::i,rank,indices(lenindex),j,T1index(lenindex)
		integer::rule1,rule2,SymRule1,SymRule2
		T2=T2_.pf.T2index(1:lenindex)
		rank=T1_%getRank()
		call T1%empty()
		do i=1,lenindex
			indices(i)=rank-i+1
		end do
		T1index=T1index_(1:lenindex)
		do i=lenindex,1,-1
			if(T1%getFlag())then
				call T1%backWard(T1index(i))
			else
				T1=T1_.pb.T1index(i)
			end if
			do j=i-1,1,-1
				if(T1index(j).gt.T1index(i))then
					T1index(j)=T1index(j)-1
				end if
			end do
			rule1=T1%getFermiArrow(rank)
			rule2=T2%getFermiArrow(i)
			call monitor_check(rule1,rule2)
			SymRule1=T1%getRule(rank)
			SymRule2=T2%getRule(i)
			call checkSymmetryRule(SymRule1,SymRule2,T1%getName(rank),T2%getName(i))
			if((rule2.lt.0).and.(rule1.gt.0))then
				call T1%external(externalContractSginOrder)
			end if
		end do
		call T1%SymTensor%backWard(indices)
		return
	end subroutine
	
	subroutine RuleSignCorrect2(T1,T1_,T1Name,T2,T2_,T2Name)
		type(fTensor),intent(inout)::T1,T2
		type(fTensor),intent(in)::T1_,T2_
		character(len=*),intent(in)::T2Name,T1Name
		integer::i,rank
		integer::rule1,rule2,SymRule1,SymRule2
		T1=T1_.pb.T1Name
		T2=T2_.pf.T2Name
		rule1=T1%getFermiArrow(T1%getrank())
		rule2=T2%getFermiArrow(1)
		call monitor_check(rule1,rule2)
		SymRule1=T1%getRule(T1%getrank())
		SymRule2=T2%getRule(1)
		call checkSymmetryRule(SymRule1,SymRule2,T1%getName(T1%getrank()),T2%getName(1))
		if((rule2.lt.0).and.(rule1.gt.0))then
			call T1%external(externalContractSginOrder)
		end if
		return
	end subroutine
	
	subroutine RuleSignCorrect_int2(T1,T1_,T1index,T2,T2_,T2index)
		type(fTensor),intent(inout)::T1,T2
		type(fTensor),intent(in)::T1_,T2_
		integer,intent(in)::T1index,T2index
		integer::i,rank
		integer::rule1,rule2,SymRule1,SymRule2
		T1=T1_.pb.T1index
		T2=T2_.pf.T2index
		rule1=T1%getFermiArrow(T1%getrank())
		rule2=T2%getFermiArrow(1)
		call monitor_check(rule1,rule2)
		SymRule1=T1%getRule(T1%getrank())
		SymRule2=T2%getRule(1)
		call checkSymmetryRule(SymRule1,SymRule2,T1%getName(T1%getrank()),T2%getName(1))
		if((rule2.lt.0).and.(rule1.gt.0))then
			call T1%external(externalContractSginOrder)
		end if
		return
	end subroutine
	
!******************  contract  *********************
!	T1:[i1,i2,i3,i4,i5,i6,i7,i8]
!	T2:[j1,j2,j3,j4,j5,j6,j7,j8,j9,j10]
!	i1=(/5,1,2/)
!	i2=(/10,3,5/)
!	then the result will be T1'*T2'
!	where ,
!	T1'=[i3,i4,i6,i7,i8,(i5*i1*i2)]
!	T2'=[(j10*j3*j5),j1,j2,j4,j6,j7,j8,j9]
!  The sign of T1' will be determine by T1, T1 permute to (/2,1,5/) with fermi rule **********NOTE here not (/5,1,2/) *************
!         And the permute to [i3,i4,i6,i7,i8,(i5*i1*i2)] non-fermi rule
!  T2' will be determine by T1, T1 permute to (/10,3,5/)  with fermi rule 
!  example   
!       C1*C2 |11> = C1*C2 * C1^+ * C2^+ |0> = -C2*C1 * C1^+ * C2^+ |0> = -|00>
!     The order of C1*C2 shoud reoder as C2*C1
!     but the data store in memery is C1*C2
!     so the function permute to C2*C1 with fermi rule, and then permute to C1*C2 with non-fermi rule
!
! 	input Tensor should be in its original dimenison,there is no contract on it
!	if present len_of_contract, len_of_contract specify the length of  i1, and i2

	type(fTensor) function contract_noName(T1_,i1,T2_,i2,len_of_contract) result(T)
		type(fTensor),intent(in) :: T1_,T2_
		integer,intent(in) :: i1(:),i2(:)
		integer,optional,intent(in)::len_of_contract(2)
		type(fTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		rank1=T1_%getrank()
		rank2=T2_%getrank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(i1))
			leni2=min(len_of_contract(2),size(i2))
		else
			leni1=size(i1)
			leni2=size(i2)
		end if
		call RuleSignCorrect_int(T1,T1_,i1,T2,T2_,i2,leni1)
		call T1%fuse(rank1-leni1+1,rank1)
		call T2%fuse(1,leni2)
		!T1=T1%FermiFuse(rank1-leni1+1,'a',.false.)
		!T2=T2%FermiFuse(leni2,.true.)
		T=T1 * T2
		return
	end function
	type(fTensor) function contract_noName2(T1_,i1,T2_,i2) result(T)
		type(fTensor),intent(in) :: T1_,T2_
		type(fTensor) :: T1,T2
		integer,intent(in) :: i1,i2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		call RuleSignCorrect_int2(T1,T1_,i1,T2,T2_,i2)
		T=T1*T2
		return
	end function
	type(fTensor) function contract_name(T1_,name1,T2_,name2,len_of_contract) result(T)
		type(fTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:),name2(:)
		integer,optional,intent(in)::len_of_contract(2)
		type(fTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.(T1_%if_original_dim().and.T2_%if_original_dim())) then
			write(*,*)"ERROR in contract"
			write(*,*)"stop"
			call error_stop()
		end if
		rank1=T1_%getrank()
		rank2=T2_%getrank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(name1))
			leni2=min(len_of_contract(2),size(name2))
		else
			leni1=size(name1)
			leni2=size(name2)
		end if
		call RuleSignCorrect(T1,T1_,name1,T2,T2_,name2,leni1)
		call T1%fuse(rank1-leni1+1,rank1)
		call T2%fuse(1,leni2)
		!T1=T1%FermiFuse(rank1-leni1+1,'a',.false.)
		!T2=T2%FermiFuse(leni2,.true.)
		T=T1 * T2
		return
	end function
	type(fTensor) function contract_name2(T1_,name1,T2_,name2) result(T)
		type(fTensor),intent(in) :: T1_,T2_
		type(fTensor) :: T1,T2
		character(len=*),intent(in)::name1,name2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.(T1_%if_original_dim().and.T2_%if_original_dim())) then
			write(*,*)"ERROR in contract"
			write(*,*)"stop"
			call error_stop()
		end if
		call RuleSignCorrect2(T1,T1_,name1,T2,T2_,name2)
		T=T1*T2
		return
	end function

!******************  contract the same names  *********************	

	subroutine find_same_name(inname1,inname2,SameName,lenofname)
		character(len=*),intent(in)::inname1(:),inname2(:)
		character(len=*),intent(inout)::SameName(:)
		integer,intent(inout)::lenofname
		integer::i,j,k
		k=1
		lenofname=0
		do i=1,size(inname1)
			do j=1,size(inname2)
				if(inname1(i).equ.inname2(j))then
					SameName(k)=inname1(i)
					lenofname=lenofname+1
					k=k+1
				end if
			end do
		end do
		return
	end subroutine

	type(fTensor) function contract_Same_name(T1_,T2_) result(T)
		type(fTensor),intent(in) :: T1_,T2_
		type(fTensor)::T1,T2
		character(len=len_of_Name+len_of_Name),allocatable::Samename(:),name1(:),name2(:)
		integer::lenofname,rank1,rank2,i
		integer,allocatable::fermiindex(:),contractindex(:)
		call writemess('Do not finished this part yet, (contract)',-1)
		call error_stop()
		
		
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.(T1_%if_original_dim().and.T2_%if_original_dim())) then
			write(*,*)"ERROR in contract"
			write(*,*)"stop"
			call error_stop()
		end if
		rank1=T1_%getRank()
		rank2=T2_%getRank()
		allocate(name1(rank1))
		allocate(name2(rank2))
		allocate(Samename(max(rank1,rank2) ))
		name1=.iname.T1_
		name2=.iname.T2_
		call find_same_name(name1,name2,SameName,lenofname) 
		allocate(fermiindex(lenofname))
		allocate(contractindex(lenofname))
		do i=1,lenofname
			fermiindex(i)=T1_%FindOrder(SameName(i))
			contractindex(i)=rank1-i+1
		end do
		T1=T1_.pb.fermiindex
		call T1%SymTensor%backward(contractindex)
		!T1=T1_.pb.SameName(1:lenofname)
		call T1%fuse(rank1-lenofname+1,rank1)
		T2=T2_.pf.SameName(1:lenofname)
		call T2%fuse(1,lenofname)
		T=T1 * T2
		return
	end function
	type(fTensor) function contract_name_int(T1_,name1,T2_,i2,len_of_contract) result(T)
		type(fTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:)
		integer,intent(in)::i2(:)
		integer :: i1(size(name1))
		integer,optional,intent(in)::len_of_contract(2)
		type(fTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2,i
		integer,allocatable::fermiindex(:),contractindex(:)
		call writemess('Do not finished this part yet, (contract)',-1)
		call error_stop()
		
		
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T1_%if_original_dim()) then
			write(*,*)"ERROR in contract"
			write(*,*)"stop"
			call error_stop()
		end if
		i1=T1_%FindOrder(name1)
		rank1=T1_%getrank()
		rank2=T2_%getrank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(i1))
			leni2=min(len_of_contract(2),size(i2))
		else
			leni1=size(i1)
			leni2=size(i2)
		end if
		allocate(fermiindex(leni1))
		allocate(contractindex(leni1))
		do i=1,leni1
			fermiindex(i)=i1(leni1-i+1)
			contractindex(i)=rank1-i+1
		end do
		T1=T1_.pb.fermiindex
		call T1%SymTensor%backward(contractindex)
	!	T1=T1_.pb.i1(1:leni1)
		call T1%fuse(rank1-leni1+1,rank1)
		T2=T2_.pf.i2(1:leni2)
		call T2%fuse(1,leni2)
		T=T1 * T2
		return
	end function
	type(fTensor) function contract_name_int2(T1_,name1,T2_,i2) result(T)
		type(fTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1
		integer,intent(in)::i2
		call writemess('Do not finished this part yet, (contract)',-1)
		call error_stop()
		
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T1_%if_original_dim()) then
			write(*,*)"ERROR in contract"
			write(*,*)"stop"
			call error_stop()
		end if
		T = (T1_.pb.name1) * (T2_.pf.i2)
		return
	end function

	type(fTensor) function contract_int_name(T1_,i1,T2_,name2,len_of_contract) result(T)
		type(fTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name2(:)
		integer,intent(in)::i1(:)
		integer,optional,intent(in)::len_of_contract(2)
		integer :: i2(size(name2))
		type(fTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2,i
		integer,allocatable::fermiindex(:),contractindex(:)
		call writemess('Do not finished this part yet, (contract)',-1)
		call error_stop()
		
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%if_original_dim()) then
			write(*,*)"ERROR in contract"
			write(*,*)"stop"
			call error_stop()
		end if
		i2=T2_%FindOrder(name2)
		rank1=T1_%getrank()
		rank2=T2_%getrank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(i1))
			leni2=min(len_of_contract(2),size(i2))
		else
			leni1=size(i1)
			leni2=size(i2)
		end if
		allocate(fermiindex(leni1))
		allocate(contractindex(leni1))
		do i=1,leni1
			fermiindex(i)=i1(leni1-i+1)
			contractindex(i)=rank1-i+1
		end do
		T1=T1_.pb.fermiindex
		call T1%SymTensor%backward(contractindex)
		!T1=T1_.pb.i1(1:leni1)
		call T1%fuse(rank1-leni1+1,rank1)
		T2=T2_.pf.i2(1:leni2)
		call T2%fuse(1,leni2)
		T=T1 * T2
		return
	end function
	type(fTensor) function contract_int_name2(T1_,i1,T2_,name2) result(T)
		type(fTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name2
		integer,intent(in)::i1
		
		call writemess('Do not finished this part yet, (contract)',-1)
		call error_stop()
		
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T2_%if_original_dim()) then
			write(*,*)"ERROR in contract"
			write(*,*)"stop"
			call error_stop()
		end if
		T = (T1_.pb.i1) * (T2_.pf.name2)
		return
	end function
	
	
	subroutine contract_name_routine(T,T1,name1,T2,name2,len_of_contract) 
		class(fTensor),target::T
		class(SymTensor),target:: T1,T2
		character(len=*),intent(in)::name1(:),name2(:)
		integer :: i1(size(name1)),i2(size(name2))
		integer,optional,intent(in)::len_of_contract(2)
		integer::leni1,leni2,rank1,rank2
		class(fTensor),pointer::pT,pT1,pT2
		if(.not.T1%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(if_original_dim(T1%SymDimension).and.if_original_dim(T2%SymDimension))) then
			call writemess("ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		select type(T1)
			type is (fTensor)
				select type(T2)
					type is (fTensor)
						if(old_verion_gcc_flag)then
							call print_contract_message()
							T=contract(T1,name1,T2,name2,len_of_contract) 
							return
						end if
						pT=>T
						pT1=>T1
						pT2=>T2
						if(associated(pT,pT1).or.associated(pT,pT2).or.associated(pT1,pT2))then
							call writemess('input Tensors can not be the same variable',-1)
							call writemess('error in call T%contract(A,[names],B,[names])')
							call writemess('T, A and B can not be a same variable')
							call error_stop
						end if
						pT=>null()
						pT1=>null()
						pT2=>null()

						i1=T1%FindOrder(name1)
						i2=T2%FindOrder(name2)
						rank1=T1%Getrank()
						rank2=T2%Getrank()
						if(present(len_of_contract))then
							leni1=min(len_of_contract(1),size(name1))
							leni2=min(len_of_contract(2),size(name2))
						else
							leni1=size(i1)
							leni2=size(i2)
						end if
						call RuleSignCorrectinout1(T,T1,name1,T2,name2,leni1)
						call T%fuse(rank1-leni1+1,rank1)
						call T2%fuse(1,leni2)
						T=T*T2  
						call T1%split( )
						call T2%split( )
						return
					class default
						call writemess('ERRROR type in fTensor%contract')
						call error_stop
				end select
			class default
				call writemess('ERRROR type in fTensor%contract')
				call error_stop
		end select
	end subroutine
	
	subroutine contract_name_routine1(T,name1,T2,name2,len_of_contract) 
		class(fTensor),target::T
		class(SymTensor),target:: T2
		character(len=*),intent(in)::name1(:),name2(:)
		integer,optional,intent(in)::len_of_contract(2)
		integer::leni1,leni2,rank1,rank2
		class(fTensor),pointer::pT,pT2
		if(.not.T%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(if_original_dim(T%SymDimension).and.if_original_dim(T2%SymDimension))) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		select type(T2)
			type is (fTensor)
				if(old_verion_gcc_flag)then
					call print_contract_message()
					T=contract(T,name1,T2,name2,len_of_contract) 
					return
				end if

				pT=>T
				pT2=>T2
				if(associated(pT,pT2))then
					call writemess('input Tensors can not be the same variable',-1)
					call writemess('error in call T%contract([names],B,[names])')
					call writemess('T and B can not be a same variable')
					call error_stop
				end if
				pT=>null()
				pT2=>null()
				rank1=T%getRank()
				rank2=T2%getRank()
				if(present(len_of_contract))then
					leni1=min(len_of_contract(1),size(name1))
					leni2=min(len_of_contract(2),size(name2))
				else
					leni1=size(name1)
					leni2=size(name2)
				end if
				
				call RuleSignCorrectinout2(T,name1,T2,name2,leni1)
				call T%fuse(rank1-leni1+1,rank1)
				call T2%fuse(1,leni2)
				T=T*T2
				call T2%split()
				return
			class default
				call writemess('ERRROR type in fTensor%contract')
				call error_stop
		end select
	end subroutine
	subroutine contract_name_routine2(T,T1,name1,name2,len_of_contract) 
		class(fTensor),target::T
		class(SymTensor),target:: T1
		character(len=*),intent(in)::name1(:),name2(:)
		integer,optional,intent(in)::len_of_contract(2)
		integer::leni1,leni2,rank1,rank2
		class(fTensor),pointer::pT,pT1
		if(.not.T1%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(if_original_dim(T1%SymDimension).and.if_original_dim(T%SymDimension))) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		select type(T1)
			type is (fTensor)
				if(old_verion_gcc_flag)then
					call print_contract_message()
					T=contract(T1,name1,T,name2,len_of_contract) 
					return
				end if

				pT=>T
				pT1=>T1
				if(associated(pT,pT1))then
					call writemess('input Tensors can not be the same variable',-1)
					call writemess('error in call T%contract(A,[names],[names])')
					call writemess('T and A can not be a same variable')
					call error_stop
				end if
				pT=>null()
				pT1=>null()
				rank1=T1%getrank()
				rank2=T%getrank()
				if(present(len_of_contract))then
					leni1=min(len_of_contract(1),size(name1))
					leni2=min(len_of_contract(2),size(name2))
				else
					leni1=size(name1)
					leni2=size(name2)
				end if
				call RuleSignCorrectinout3(T1,name1,T,name2,leni1)
				call T1%fuse(rank1-leni1+1,rank1)
				call T%fuse(1,leni2)
				T=T1*T
				call T1%split( )
				return
			class default
				call writemess('ERRROR type in fTensor%contract')
				call error_stop
		end select
	end subroutine
	subroutine contract_name_routine4(T,T1,name1,T2,name2) 
		class(fTensor),target::T
		class(SymTensor),target:: T1,T2
		character(len=*),intent(in)::name1,name2
		class(fTensor),pointer::pT,pT1,pT2
		if(.not.T1%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(if_original_dim(T1%SymDimension).and.if_original_dim(T2%SymDimension))) then
			call writemess("ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		select type(T1)
			type is (fTensor)
				select type(T2)
					type is (fTensor)
						if(old_verion_gcc_flag)then
							call print_contract_message()
							T=contract(T1,name1,T2,name2) 
							return
						end if
						pT=>T
						pT1=>T1
						pT2=>T2
						if(associated(pT,pT1).or.associated(pT,pT2).or.associated(pT1,pT2))then
							call writemess('input Tensors can not be the same variable',-1)
							call writemess('error in call T%contract(A,[names],B,[names])')
							call writemess('T, A and B can not be a same variable')
							call error_stop
						end if
						pT=>null()
						pT1=>null()
						pT2=>null()
						call RuleSignCorrectinout1(T,T1,[name1],T2,[name2],1)
						T=T*T2  
						return
					class default
						call writemess('ERRROR type in fTensor%contract')
						call error_stop
				end select
			class default
				call writemess('ERRROR type in fTensor%contract')
				call error_stop
		end select
	end subroutine
	subroutine contract_name_routine5(T,name1,T2,name2) 
		class(fTensor),target::T
		class(SymTensor),target:: T2
		character(len=*),intent(in)::name1,name2
		class(fTensor),pointer::pT,pT2
		if(.not.T%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(if_original_dim(T2%SymDimension))) then
			call writemess("ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		select type(T)
			type is (fTensor)
				select type(T2)
					type is (fTensor)
						if(old_verion_gcc_flag)then
							call print_contract_message()
							T=contract(T,name1,T2,name2)
							return 
						end if
						pT=>T
						pT2=>T2
						if(associated(pT,pT2))then
							call writemess('input Tensors can not be the same variable',-1)
							call writemess('error in call T%contract(A,[names],B,[names])')
							call writemess('T, A and B can not be a same variable')
							call error_stop
						end if
						pT=>null()
						pT2=>null()
						call RuleSignCorrectinout2(T,[name1],T2,[name2],1)
						T=T*T2  
						return
					class default
						call writemess('ERRROR type in fTensor%contract')
						call error_stop
				end select
			class default
				call writemess('ERRROR type in fTensor%contract')
				call error_stop
		end select
	end subroutine
	subroutine contract_name_routine6(T,T1,name1,name2) 
		class(fTensor),target::T
		class(SymTensor),target :: T1
		character(len=*),intent(in)::name1,name2
		class(fTensor),pointer::pT,pT1
		if(.not.T1%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.if_original_dim(T1%SymDimension)) then
			call writemess("ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		select type(T1)
			type is (fTensor)
				select type(T)
					type is (fTensor)
						if(old_verion_gcc_flag)then
							call print_contract_message()
							T=contract(T1,name1,T,name2) 
							return
						end if
						pT=>T
						pT1=>T1
						if(associated(pT,pT1))then
							call writemess('input Tensors can not be the same variable',-1)
							call writemess('error in call T%contract(name1,B,name2)')
							call writemess('T and B can not be a same variable')
							call error_stop
						end if
						pT=>null()
						pT1=>null()
						call RuleSignCorrectinout3(T1,[name1],T,[name2],1)
						T=T1*T
						return
					class default
						call writemess('ERRROR type in fTensor%contract')
						call error_stop
				end select
			class default
				call writemess('ERRROR type in fTensor%contract')
				call error_stop
		end select
	end subroutine
!*********************************************
 
	type(fTensor) function Htranspose(T)
		type(fTensor),intent(in) :: T
		integer :: rank,i,charlen
		integer,allocatable::indices(:)
		character(len=max_len_of_char_in_TData)::Tname
		rank=T%getRank()
		allocate(indices(rank))
		do i=1,rank
			indices(i)=rank-i+1
		end do
		Htranspose%SymTensor=(.con.T%SymTensor).p.indices
		do i=1,rank
			call Htranspose%setFermiArrow(i,-1*Htranspose%getFermiArrow(i))
			call Htranspose%setRule(i,-1*Htranspose%getRule(i))
		end do
		return
	end function
	
	type(fTensor) function Htranspose2(T)result(Htranspose)
		type(fTensor),intent(in) :: T
		integer :: rank,i,charlen
		integer,allocatable::indices(:)
		character(len=max_len_of_char_in_TData)::Tname
		rank=T%getRank()
		allocate(indices(rank))
		do i=1,rank
			indices(i)=rank-i+1
		end do
		Htranspose%SymTensor=(.con.T%SymTensor).p.indices
		do i=1,rank
			call Htranspose%setFermiArrow(i,-1*Htranspose%getFermiArrow(i))
			call Htranspose%setRule(i,-1*Htranspose%getRule(i))
		end do
		do i=1,rank
			Tname=Htranspose%outTensorName(i)
			charlen=len(trim(Tname))
			if(Tname(charlen:charlen).eq.dag_mark) then
				call Htranspose%setName(i,Tname(1:charlen-1))
			else
				call Htranspose%setName(i,Tname+dag_mark)
			end if
		end do
		return
	end function
	
	type(fTensor) function conjugate(T)
		type(fTensor),intent(in) :: T
		integer :: i
		conjugate%SymTensor=.con.T%SymTensor
		return
	end function
	
	
	
	
	type(fTensor) function expmfTensor(T)
		type(fTensor),intent(in) :: T
		expmfTensor%SymTensor=expm(T%SymTensor)
		return
	end function
	type(fTensor) function FermiFusefTensorRule(B,outDim,Order,NewRule,row_)result(FermiFusefTensor)
		class(fTensor),intent(in)::B
		type(SymDimension),intent(inout)::outDim
		integer,intent(in)::NewRule
		type(Tensor),intent(inout)::Order
		logical,optional,intent(in)::row_
		logical::row
		integer::i,rank
		type(QuanNum)::NewQuanNum
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		if(row)then
			outDim=B%SymDimension.subdim.[1,2]
			NewQuanNum=fuseOrder(order,B%quantumnumber(1),B%quantumnumber(2),NewRule)
			FermiFusefTensor=B%QuanFuse(NewQuanNum,Order,.true.)
		else
			rank=B%getRank()
			i=rank-1
			outDim=B%SymDimension.subdim.[i,rank]
			NewQuanNum=fuseOrder(order,B%quantumnumber(i),B%quantumnumber(rank),NewRule)
			FermiFusefTensor=B.pb.i
			call FermiFusefTensor%SymTensor%backward(i)
			FermiFusefTensor=FermiFusefTensor%QuanFuse(NewQuanNum,Order,.false.)
		end if
		return
	end function


	type(fTensor) function FermiFusefTensor(B,outDim,Order,row_)
		class(fTensor),intent(in)::B
		type(SymDimension),intent(inout)::outDim
		type(Tensor),intent(inout)::Order
		logical,optional,intent(in)::row_
		logical::row
		integer::i,rank,rule
		type(QuanNum)::NewQuanNum
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		if(row)then
			rule=1
			outDim=B%SymDimension.subdim.[1,2]
			NewQuanNum=fuseOrder(order,B%quantumnumber(1),B%quantumnumber(2),rule)
			FermiFusefTensor=B%QuanFuse(NewQuanNum,Order,.true.)
		else
			rule=-1
			rank=B%getRank()
			i=rank-1
			outDim=B%SymDimension.subdim.[i,rank]
			NewQuanNum=fuseOrder(order,B%quantumnumber(i),B%quantumnumber(rank),rule)
			FermiFusefTensor=B.pb.i
			call FermiFusefTensor%SymTensor%backward(i)
			FermiFusefTensor=FermiFusefTensor%QuanFuse(NewQuanNum,Order,.false.)
		end if
		return
	end function

	
	type(fTensor) function FermiFusefTensorRule_(B,outDim,Order,NewRule,non_fermi_fuse,row_)result(FermiFusefTensor)
		class(fTensor),intent(in)::B
		type(SymDimension),intent(inout)::outDim
		type(Tensor),intent(inout)::Order
		integer,intent(in)::NewRule
		logical,optional,intent(in)::row_
		character*1,intent(in)::non_fermi_fuse
		logical::row
		integer::i,rank
		type(QuanNum)::NewQuanNum
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		if(row)then
			outDim=B%SymDimension.subdim.[1,2]
			NewQuanNum=fuseOrder(order,B%quantumnumber(1),B%quantumnumber(2),NewRule)
			FermiFusefTensor=B%QuanFuse(NewQuanNum,Order,.true.)
		else
			rank=B%getRank()
			i=rank-1
			outDim=B%SymDimension.subdim.[i,rank]
			NewQuanNum=fuseOrder(order,B%quantumnumber(i),B%quantumnumber(rank),NewRule)
			FermiFusefTensor=B%QuanFuse(NewQuanNum,Order,.false.)
		end if
		return
	end function

	type(fTensor) function FermiFusefTensor_(B,outDim,Order,non_fermi_fuse,row_)result(FermiFusefTensor)
		class(fTensor),intent(in)::B
		type(SymDimension),intent(inout)::outDim
		type(Tensor),intent(inout)::Order
		logical,optional,intent(in)::row_
		character*1,intent(in)::non_fermi_fuse
		logical::row
		integer::i,rank,rule
		type(QuanNum)::NewQuanNum
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		if(row)then
			rule=1
			outDim=B%SymDimension.subdim.[1,2]
			NewQuanNum=fuseOrder(order,B%quantumnumber(1),B%quantumnumber(2),rule)
			FermiFusefTensor=B%QuanFuse(NewQuanNum,Order,.true.)
		else
			rule=-1
			rank=B%getRank()
			i=rank-1
			outDim=B%SymDimension.subdim.[i,rank]
			NewQuanNum=fuseOrder(order,B%quantumnumber(i),B%quantumnumber(rank),rule)
			FermiFusefTensor=B%QuanFuse(NewQuanNum,Order,.false.)
		end if
		return
	end function
	
	
	type(fTensor) function FermiFusefTensor2(B,row_)result(FermiFusefTensor)
		class(fTensor),intent(in)::B
		type(Tensor)::Order
		logical,optional,intent(in)::row_
		logical::row
		integer::i,rank,rule
		type(QuanNum)::NewQuanNum
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		if(row)then
			rule=1
			NewQuanNum=fuseOrder(order,B%quantumnumber(1),B%quantumnumber(2),rule)
			FermiFusefTensor=B%QuanFuse(NewQuanNum,Order,.true.)
		else
			rule=-1
			rank=B%getRank()
			i=rank-1
			NewQuanNum=fuseOrder(order,B%quantumnumber(i),B%quantumnumber(rank),rule)
			FermiFusefTensor=B.pb.i
			call FermiFusefTensor%SymTensor%backward(i)
			FermiFusefTensor=FermiFusefTensor%QuanFuse(NewQuanNum,Order,.false.)
		end if
		return
	end function
	

	type(fTensor) function FermiFusefTensor2_(B,non_fermi_fuse,row_)result(FermiFusefTensor)
		class(fTensor),intent(in)::B
		character*1,intent(in)::non_fermi_fuse
		type(Tensor)::Order
		logical,optional,intent(in)::row_
		logical::row
		integer::i,rank,rule
		type(QuanNum)::NewQuanNum
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		if(row)then
			rule=1
			NewQuanNum=fuseOrder(order,B%quantumnumber(1),B%quantumnumber(2),rule)
			FermiFusefTensor=B%QuanFuse(NewQuanNum,Order,.true.)
		else
			rule=-1
			rank=B%getRank()
			i=rank-1
			NewQuanNum=fuseOrder(order,B%quantumnumber(i),B%quantumnumber(rank),rule)
			FermiFusefTensor=B%QuanFuse(NewQuanNum,Order,.false.)
		end if
		return
	end function


!The order of fusing all legs	 will store in outOrder
!orderinfo store the size(row number) of each fusing order
!example, fusing 4 leg, totaldata of orderinfo will be 3
!  fusing leg1 and leg2,result leg_1
!  order is
!     1 2 1 1 1 1 1
!     1 2 2 1 1 1 1
! outOrder is 
!     1 2 1 1 1 1 1
!     1 2 2 1 1 1 1
! outOrder(3)=2---> there are two row of data
!
!  fusing leg_1 and leg3,result leg_2
!  order is
!     2 2 1 1 1 1 1
! outOrder is 
!     2 2 1 1 1 1 1 ---> this is the neworder
!     1 2 1 1 1 1 1
!     1 2 2 1 1 1 1
! outOrder(2)=1---> there is 1 row of data
!
!  fusing leg_2 and leg4,result leg_3
!  order is
!     2 1 1 1 1 1 1
!     1 2 1 1 1 1 1
! outOrder is 
!     2 1 1 1 1 1 1 ---> this is the neworder
!     1 2 1 1 1 1 1 ---> this is the neworder
!     2 2 1 1 1 1 1 
!     1 2 1 1 1 1 1
!     1 2 2 1 1 1 1
! outOrder(1)=2---> there are 2 row of data
!
!When split data, read the outOrder, get the subTensor of order to split
	type(fTensor) function FermiFusefTensor3(B,ith,outDimen,outOrder,orderinfo,row_)result(T)
		class(fTensor),intent(in)::B
		integer,intent(in)::ith
		type(SymDimension),intent(inout)::outDimen
		type(Tensor),intent(inout)::outOrder,orderinfo
		logical,intent(in),optional::row_
		logical::row
		integer::i,rank,j,rankV
		type(Tensor)::order
		type(SymDimension)::dimen
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		call orderinfo%empty()
		if(row)then	
			if(ith.gt.1)then
				call orderinfo%allocate((/ith-1/),'integer')
				T=FermiFusefTensor(B,dimen,order,.true.)
				outDimen=dimen
				outOrder=order
				if(order%getRank().eq.0)then
					call writemess('ERROR in FermiFusefTensor3',-1)
					call error_stop
				end if
				if(order%getRank().eq.1)then
					call orderinfo%setValue(ith-1,1)
				else
					call orderinfo%setValue(ith-1,order%dim(1))
				end if
				do i=2,ith-1
					T=FermiFusefTensor(T,dimen,order,.true.)
					outDimen=outDimen+dimen
					if(order%getRank().eq.1)then
						call orderinfo%setValue(ith-i,1)
					else
						call orderinfo%setValue(ith-i,order%dim(1))
					end if
					outOrder=order.Rpaste.outOrder
					!call outOrder%paste(order,.true.)
				end do
			else
				T=B
			end if
		else
			rank=B%getRank()
			rankV=rank-ith+1
			if(rankV.gt.1)then
				call orderinfo%allocate((/rankV-1/),'integer')
				T=FermiFusefTensor(B,dimen,order,.false.)
				outDimen=dimen
				outOrder=order
				if(order%getRank().eq.0)then
					call writemess('ERROR in FermiFusefTensor3',-1)
					call error_stop
				end if
				if(order%getRank().eq.1)then
					call orderinfo%setValue(rankV-1,1)
				else
					call orderinfo%setValue(rankV-1,order%dim(1))
				end if
				do i=2,rankV-1
					T=FermiFusefTensor(T,dimen,order,.false.)
					outDimen=outDimen+dimen
					if(order%getRank().eq.1)then
						call orderinfo%setValue(rankV-i,1)
					else
						call orderinfo%setValue(rankV-i,order%dim(1))
					end if
					outOrder=order.Rpaste.outOrder
					!call outOrder%paste(order,.true.)
				end do
			else
				T=B
			end if
		end if
		
		return
	end function
	
	type(fTensor) function FermiFusefTensor4(B,ith,row_)result(T)
		class(fTensor),intent(in)::B
		integer,intent(in)::ith
		logical,intent(in),optional::row_
		logical::row
		integer::i,rank,j,rankV
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		if(row)then	
			if(ith.gt.1)then
				T=FermiFusefTensor2(B,.true.)
				do i=2,ith-1
					T=FermiFusefTensor2(T,.true.)
				end do
			else
				T=B
			end if
		else
			rank=B%getRank()
			rankV=rank-ith+1
			if(rankV.gt.1)then
				T=FermiFusefTensor2(B,.false.)
				do i=2,rankV-1
					T=FermiFusefTensor2(T,.false.)
				end do
			else
				T=B
			end if
		end if
		
		return
	end function
	
	type(fTensor) function FermiFusefTensor4_(B,ith,non_fermi_fuse,row_)result(T)
		class(fTensor),intent(in)::B
		integer,intent(in)::ith
		logical,intent(in),optional::row_
		character*1,intent(in)::non_fermi_fuse
		logical::row
		integer::i,rank,j,rankV
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		if(row)then	
			if(ith.gt.1)then
				T=FermiFusefTensor2_(B,non_fermi_fuse,.true.)
				do i=2,ith-1
					T=FermiFusefTensor2_(T,non_fermi_fuse,.true.)
				end do
			else
				T=B
			end if
		else
			rank=B%getRank()
			rankV=rank-ith+1
			if(rankV.gt.1)then
				T=FermiFusefTensor2_(B,non_fermi_fuse,.false.)
				do i=2,rankV-1
					T=FermiFusefTensor2_(T,non_fermi_fuse,.false.)
				end do
			else
				T=B
			end if
		end if
		
		return
	end function
	
	
	type(fTensor) function FermiSplitfTensor(B,NewQuanDimen,Order,row_)
		class(fTensor),intent(in)::B
		type(SymDimension),intent(in)::NewQuanDimen
		type(QuanNum)::NewQ1,NewQ2
		type(Tensor),intent(in)::Order
		logical,optional,intent(in)::row_
		logical::row
		integer::i
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		NewQ1=NewQuanDimen%QuantumNumber(1)
		NewQ2=NewQuanDimen%QuantumNumber(2)
		if(row)then
			FermiSplitfTensor=B%QuanSplit(NewQ1,NewQ2,Order,.true.)
		else
			i=B%getrank()
			FermiSplitfTensor=B%QuanSplit(NewQ1,NewQ2,Order,.false.)
			call FermiSplitfTensor%SymTensor%backward(i)
			call FermiSplitfTensor%backward(i)
		end if
		if(NewQuanDimen%outNameFlag().eq.1)then
			if(row)then
				call FermiSplitfTensor%setName(1,NewQuanDimen%outName(1))
				call FermiSplitfTensor%setName(2,NewQuanDimen%outName(2))
			else
				call FermiSplitfTensor%setName(FermiSplitfTensor%getRank()-1,NewQuanDimen%outName(1))
				call FermiSplitfTensor%setName(FermiSplitfTensor%getRank(),NewQuanDimen%outName(2))
			end if
		end if
		return
	end function
	type(fTensor) function FermiSplitfTensor2(B,inDimen,inOrder,orderinfo,row)result(T)
		class(fTensor),intent(in)::B
		type(SymDimension),intent(in)::inDimen
		type(Tensor),intent(in)::inOrder,orderinfo
		logical,intent(in),optional::row
		integer::i,rank,j,rankV,orderindex,istart,iend
		type(Tensor)::order
		type(SymDimension)::dimen
		rank=inDimen%getRank()
		if(rank.eq.0)then
			T=B
			return
		end if
		if(rank.lt.2)then
			call writemess('ERROR in Split fTensor, input parameter error',-1)
			call error_stop()
		end if
		orderindex=1
		dimen=inDimen.subdim.[rank-1,rank]
		istart=1
		iend=orderinfo%ii(orderindex)
		order=inOrder%subTensor((/-1,istart,iend/))
		istart=iend+1
		T=FermiSplitfTensor(B,dimen,Order,row)
		do i=rank-2,2,-2
			j=i-1
			orderindex=orderindex+1
			dimen=inDimen.subdim.[i-1,i]
			iend=istart+orderinfo%ii(orderindex)-1
			order=inOrder%subTensor((/-1,istart,iend/))
			istart=iend+1
			T=FermiSplitfTensor(T,dimen,order,row)
		end do
		return
	end function
	
	
	function SVD(T,QuanNumCut)
		type(fTensor),allocatable::SVD(:)
		class(fTensor),intent(in)::T
		type(QuanNum),optional,intent(in)::QuanNumCut
		allocate(SVD(3))
		call T%SVDSymTensorRoutine(SVD(1)%SymTensor,SVD(2)%SymTensor,SVD(3)%SymTensor,QuanNumCut)
		call SVD(1)%setFermiArrow(2,-1)
		call SVD(2)%setFermiArrow(1,1)
		call SVD(2)%setFermiArrow(2,-1)
		call SVD(3)%setFermiArrow(1,1)

		call SVD(1)%setRule(2,T%getRule(2))
		call SVD(2)%setRule(1,T%getRule(1))
		call SVD(2)%setRule(2,T%getRule(2))
		call SVD(3)%setRule(1,T%getRule(1))
		return
	end function

	function SVDNum(T,NumSave)result(SVD)
		type(fTensor),allocatable::SVD(:)
		class(fTensor),intent(in)::T
		integer,intent(in)::NumSave
		allocate(SVD(3))
		call T%SVDSymTensorRoutine(SVD(1)%SymTensor,SVD(2)%SymTensor,SVD(3)%SymTensor,NumSave)
		call SVD(1)%setFermiArrow(2,-1)
		call SVD(2)%setFermiArrow(1,1)
		call SVD(2)%setFermiArrow(2,-1)
		call SVD(3)%setFermiArrow(1,1)

		call SVD(1)%setRule(2,T%getRule(2))
		call SVD(2)%setRule(1,T%getRule(1))
		call SVD(2)%setRule(2,T%getRule(2))
		call SVD(3)%setRule(1,T%getRule(1))
		return
	end function
	
	subroutine SVDRoutineu1f(T,U,S,V,QuanNumCut)
		class(fTensor),intent(in)::T
		type(fTensor),intent(inout)::U,S,V
		type(QuanNum),optional,intent(in)::QuanNumCut
		call T%SVDSymTensorRoutine(U%SymTensor,S%SymTensor,V%SymTensor,QuanNumCut)
		call U%setFermiArrow(2,-1)
		call S%setFermiArrow(1,1)
		call S%setFermiArrow(2,-1)
		call V%setFermiArrow(1,1)

		call U%setRule(2,T%getRule(2))
		call S%setRule(1,T%getRule(1))
		call S%setRule(2,T%getRule(2))
		call V%setRule(1,T%getRule(1))
		return
	end subroutine
	subroutine SVDRoutineNumu1f(T,U,S,V,NumSave)
		class(fTensor),intent(in)::T
		type(fTensor),intent(inout)::U,S,V
		integer,intent(in)::NumSave
		call T%SVDSymTensorRoutine(U%SymTensor,S%SymTensor,V%SymTensor,NumSave)
		call U%setFermiArrow(2,-1)
		call S%setFermiArrow(1,1)
		call S%setFermiArrow(2,-1)
		call V%setFermiArrow(1,1)

		call U%setRule(2,T%getRule(2))
		call S%setRule(1,T%getRule(1))
		call S%setRule(2,T%getRule(2))
		call V%setRule(1,T%getRule(1))
		return
	end subroutine
	function SVDName(Tin,nameU,nameV,QuanNumCut)
		type(fTensor),allocatable::SVDName(:)
		class(fTensor),intent(in)::Tin
		type(QuanNum),optional,intent(in)::QuanNumCut
		character(len=*)::nameU,nameV
		type(fTensor)::T
		integer::rankU,rankV,i,j,rank
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		logical::reverFlag
		rank=Tin%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if(Tin%outTensorName(i).equ.nameU) rankU=rankU+1
			if(Tin%outTensorName(i).equ.nameV) rankV=rankV+1
		end do
		if(rankU.le.rankV)then
			T=Tin.pf.nameU
		else
			T=Tin.pb.nameV
		end if
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in SVDcutoff_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameU,-1)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameV,-1)
			call error_stop()
		end if
		T=T%FermiFuse(rankU,dimen(1),order(1),orderinfo(1),.true.)
		T=T%FermiFuse(2,dimen(2),order(2),orderinfo(2),.false.)
		if(T%getRule(1)*T%getRule(2).ge.0)then
			call reverseSymmetryRule(T%SymTensor)
			reverFlag=.true.
		else
			reverFlag=.false.
		end if
		allocate(SVDName(3))
		SVDName=SVD(T,QuanNumCut)
		if(reverFlag)then
			call reverseSymmetryRule(SVDName(3)%SymTensor)
		end if
		if(Tin%getNameFlag().ne.0)then
			call SVDName(1)%setName(2,'SVD.U')
			call SVDName(2)%setName(1,'SVD.s1')
			call SVDName(2)%setName(2,'SVD.s2')
			call SVDName(3)%setName(1,'SVD.V')
		end if
		SVDName(1)=SVDName(1)%FermiSplit(dimen(1),order(1),orderinfo(1),.true.)
		SVDName(3)=SVDName(3)%FermiSplit(dimen(2),order(2),orderinfo(2),.false.)
		return
	end function 

	function SVDNameNum(Tin,nameU,nameV,NumSave)Result(SVDName)
		type(fTensor),allocatable::SVDName(:)
		class(fTensor),intent(in)::Tin
		integer,intent(in)::NumSave
		character(len=*)::nameU,nameV
		type(fTensor)::T
		integer::rankU,rankV,i,j,rank
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		logical::reverFlag
		rank=Tin%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if(Tin%outTensorName(i).equ.nameU) rankU=rankU+1
			if(Tin%outTensorName(i).equ.nameV) rankV=rankV+1
		end do
		if(rankU.le.rankV)then
			T=Tin.pf.nameU
		else
			T=Tin.pb.nameV
		end if
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in SVDcutoff_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameU,-1)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameV,-1)
			call error_stop()
		end if
		T=T%FermiFuse(rankU,dimen(1),order(1),orderinfo(1),.true.)
		T=T%FermiFuse(2,dimen(2),order(2),orderinfo(2),.false.)
		if(T%getRule(1)*T%getRule(2).ge.0)then
			call reverseSymmetryRule(T%SymTensor)
			reverFlag=.true.
		else
			reverFlag=.false.
		end if
		allocate(SVDName(3))
		call T%SVDSymTensorRoutine(SVDName(1)%SymTensor,SVDName(2)%SymTensor,SVDName(3)%SymTensor,NumSave)
		call SVDName(1)%setFermiArrow(2,-1)
		call SVDName(2)%setFermiArrow(1,1)
		call SVDName(2)%setFermiArrow(2,-1)
		call SVDName(3)%setFermiArrow(1,1)

		call SVDName(1)%setRule(2,T%getRule(2))
		call SVDName(2)%setRule(1,T%getRule(1))
		call SVDName(2)%setRule(2,T%getRule(2))
		call SVDName(3)%setRule(1,T%getRule(1))

		if(Tin%getNameFlag().ne.0)then
			call SVDName(1)%setName(2,'SVD.U')
			call SVDName(2)%setName(1,'SVD.s1')
			call SVDName(2)%setName(2,'SVD.s2')
			call SVDName(3)%setName(1,'SVD.V')
		end if
		if(reverFlag)then
			call reverseSymmetryRule(SVDName(3)%SymTensor)
		end if
		SVDName(1)=SVDName(1)%FermiSplit(dimen(1),order(1),orderinfo(1),.true.)
		SVDName(3)=SVDName(3)%FermiSplit(dimen(2),order(2),orderinfo(2),.false.)
		return
	end function 

	function SVDNameNumMaxValue(Tin,nameU,nameV,minNum,NumSave,maxValue,Vtype)Result(SVDName)
		type(fTensor),allocatable::SVDName(:)
		class(fTensor),intent(in)::Tin
		integer,intent(in)::minNum
		integer,intent(inout)::NumSave
		real*8,intent(in)::maxValue
		character(len=*),optional,intent(in)::Vtype
		character(len=*)::nameU,nameV
		type(fTensor)::T
		integer::rankU,rankV,i,j,rank
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		logical::reverFlag
		rank=Tin%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if(Tin%outTensorName(i).equ.nameU) rankU=rankU+1
			if(Tin%outTensorName(i).equ.nameV) rankV=rankV+1
		end do
		if(rankU.le.rankV)then
			T=Tin.pf.nameU
		else
			T=Tin.pb.nameV
		end if
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in SVDcutoff_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameU,-1)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameV,-1)
			call error_stop()
		end if
		T=T%FermiFuse(rankU,dimen(1),order(1),orderinfo(1),.true.)
		T=T%FermiFuse(2,dimen(2),order(2),orderinfo(2),.false.)
		if(T%getRule(1)*T%getRule(2).ge.0)then
			call reverseSymmetryRule(T%SymTensor)
			reverFlag=.true.
		else
			reverFlag=.false.
		end if
		allocate(SVDName(3))
		call T%SVDSymTensorRoutine(SVDName(1)%SymTensor,SVDName(2)%SymTensor,SVDName(3)%SymTensor,&
			minNum,NumSave,maxValue,Vtype)
		call SVDName(1)%setFermiArrow(2,-1)
		call SVDName(2)%setFermiArrow(1,1)
		call SVDName(2)%setFermiArrow(2,-1)
		call SVDName(3)%setFermiArrow(1,1)

		call SVDName(1)%setRule(2,T%getRule(2))
		call SVDName(2)%setRule(1,T%getRule(1))
		call SVDName(2)%setRule(2,T%getRule(2))
		call SVDName(3)%setRule(1,T%getRule(1))

		if(Tin%getNameFlag().ne.0)then
			call SVDName(1)%setName(2,'SVD.U')
			call SVDName(2)%setName(1,'SVD.s1')
			call SVDName(2)%setName(2,'SVD.s2')
			call SVDName(3)%setName(1,'SVD.V')
		end if
		if(reverFlag)then
			call reverseSymmetryRule(SVDName(3)%SymTensor)
		end if
		SVDName(1)=SVDName(1)%FermiSplit(dimen(1),order(1),orderinfo(1),.true.)
		SVDName(3)=SVDName(3)%FermiSplit(dimen(2),order(2),orderinfo(2),.false.)
		return
	end function 

	subroutine SVDNameNumRoutine(Tin,U,S,V,SVDname,NumSave,leftFlag)
		class(fTensor),intent(in)::Tin
		type(fTensor),intent(inout)::U,S,V
		integer,intent(in)::NumSave
		character(len=*)::SVDname(:)
		logical,intent(in)::leftFlag
		type(fTensor)::T
		integer::lenSVDname,rank,lenLeft
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		logical::reverFlag
		rank=Tin%getRank()
		lenSVDname=size(SVDname)
		if(lenSVDname.ge.rank)then
			call writemess("ERROR in SVDNameNumRoutine, input name larger than the rank of the Tensor",-1)
			call error_stop()
		end if
		if(leftFlag)then
			T=Tin.pf.SVDname
			lenLeft=lenSVDname
		else
			T=Tin.pb.SVDname
			lenLeft=rank-lenSVDname
		end if
		T=T%FermiFuse(lenLeft,dimen(1),order(1),orderinfo(1),.true.)
		T=T%FermiFuse(2,dimen(2),order(2),orderinfo(2),.false.)

		if(T%getRule(1)*T%getRule(2).ge.0)then
			call reverseSymmetryRule(T%SymTensor)
			reverFlag=.true.
		else
			reverFlag=.false.
		end if
		call T%SVDSymTensorRoutine(U%SymTensor,S%SymTensor,V%SymTensor,NumSave)
		call U%setFermiArrow(2,-1)
		call S%setFermiArrow(1,1)
		call S%setFermiArrow(2,-1)
		call V%setFermiArrow(1,1)

		call U%setRule(2,T%getRule(2))
		call S%setRule(1,T%getRule(1))
		call S%setRule(2,T%getRule(2))
		call V%setRule(1,T%getRule(1))

		if(Tin%getNameFlag().ne.0)then
			call U%setName(2,'SVD.U')
			call S%setName(1,'SVD.s1')
			call S%setName(2,'SVD.s2')
			call V%setName(1,'SVD.V')
		end if
		if(reverFlag)then
			call reverseSymmetryRule(V%SymTensor)
		end if
		U=U%FermiSplit(dimen(1),order(1),orderinfo(1),.true.)
		V=V%FermiSplit(dimen(2),order(2),orderinfo(2),.false.)
		return
	end subroutine 

	subroutine SVDNameNumMaxValueRoutine(Tin,U,S,V,maxValue,Vtype,SVDname,minNum,NumSave,leftFlag)
		class(fTensor),intent(in)::Tin
		type(fTensor),intent(inout)::U,S,V
		real*8,intent(in)::maxValue
		character(len=*),intent(in)::Vtype
		integer,intent(in)::minNum
		integer,intent(inout)::NumSave
		character(len=*)::SVDname(:)
		logical,intent(in)::leftFlag
		type(fTensor)::T
		integer::lenSVDname,rank,lenLeft
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		logical::reverFlag
		rank=Tin%getRank()
		lenSVDname=size(SVDname)
		if(lenSVDname.ge.rank)then
			call writemess("ERROR in SVDNameNumRoutine, input name larger than the rank of the Tensor",-1)
			call error_stop()
		end if
		if(leftFlag)then
			T=Tin.pf.SVDname
			lenLeft=lenSVDname
		else
			T=Tin.pb.SVDname
			lenLeft=rank-lenSVDname
		end if
		T=T%FermiFuse(lenLeft,dimen(1),order(1),orderinfo(1),.true.)
		T=T%FermiFuse(2,dimen(2),order(2),orderinfo(2),.false.)

		if(T%getRule(1)*T%getRule(2).ge.0)then
			call reverseSymmetryRule(T%SymTensor)
			reverFlag=.true.
		else
			reverFlag=.false.
		end if
		call T%SVDSymTensorRoutine(U%SymTensor,S%SymTensor,V%SymTensor,minNum,NumSave,maxValue,Vtype)
		call U%setFermiArrow(2,-1)
		call S%setFermiArrow(1,1)
		call S%setFermiArrow(2,-1)
		call V%setFermiArrow(1,1)

		call U%setRule(2,T%getRule(2))
		call S%setRule(1,T%getRule(1))
		call S%setRule(2,T%getRule(2))
		call V%setRule(1,T%getRule(1))

		if(Tin%getNameFlag().ne.0)then
			call U%setName(2,'SVD.U')
			call S%setName(1,'SVD.s1')
			call S%setName(2,'SVD.s2')
			call V%setName(1,'SVD.V')
		end if
		if(reverFlag)then
			call reverseSymmetryRule(V%SymTensor)
		end if
		U=U%FermiSplit(dimen(1),order(1),orderinfo(1),.true.)
		V=V%FermiSplit(dimen(2),order(2),orderinfo(2),.false.)
		return
	end subroutine 

	
	
	subroutine QRdecompositionfTensor(T,Q,R) 
		class(fTensor),intent(in)::T
		type(fTensor),intent(inout)::Q,R
		call T%SymTensor%SymQRRoutine(Q%SymTensor,R%SymTensor)
		call Q%setFermiArrow(2,-1)
		call R%setFermiArrow(1,1)
		call Q%setRule(2,T%getRule(2))
		call R%setRule(1,T%getRule(1))
		return
	end subroutine
	function QRfTensor_noname(T)result(Res)
		type(fTensor),allocatable::Res(:)
		class(fTensor),intent(in)::T
		allocate(Res(2))
		call QRdecompositionfTensor(T,Res(1),Res(2)) 
		return
	end function
	
	function QRfTensor_name(Tin,nameU,nameV)result(Res)
		type(fTensor),allocatable::Res(:)
		character(len=*),intent(in)::nameU,nameV
		class(fTensor),intent(in)::Tin
		integer::rankU,rankV,i,j,rank
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		Type(fTensor)::T
		logical::reverFlag
		rank=Tin%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if(Tin%outTensorName(i).equ.nameU) rankU=rankU+1
			if(Tin%outTensorName(i).equ.nameV) rankV=rankV+1
		end do
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in QRfTensor_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call Tin%diminfo()
			call writemess(nameU+','+nameV,-1)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			T=Tin.pf.nameU
		else
			T=Tin.pb.nameV
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in QRfTensor_name,no such name",-1)
			call writemess(nameU,-1)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in QRfTensor_name,no such name",-1)
			call writemess(nameV,-1)
			call error_stop()
		end if
		T=T%FermiFuse(rankU,dimen(1),order(1),orderinfo(1),.true.)
		T=T%FermiFuse(2,dimen(2),order(2),orderinfo(2),.false.)
		if(T%getRule(1)*T%getRule(2).ge.0)then
			call reverseSymmetryRule(T%SymTensor)
			reverFlag=.true.
		else
			reverFlag=.false.
		end if
		allocate(Res(2))
		call QRdecompositionfTensor(T,Res(1),Res(2)) 
		if(Tin%getNameFlag().ne.0)then
			call Res(1)%setName(2,'QR.Q')
			call Res(2)%setName(1,'QR.R')
		end if
		if(reverFlag)then
			call reverseSymmetryRule(Res(2)%SymTensor)
		end if
		Res(1)=Res(1)%FermiSplit(dimen(1),order(1),orderinfo(1),.true.)
		Res(2)=Res(2)%FermiSplit(dimen(2),order(2),orderinfo(2),.false.)
	end function

	function QRfTensor_name2(Tin,QRName,LeftFlag_)result(Res)
		type(fTensor),allocatable::Res(:)
		character(len=*),intent(in)::QRName(:)
		logical,optional,intent(in)::LeftFlag_
		class(fTensor),intent(in)::Tin
		logical::LeftFlag,reverFlag
		integer::lenQRName,i,j,rank
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		Type(fTensor)::T
		integer::lenLeft
		if(present(LeftFlag_))then
			LeftFlag=LeftFlag_
		else
			LeftFlag=.true.
		end if
		rank=Tin%getRank()
		lenQRName=size(QRName)
		if(lenQRName.ge.rank)then
			call writemess("ERROR in QRfTensor, input name larger than the rank of the Tensor",-1)
			call error_stop()
		end if
		if(LeftFlag)then
			T=Tin.pf.QRName
			lenLeft=lenQRName
		else
			T=Tin.pb.QRName
			lenLeft=rank-lenQRName
		end if
		T=T%FermiFuse(lenLeft,dimen(1),order(1),orderinfo(1),.true.)
		T=T%FermiFuse(2,dimen(2),order(2),orderinfo(2),.false.)

		if(T%getRule(1)*T%getRule(2).ge.0)then
			call reverseSymmetryRule(T%SymTensor)
			reverFlag=.true.
		else
			reverFlag=.false.
		end if
		allocate(Res(2))
		call QRdecompositionfTensor(T,Res(1),Res(2)) 
		if(Tin%getNameFlag().ne.0)then
			call Res(1)%setName(2,'QR.Q')
			call Res(2)%setName(1,'QR.R')
		end if
		if(reverFlag)then
			call reverseSymmetryRule(Res(2)%SymTensor)
		end if
		Res(1)=Res(1)%FermiSplit(dimen(1),order(1),orderinfo(1),.true.)
		Res(2)=Res(2)%FermiSplit(dimen(2),order(2),orderinfo(2),.false.)
	end function

	subroutine QRfTensorRoutine_name2(Tin,Q,R,QRName,LeftFlag_)
		character(len=*),intent(in)::QRName(:)
		logical,optional,intent(in)::LeftFlag_
		class(fTensor),intent(in)::Tin
		type(fTensor),intent(inout)::Q,R
		logical::LeftFlag,reverFlag
		integer::lenQRName,i,j,rank
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		Type(fTensor)::T
		integer::lenLeft
		if(present(LeftFlag_))then
			LeftFlag=LeftFlag_
		else
			LeftFlag=.true.
		end if
		rank=Tin%getRank()
		lenQRName=size(QRName)
		if(lenQRName.ge.rank)then
			call writemess("ERROR in QRfTensor, input name larger than the rank of the Tensor",-1)
			call error_stop()
		end if
		if(LeftFlag)then
			T=Tin.pf.QRName
			lenLeft=lenQRName
		else
			T=Tin.pb.QRName
			lenLeft=rank-lenQRName
		end if
		T=T%FermiFuse(lenLeft,dimen(1),order(1),orderinfo(1),.true.)
		T=T%FermiFuse(2,dimen(2),order(2),orderinfo(2),.false.)

		if(T%getRule(1)*T%getRule(2).ge.0)then
			call reverseSymmetryRule(T%SymTensor)
			reverFlag=.true.
		else
			reverFlag=.false.
		end if
		call QRdecompositionfTensor(T,Q,R) 
		if(Tin%getNameFlag().ne.0)then
			call Q%setName(2,'QR.Q')
			call R%setName(1,'QR.R')
		end if
		if(reverFlag)then
			call reverseSymmetryRule(R%SymTensor)
		end if
		Q=Q%FermiSplit(dimen(1),order(1),orderinfo(1),.true.)
		R=R%FermiSplit(dimen(2),order(2),orderinfo(2),.false.)
	end subroutine
	
	subroutine LQdecompositionfTensor(T,L,Q) 
		class(fTensor),intent(in)::T
		type(fTensor),intent(inout)::L,Q
		call T%SymTensor%SymLQRoutine(L%SymTensor,Q%SymTensor)
		call L%setFermiArrow(2,-1)
		call Q%setFermiArrow(1,1)
		call L%setRule(2,T%getRule(2))
		call Q%setRule(1,T%getRule(1))
		return
	end subroutine
	function LQfTensor_noname(T)result(Res)
		type(fTensor),allocatable::Res(:)
		class(fTensor),intent(in)::T
		allocate(Res(2))
		call LQdecompositionfTensor(T,Res(1),Res(2)) 
		return
	end function
	
	function LQfTensor_name(Tin,nameU,nameV)result(Res)
		type(fTensor),allocatable::Res(:)
		character(len=*),intent(in)::nameU,nameV
		class(fTensor),intent(in)::Tin
		integer::rankU,rankV,i,j,rank
		Type(fTensor)::T
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		logical::reverFlag
		rank=Tin%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if(Tin%outTensorName(i).equ.nameU) rankU=rankU+1
			if(Tin%outTensorName(i).equ.nameV) rankV=rankV+1
		end do
		if(rankU.le.rankV)then
			T=Tin.pf.nameU
		else
			T=Tin.pb.nameV
		end if
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in LQfTensor_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in LQfTensor_name,no such name",-1)
			call writemess(nameU,-1)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in LQfTensor_name,no such name",-1)
			call writemess(nameV,-1)
			call error_stop()
		end if
		T=T%FermiFuse(rankU,dimen(1),order(1),orderinfo(1),.true.)
		T=T%FermiFuse(2,dimen(2),order(2),orderinfo(2),.false.)
		if(T%getRule(1)*T%getRule(2).ge.0)then
			call reverseSymmetryRule(T%SymTensor)
			reverFlag=.true.
		else
			reverFlag=.false.
		end if
		allocate(Res(2))
		call LQdecompositionfTensor(T,Res(1),Res(2)) 
		if(Tin%getNameFlag().ne.0)then
			call Res(1)%setName(2,'LQ.L')
			call Res(2)%setName(1,'LQ.Q')
		end if
		if(reverFlag)then
			call reverseSymmetryRule(Res(2)%SymTensor)
		end if
		Res(1)=Res(1)%FermiSplit(dimen(1),order(1),orderinfo(1),.true.)
		Res(2)=Res(2)%FermiSplit(dimen(2),order(2),orderinfo(2),.false.)
		return
	end function
	
	function LQfTensor_name2(Tin,LQName,LeftFlag_)result(Res)
		type(fTensor),allocatable::Res(:)
		character(len=*),intent(in)::LQName(:)
		logical,optional,intent(in)::LeftFlag_
		class(fTensor),intent(in)::Tin
		logical::LeftFlag,reverFlag
		integer::lenLQName,i,j,rank
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		Type(fTensor)::T
		integer::lenLeft
		if(present(LeftFlag_))then
			LeftFlag=LeftFlag_
		else
			LeftFlag=.true.
		end if
		rank=Tin%getRank()
		lenLQName=size(LQName)
		if(lenLQName.gt.rank)then
			call writemess("ERROR in QRfTensor, input name larger than the rank of the Tensor",-1)
			call error_stop()
		end if
		if(LeftFlag)then
			T=Tin.pf.LQName
			lenLeft=lenLQName
		else
			T=Tin.pb.LQName
			lenLeft=rank-lenLQName
		end if

		T=T%FermiFuse(lenLeft,dimen(1),order(1),orderinfo(1),.true.)
		T=T%FermiFuse(2,dimen(2),order(2),orderinfo(2),.false.)


		if(T%getRule(1)*T%getRule(2).ge.0)then
			call reverseSymmetryRule(T%SymTensor)
			reverFlag=.true.
		else
			reverFlag=.false.
		end if
		allocate(Res(2))
		call LQdecompositionfTensor(T,Res(1),Res(2)) 
		if(Tin%getNameFlag().ne.0)then
			call Res(1)%setName(2,'LQ.L')
			call Res(2)%setName(1,'LQ.Q')
		end if
		if(reverFlag)then
			call reverseSymmetryRule(Res(2)%SymTensor)
		end if
		Res(1)=Res(1)%FermiSplit(dimen(1),order(1),orderinfo(1),.true.)
		Res(2)=Res(2)%FermiSplit(dimen(2),order(2),orderinfo(2),.false.)
		return
	end function
		
	
	type(fTensor) function directProductTensor(T1,T2)result(Res)
		type(fTensor),intent(in) :: T1,T2
		Res%SymTensor=T1%SymTensor.kron.T2%SymTensor
		return
	end function
	
	
	type(Tensor) function TdotTensor(phi1,phi2)result(dotTensor)
		Type(fTensor),intent(in)::phi1,phi2
		dotTensor=phi1%SymTensor.x.phi2%SymTensor
		RETURN
	end function
	
	


! reorder the Tensor, order the leg as the order in allTensorName(:)

	subroutine  DimOrder(inoutT,allTensorName)
		class(fTensor),intent(inout)::inoutT
		character(len=*),intent(in)::allTensorName(:)
		character(len=max_len_of_char_in_TData),allocatable :: TensorName(:)
		character(len=max_len_of_char_in_TData),allocatable :: Names(:)
		integer::rank,i
		if(.not.inoutT%getFlag())then
			call writemess('There is no data in the fTensor',-1)
			call error_stop
		endif
		rank=inoutT%getRank()
		allocate(TensorName(rank))
		allocate(Names(rank))
		do i=1,rank
			TensorName(i)=inoutT%outTensorName(i)
			Names(i)=inoutT%getName(i)
		end do
		call NameOrder_sort(TensorName,Names,allTensorName,rank)
		call inoutT%permute(Names)
		return
	end subroutine

!Do not Check if There is no such name

	subroutine NameOrder_sort(TensorName,inoutName,allTensorName,n)
		integer,intent(in)::n
		character(len=*),intent(in)::allTensorName(:)
		character(len=*),intent(inout) :: inoutName(n),TensorName(n)
		character(len=max_len_of_char_in_TData) :: temp
		integer :: i,j
		do i=1,n-1
		 do j=i+1,n
			  if ( checkOrder(TensorName(i),TensorName(j),allTensorName) ) then
					temp = inoutName(i)
					inoutName(i) = inoutName(j)
					inoutName(j) = temp
					temp = TensorName(i)
					TensorName(i) = TensorName(j)
					TensorName(j) = temp
			  endif
		 enddo
		enddo
		return
	end subroutine


	logical function checkOrder(name1,name2,allTensorName)
		character(len=*),intent(in)::name1,name2,allTensorName(:)
		integer::i,lenName
		lenName=size(allTensorName)
		do i=1,lenName
			if(name1.equ.allTensorName(i))then
				checkOrder=.false.
				return
			end if
			if(name2.equ.allTensorName(i))then
				checkOrder=.true.
				return
			end if
		end do
		call writemess('ERROR in checkOrder',-1)
		call error_stop
	end function



!**********************************************************************
!**********************************************************************
!	the code below is for MPI
!**********************************************************************

	subroutine MPI_send_fTensor(Ten1,Ten2,ID1,ID2,ierr,MPIcommon)
		type(fTensor),intent(in)::Ten1
		type(fTensor),intent(inout)::Ten2
		integer,intent(in)::ID1,ID2
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		call MPI_send_SymTensor(Ten1%SymTensor,Ten2%SymTensor,ID1,ID2,ierr,MPIcommon)
		return
	end subroutine
	subroutine MPI_BCAST_fTensor(Ten1,ID,ierr,MPIcommon)
		type(fTensor),intent(inout)::Ten1
		integer,intent(in)::ID
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		call MPI_BCAST_SymTensor(Ten1%SymTensor,ID,ierr,MPIcommon)
		return
	end subroutine



	subroutine MPI_SUM_fTensor1(inoutTensor,ierr,MPIcommon)
		type(fTensor),intent(inout)::inoutTensor
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		call MPI_SUM_SymTensor(inoutTensor%SymTensor,ierr,MPIcommon)
		return
	end subroutine



end module
