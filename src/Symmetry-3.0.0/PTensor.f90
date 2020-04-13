module ParityTool!The parity of is +1
	use SymDimension_typede
	use usefull_function
	use Tensor_type
	implicit none
	private
	public::NonZeroElement,rule,fuseOrder,quantumnumber
	interface rule
		module procedure Rule1
		module procedure Rule2
		module procedure PRule1
		module procedure PRule2
		module procedure PRule3
	end interface
	interface fuseOrder
		module procedure QuanFuseOrderFunction
	end interface
	interface quantumnumber
		module procedure Parity
		module procedure Parity2
		module procedure Parity3
		module procedure Parity4
	end interface
	
	
	
contains
	
	
	
	logical function Rule1(i,j,k,QN1,QN2,QN3)
		integer,intent(in)::i,j,k
		type(QuanNum),intent(in)::QN1,QN2,QN3
		real*4::Q1,Q2,Q3
		Q1=QN1%getQN(i)
		Q2=QN2%getQN(j)
		Q3=QN3%getQN(k)
		Rule1=PRule1(Q1,Q2,Q3)
	end function
	logical function Rule2(indices,QN)
		integer,intent(in)::indices(:)
		type(QuanNum),intent(in)::QN(:)
		real*4::Quan(size(QN))
		integer::i
		do i=1,size(QN)
			Quan(i)=QN(i)%getQN(indices(i))
		end do
		Rule2=PRule2(Quan)
	end function
	logical function PRule1(Q1,Q2,Q3)
		real*4,intent(in)::Q1,Q2,Q3
		if(Q1.le.0)then
			PRule1=(Q2*Q3).le.0
		else
			PRule1=(Q2*Q3).gt.0
		end if
		return
	end function
	logical function PRule2(Q)
		real*4,intent(in)::Q(:)
		real*4::temp
		integer::i
		temp=Q(1)
		do i=2,size(Q)
			temp=temp*Q(i)
		end do
		PRule2=temp.gt.0
		return
	end function
	
	logical function PRule3(dimen,indices)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::indices(:)
		real*4::Q(size(indices))
		integer::i
		real*4::temp
		do i=1,size(indices)
			Q(i)=dimen%getQN(i,indices(i))
		end do
		temp=Q(1)
		do i=2,size(Q)
			temp=temp*Q(i)
		end do
		PRule3=temp.gt.0
		return
	end function
	
	
	type(QuanNum) function Parity(deg)
		integer,intent(in)::deg(2)
		call Parity%setRule(0)
		call Parity%setMaxQN(1.)
		call Parity%setQN((/-1.,1./))
		call Parity%setDeg(Deg)
		return
	end function
	type(QuanNum) function Parity2(deg,rule)result(Parity)
		integer,intent(in)::deg(2),rule
		call Parity%setRule(rule)
		call Parity%setMaxQN(1.)
		call Parity%setQN((/-1.,1./))
		call Parity%setDeg(Deg)
		return
	end function
	type(QuanNum) function Parity3(QN,deg)result(Parity)
		real*4::QN
		integer,intent(in)::deg
		call Parity%setRule(0)
		call Parity%setMaxQN(QN)
		call Parity%setQN((/QN/))
		call Parity%setDeg(1,Deg)
		return
	end function
	type(QuanNum) function Parity4(QN,deg,rule)result(Parity)
		real*4::QN
		integer,intent(in)::deg,rule
		call Parity%setRule(rule)
		call Parity%setMaxQN(QN)
		call Parity%setQN((/QN/))
		call Parity%setDeg(1,Deg)
		return
	end function
	
!outQ=Q1+Q2
	type(QuanNum) function QuanFuseOrderFunction(order,Q1,Q2)
		type(QuanNum),intent(in)::Q1,Q2
		type(Tensor),intent(inout)::order
		call QuanFuseOrder(order,QuanFuseOrderFunction,Q1,Q2)
		return
	end function

!
!order of codes for do
! outQ=Q1+Q2 
!  -1       -1  -1   false
!  -1       +1  -1
!  -1       -1  +1
!  -1       +1  +1   false
!  +1       -1  -1
!  +1       +1  -1   false
!  +1       -1  +1   false
!  +1       +1  +1
!order output outi,Q1i,Q2i,degQ1,degQ2,degstart,DegEnd
!Use in fuse and split Tensor
!example       Q1              Q2
! QN          -   +         -    +
! Deg         2   3         4    5
! index       1   2         1    2
! Q3=Q1+Q2
!                    Q3
!             -                  +
!        (+ -)  (- +)         (- -)  (+ +)
!         3*4    2*5           2*4   3*5 
! index       1                   2
!output(1): -  --->  (+ -)
!      outi    =1
!      Q1i     =2
!      Q2i     =1
!      degQ1   =3
!      degQ2   =4
!      degstart=1
!      DegEnd  =1+3*4
!output(2): -  --->  (- +)
!      outi    =1
!      Q1i     =1
!      Q2i     =2
!      degQ1   =2
!      degQ2   =5
!      degstart=(1+3*4)+1
!      DegEnd  =(1+3*4+1)+2*5
!Do not store the data that degQ1=0 or degQ2=0
! 
	subroutine QuanFuseOrder(order,outQ,Q1,Q2)
		type(QuanNum),intent(inout)::outQ
		type(QuanNum),intent(in)::Q1,Q2
		type(Tensor),intent(inout)::order
		integer::i,j,k,deg(2),deg1,deg2,degstart,maxi,maxj,maxk
		real*4::outQN,QN1,QN2
		logical::goon,deg_not_zero_flag
		maxj=Q2%getQNlength()
		maxk=Q1%getQNlength()
		maxi=max(maxj,maxk)
		if(maxi.eq.1)then
			call QuanFuseOrder_one_QN(order,outQ,Q1,Q2)
			return
		end if
		call outQ%setRule(Q1%getRule())
		call outQ%setMaxQN(1.)
		call outQ%setQN((/-1.,1./))
		deg=0
		call order%empty()
		do i=1,maxi
			degstart=0
			do j=1,maxj
				do k=1,maxk
					outQN=outQ%getQN(i)
					QN1=Q1%getQN(k)
					QN2=Q2%getQN(j)
					goon=PRule1(outQN,QN1,QN2)
					if(goon)then
						deg1=Q1%getDeg(k)
						deg2=Q2%getDeg(j)
						deg_not_zero_flag=(deg1.ne.0).and.(deg2.ne.0)
						if(deg_not_zero_flag)then
							deg(i)=deg(i)+deg1*deg2
							if(order%getFlag())then
								call order%addrow(Tensor( (/i,k,j,deg1,deg2,degstart+1,deg(i)/)))
								degstart=deg(i)
							else
								order=(/i,k,j,deg1,deg2,1,deg(i)/)
								degstart=deg(i)
							end if
						end if
					end if
				end do
			end do
		end do
		call outQ%setDeg(Deg)
		if(order%getRank().eq.1)call order%resetDim((/1,order%getTotalData()/))
		return
	end subroutine
	subroutine QuanFuseOrder_one_QN(order,outQ,Q1,Q2)
		type(QuanNum),intent(inout)::outQ
		type(QuanNum),intent(in)::Q1,Q2
		type(Tensor),intent(inout)::order
		integer::i,j,k,deg,deg1,deg2
		real*4::outQN,QN1,QN2
		logical::deg_not_zero_flag
		QN1=Q1%getQN(1)
		QN2=Q2%getQN(1)
		if(QN1*QN2.gt.0)then
			outQN=1.
		else
			outQN=-1.
		end if
		call outQ%setRule(1)
		call outQ%setMaxQN(outQN)
		call outQ%setQN((/outQN/))
		deg=0
		deg1=Q1%getDeg(1)
		deg2=Q2%getDeg(1)
		deg_not_zero_flag=(deg1.ne.0).and.(deg2.ne.0)
		if(deg_not_zero_flag)then
			deg=deg1*deg2
			order=(/1,1,1,deg1,deg2,1,deg/)
		end if
		call outQ%setDeg(1,Deg)
		call order%resetDim((/1,order%getTotalData()/))
		return
	end subroutine
		

	type(Tensor) function NonZeroElement(dimen)
		class(SymDimension),intent(in)::dimen
		call NonZeroElementParity(NonZeroElement,dimen)
		if(NonZeroElement%getRank().eq.1)then
			call NonZeroElement%resetDim((/1,NonZeroElement%getTotalData()/))
		end if
		return
	end function
	subroutine NonZeroElementParity(Res,dimen)
		type(Tensor),intent(inout)::Res
		class(SymDimension),intent(in)::dimen
		integer,allocatable::indices(:),maxdim(:),mindim(:)
		integer::Rank
		logical::goon,rulefit
		real*4,allocatable::QN(:)
		rank=dimen%getRank()
		allocate(indices(rank))
		allocate(maxdim(rank))
		allocate(mindim(rank))
		allocate(QN(rank))
		mindim=1
		indices=mindim
		maxdim=dimen%dim()
		goon=.true.
		call Res%empty()
		do while(goon)
			QN=dimen%QNDim(indices)
			rulefit=Prule2(QN)
			if(rulefit)then
				if(Res%getFlag())then
					call Res%addrow(Tensor(indices))
				else
					Res=indices
				end if
			end if
			goon=inde_counter(indices,mindim,maxdim,1)
		end do
		return
	end subroutine
	
	
end module


module ParityTensor
	use SymTensor_type
	use SymDimension_typede
	use Tensor_type
	use Dimension_typede
	use usefull_function
	use ParityTool
	implicit none
	private
	
	public::PTensor
	type,extends (SymTensor) :: PTensor
	contains
		generic,public::allocatePBlock=>allocatePBlock1,allocatePBlock2,allocatePBlock3,allocatePBlockAll
		generic,public::randomPTensor=>set_random_PTensor_val,set_random_PTensor_vec,set_random_PTensor_all	
		generic,public::setPValue=>set_Symelement_vec,set_Symelement_val!there is symmetry check in the input
		generic,public::SVDPTensor=>SVD,SVDName
		procedure,public::SVDPRoutine=>SVDRoutine
		procedure,public::SymmetryCheck
		generic,public::ParitySplit=>ParitySplitPTensor2,ParitySplitPTensor
		generic,public::ParityFuse=>ParityFusePTensor,ParityFusePTensor2,ParityFusePTensor3,ParityFusePTensor4
		generic,public::QRPTensor=>QRPTensor_noName,QRPTensor_name!It use QuanFuse not ParityFuse
		procedure,public::QRPRoutine=>QRdecompositionPTensor
		generic,public::LQPTensor=>LQPTensor_noName,LQPTensor_name!It use QuanFuse not ParityFuse
		procedure,public::LQPRoutine=>LQdecompositionPTensor
		generic,public::eye=>eye_Tensor3,eye_Tensor4,eye_Tensor5
		procedure,public::reorder=>DimOrder
		! reorder the Tensor, order the leg as the order in allTensorName(:)
		
		procedure::allocateTensor2
		procedure::allocateTensor5
		procedure::allocateTensor8
		procedure::allocatePBlock1
		procedure::allocatePBlock2
		procedure::allocatePBlock3
		procedure::allocatePBlockAll
		procedure::set_random_PTensor_val
		procedure::set_random_PTensor_vec
		procedure::set_random_PTensor_all	
		procedure::set_Symelement_vec
		procedure::set_Symelement_val!As the dimension of the block is fix by the degeneracy of the symDimension
		                                           !This subroutine will ignore the input Tensor dimension
		procedure::SVD
		procedure::SVDName
		procedure::ParitySplitPTensor
		procedure::ParitySplitPTensor2
		procedure::ParityFusePTensor
		procedure::ParityFusePTensor2
		procedure::ParityFusePTensor3
		procedure::ParityFusePTensor4
		procedure::QRPTensor_noName
		procedure::QRPTensor_name
		procedure::LQPTensor_noName
		procedure::LQPTensor_name
		procedure::eye_Tensor3
		procedure::eye_Tensor4
		procedure::eye_Tensor5
		procedure::eye_Tensor6
	end type PTensor
	
	
	public::expm
	interface expm
		module procedure  expmPTensor
	end interface
	
	public::operator(.pf.)
	interface operator(.pf.)
		module procedure ParityPermuteForWard
		module procedure ParityPermuteForWard_vec
		module procedure ParityPermuteForWard_name
		module procedure ParityPermuteForWard_name_vec
		module procedure ParityPermuteForWard_Tensor
	end interface
	public::operator(.pb.)
	interface operator(.pb.)
		module procedure ParityPermuteBackWard
		module procedure ParityPermuteBackWard_vec
		module procedure ParityPermuteBackWard_name
		module procedure ParityPermuteBackWard_name_vec
		module procedure ParityPermuteBackWard_Tensor
	end interface
	public::operator(.p.)
	interface operator(.p.)
		module procedure ParityPermutation
		module procedure ParityPermutation_name
		module procedure ParityPermutation_Tensor
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
		module procedure add_PTensor
	end interface
	
	public::operator(-)
	interface operator(-)
		module procedure minu_PTensor
	end interface
	
	public::operator(.H.)
	interface operator(.H.)
		module procedure Htranspose! Htranspose, and the symmetry rule will go inverse
	end interface
	public::operator(.Hn.)
	interface operator(.Hn.)
		module procedure Htranspose2! Htranspose, and the symmetry rule will go inverse
	end interface
	
	public::operator(.subdim.)!overwrite the function in type(Dimension)
	interface operator(.subdim.)
		module procedure getSubDim2
		module procedure getSubDim3
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
		module procedure PTensorToSymTensor
		module procedure PTensorToTensor
		module procedure SymTensorToPTensor!Do not check the Symmetry rule
		module procedure TensorToPTensor!Do not check the Symmetry rule
	end interface
	
	public::MPI_BCAST_PTensor,MPI_send_PTensor
	public::MPI_SUM_PTensor!Do not check the input PTensor
	interface MPI_SUM_PTensor
		module procedure MPI_SUM_PTensor1
	end interface
contains
	logical function SymmetryCheck(T)
		class(PTensor),intent(in)::T
		integer::i
		integer,allocatable::Tdim(:),Maxdim(:)
		allocate(TDim(T%GetRank()))
		allocate(Maxdim(T%GetRank()))
		Maxdim=T%dim()
		SymmetryCheck=.true.
		do i=1,T%getTotalData()
			if(T%getFlag(i))then
				call IndesToaddress(Maxdim,Tdim,i)	
				SymmetryCheck=rule(T%SymDimension,TDim)
				if(.not.SymmetryCheck)return
			end if
		end do
	end function
	
	subroutine assignmentTen(T,T2)
		type(PTensor),intent(inout) ::T
		type(PTensor),intent(in) :: T2
		T%SymTensor=T2%SymTensor
		return
	end subroutine
	subroutine assignmentreal4(r,T2)
		real*4,intent(inout)::r
		type(PTensor),intent(in) :: T2
		r=T2%SymTensor
		return
	end subroutine
	subroutine assignmentreal8(r,T2)
		real*8,intent(inout)::r
		type(PTensor),intent(in) :: T2
		r=T2%SymTensor
		return
	end subroutine
	subroutine assignmentcom4(r,T2)
		complex(kind=4),intent(inout)::r
		type(PTensor),intent(in) :: T2
		r=T2%SymTensor
		return
	end subroutine
	subroutine assignmentcom8(r,T2)
		complex(kind=8),intent(inout)::r
		type(PTensor),intent(in) :: T2
		r=T2%SymTensor
		return
	end subroutine
	subroutine PTensorToSymTensor(T,T2)
		type(SymTensor),intent(inout) ::T
		type(PTensor),intent(in) :: T2
		T=T2%SymTensor
		return
	end subroutine
	subroutine PTensorToTensor(T,T2)
		type(Tensor),intent(inout) ::T
		type(PTensor),intent(in) :: T2
		T=T2%SymTensor
		return
	end subroutine
	subroutine SymTensorToPTensor(T,T2)
		type(PTensor),intent(inout) ::T
		type(SymTensor),intent(in) :: T2
		T%SymTensor=T2
		return
	end subroutine
	subroutine TensorToPTensor(T,T2)
		type(PTensor),intent(inout) ::T
		type(Tensor),intent(in) :: T2
		if(.not.T%getFlag())then
			call writemess('Allocate PTensor before setting value')
			call error_stop()
		end if
		T%SymTensor=T2
		return
	end subroutine
	
	subroutine allocatePBlockAll(T)
		class(PTensor),intent(inout)::T
		type(Tensor)::elementindex
		elementindex=NonZeroElement(T)
		call T%SymTensor%allocateBlock(elementindex)
		return
	end subroutine
	subroutine allocatePBlock1(T,vec)
		class(PTensor),intent(inout)::T
		integer,intent(in)::vec(:)
		logical::goon
		integer::i
		character*500,allocatable::w
		real*4,allocatable::qnum(:)
		goon=rule(T%SymDimension,vec)
		if(goon)then
			call T%SymTensor%allocateBlock(vec)
			return
		end if
		call writemess('Can not allocate the element of a PTensor')
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
	
	subroutine allocatePBlock2(T,element)
		class(PTensor),intent(inout)::T
		type(Tensor),intent(in)::element
		integer,allocatable::Adim(:)
		if(element.gt.T%getTotalData())Then
			write(*,*)"Index is larger than the len of SymTensor"
			call error_stop()
		end if
		allocate(Adim(T%getRank()))
		call IndesToaddress(T%dim(),Adim,element%ii(1))
		call allocatePBlock1(T,Adim)
		return
	end subroutine
	subroutine allocatePBlock3(T,element)
		class(PTensor),intent(inout)::T
		integer,intent(in)::element
		integer,allocatable::Adim(:)
		if(element.gt.T%getTotalData())Then
			write(*,*)"Index is larger than the len of SymTensor"
			call error_stop()
		end if
		allocate(Adim(T%getRank()))
		call IndesToaddress(T%dim(),Adim,element)
		call allocatePBlock1(T,Adim)
		return
	end subroutine
	
	subroutine set_Symelement_vec(T,vec,element) 
		class(PTensor),intent(inout)::T
		integer,intent(in)::vec(:)
		type(Tensor),intent(in)::element
		logical::goon
		integer::i
		character*500,allocatable::w
		real*4,allocatable::qnum(:)
		goon=rule(T%SymDimension,vec)
		if(goon)then
			call T%SymTensor%setValue(vec,element)
			return
		end if
		call writemess('Can not set the value to the element of a PTensor')
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
		class(PTensor),intent(inout)::T
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
	
	subroutine set_random_PTensor_all(T,region)
		class(PTensor),intent(inout)::T
		real*8,optional,intent(in)::region(2)
		integer::I
		type(Tensor)::elementindex
		elementindex=NonZeroElement(T)
		call T%SymTensor%allocateBlock(elementindex)
		call T%SymTensor%random(region)
		return
	end subroutine
	subroutine set_random_PTensor_vec(T,vec,region)
		class(PTensor),intent(inout)::T
		integer,intent(in)::vec(:)
		real*8,optional,intent(in)::region(2)
		logical::goon
		integer::i
		character*500,allocatable::w
		real*4,allocatable::qnum(:)
		goon=rule(T%SymDimension,vec)
		if(goon)then
			call T%SymTensor%random(vec,region)
			return
		end if
		call writemess('Can not set random value to the element of a PTensor')
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
	subroutine set_random_PTensor_val(T,ith,region)
		class(PTensor),intent(inout)::T
		integer,intent(in)::ith
		real*8,optional,intent(in)::region(2)
		integer,allocatable::Adim(:)
		if(ith.gt.T%getTotalData())Then
			write(*,*)"Index is larger than the len of SymTensor"
			call error_stop()
		end if
		allocate(Adim(T%getRank()))
		call IndesToaddress(T%dim(),Adim,ith)
		call set_random_PTensor_vec(T,Adim,region)
		return
	end subroutine
		
	subroutine allocateTensor2(T,dimen,typede)
		class(PTensor),intent(inout) ::T
		integer,intent(in)::dimen(:)
		integer,intent(in)::typede
		integer::length,i
		type(QuanNum),allocatable::Q(:)
		type(SymDimension)::QDimen
		length=size(dimen)
		allocate(Q(length))
		do i=1,length
			Q(i)=quantumnumber((/dimen(i),dimen(i)/))
		end do
		QDimen=Q
		call T%SymTensor%allocate(QDimen,typede)
		return
	end subroutine
	
	
	subroutine allocateTensor5(T,dimen,typede)
		class(PTensor),intent(inout) ::T
		integer,intent(in)::dimen(:)
		character(len=*),intent(in)::typede
		integer::length,i
		type(QuanNum),allocatable::Q(:)
		type(SymDimension)::QDimen
		length=size(dimen)
		allocate(Q(length))
		do i=1,length
			Q(i)=quantumnumber((/dimen(i),dimen(i)/))
		end do
		QDimen=Q
		call T%SymTensor%allocate(QDimen,typede)
		return
	end subroutine
	
	subroutine allocateTensor8(T,dimen)
		class(PTensor),intent(inout) ::T
		integer,intent(in)::dimen(:)
		integer::length,i
		type(QuanNum),allocatable::Q(:)
		type(SymDimension)::QDimen
		length=size(dimen)
		allocate(Q(length))
		do i=1,length
			Q(i)=quantumnumber((/dimen(i),dimen(i)/))
		end do
		QDimen=Q
		call T%SymTensor%allocate(QDimen)
		return
	end subroutine
	
	
	
	
	
	
	
	type(Symdimension) function  getSubDim2(T,inde)
		type(PTensor),intent(in) :: T
		integer,intent(in)::inde
		getSubDim2=T%SymDimension.subdim.inde
		return
	end function
	type(Symdimension) function  getSubDim3(T)
		type(PTensor),intent(in) :: T
		getSubDim3=T%SymDimension
		return
	end function
	type(Symdimension) function  getSubDim2_name(T,w)
		type(PTensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::w
		integer::inde
		inde=T%Nameorder(w)
		getSubDim2_name=T%SymDimension.subdim.inde
		return
	end function
	
	
	
	
	
	
	
	
	
	
	
	type(PTensor) function ParityPermutation(T,newOrder)
		type(PTensor),intent(in) :: T
		integer,intent(in)::newOrder(:)
		ParityPermutation%SymTensor=T%SymTensor.p.newOrder
		return
	end function
	type(PTensor) function ParityPermutation_name(T,newOrderchar)result(Res)
		type(PTensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::newOrderchar(:)
		Res%SymTensor=T%SymTensor.p.newOrderchar
		return
	end function
	type(PTensor) function ParityPermutation_Tensor(T,Order)result(Res)
		type(PTensor),intent(in) :: T
		type(Tensor),intent(in)::Order
		select case(Order%getType())
			case (1)
				Res=ParityPermutation(T,Order%ii())
			case (7)
				Res=ParityPermutation_name(T,Order%ai())
			case default
				call writemess('error in permutation, the data type of order')
				call error_Stop()
		end select
		return
	end function


	Type(PTensor) function ParityPermuteForWard(T,ith)Result(Res)
		Type(PTensor),intent(in)::T
		integer,intent(in)::ith
		Res%SymTensor=T%SymTensor.pf.ith
		return
	end function
	
	Type(PTensor) function ParityPermuteForWard_vec(T,vec)Result(Res)
		Type(PTensor),intent(in)::T
		integer,intent(in)::vec(:)
		Res%SymTensor=T%SymTensor.pf.vec
		return
	end function
	Type(PTensor) function ParityPermuteForWard_name(T,a)Result(Res)
		Type(PTensor),intent(in)::T
		character(len=*),intent(in)::a
		Res%SymTensor=T%SymTensor.pf.a
		return
	end function
	Type(PTensor) function ParityPermuteForWard_name_vec(T,a)Result(Res)
		Type(PTensor),intent(in)::T
		character(len=*),intent(in)::a(:)
		Res%SymTensor=T%SymTensor.pf.a
		return
	end function
	
	type(PTensor) function ParityPermuteForWard_Tensor(T,Tenorder)result(permutefo)
		type(PTensor),intent(in)::T
		type(Tensor),intent(in)::Tenorder
		permutefo%SymTensor=T.pf.Tenorder
		return
	end function
	
	
	
	
	
	Type(PTensor) function ParityPermuteBackWard(T,ith)Result(Res)
		class(PTensor),intent(in)::T
		integer,intent(in)::ith
		Res%SymTensor=T%SymTensor.pb.ith
		return
	end function
	Type(PTensor) function ParityPermuteBackWard_name(T,a)Result(Res)
		class(PTensor),intent(in)::T
		character(len=*),intent(in)::a
		integer::ith
		Res%SymTensor=T%SymTensor.pb.a
		return
	end function
	Type(PTensor) function ParityPermuteBackWard_name_vec(T,a)Result(Res)
		class(PTensor),intent(in)::T
		character(len=*),intent(in)::a(:)
		Res%SymTensor=T%SymTensor.pb.a
		return
	end function
	Type(PTensor) function ParityPermuteBackWard_vec(T,vec)Result(Res)
		class(PTensor),intent(in)::T
		integer,intent(in)::vec(:)
		Res%SymTensor=T%SymTensor.pb.vec
		return
	end function
	
	type(PTensor) function ParityPermuteBackWard_Tensor(T,Tenorder)result(permutefo)
		type(PTensor),intent(in)::T
		type(Tensor),intent(in)::Tenorder
		permutefo=T.pb.Tenorder
		return
	end function
!************************************************************************
	type(PTensor) function ProductTensor (T1,T2)
		type(PTensor),intent(in) :: T1,T2
		call ProductTensor%empty()
		call ProductTensor%ProductTensorRoutine(T1,T2)
		return
	end function
	
	type(PTensor) function divide_real8(T1,num) result(Res)
		type(PTensor),intent(in) :: T1
		real(kind=8),intent(in) ::   num
		Res%SymTensor=T1%SymTensor/num
		return
	end function
	type(PTensor) function divide_real4(T1,num) result(Res)
		type(PTensor),intent(in) :: T1
		real(kind=4),intent(in) ::   num
		Res%SymTensor=T1%SymTensor/num
		return
	end function
	type(PTensor) function divide_Tensor(T1,num) result(Res)
		type(PTensor),intent(in) :: T1
		type(Tensor),intent(in) ::   num
		Res%SymTensor=T1%SymTensor/num
		return
	end function
	
	type(PTensor) function add_PTensor(T1,T2) result(Res)
		type(PTensor),intent(in) :: T1,T2
		Res%SymTensor=T1%SymTensor+T2%SymTensor
		return
	end function
	
	type(PTensor) function minu_PTensor(T1,T2) result(Res)
		type(PTensor),intent(in) :: T1,T2
		Res%SymTensor=T1%SymTensor-T2%SymTensor
		return
	end function
	
	type(PTensor) function multiply_number_real8(T1,num) result(Res)
		type(PTensor),intent(in) :: T1
		real(kind=8),intent(in) ::   num
		Res%SymTensor=T1%SymTensor*num
		return
	end function
	type(PTensor) function multiply_number_real4(T1,num) result(Res)
		type(PTensor),intent(in) :: T1
		real(kind=4),intent(in) ::   num
		Res%SymTensor=T1%SymTensor*num
		return
	end function
	type(PTensor) function multiply_number_com4(T1,num) result(Res)
		type(PTensor),intent(in) :: T1
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
	
	type(PTensor) function contract_noName(T1_,i1,T2_,i2,len_of_contract) result(T)
		type(PTensor),intent(in) :: T1_,T2_
		integer,intent(in) :: i1(:),i2(:)
		integer,optional,intent(in)::len_of_contract(2)
		T%SymTensor=contract(T1_%SymTensor,i1,T2_%SymTensor,i2,len_of_contract)
		return
	end function
	type(PTensor) function contract_noName2(T1_,i1,T2_,i2) result(T)
		type(PTensor),intent(in) :: T1_,T2_
		integer,intent(in) :: i1,i2
		T%SymTensor=contract(T1_%SymTensor,i1,T2_%SymTensor,i2)
		return
	end function
	type(PTensor) function contract_name(T1_,name1,T2_,name2,len_of_contract) result(T)
		type(PTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:),name2(:)
		integer,optional,intent(in)::len_of_contract(2)
		T%SymTensor=contract(T1_%SymTensor,name1,T2_%SymTensor,name2,len_of_contract)
		return
	end function
	type(PTensor) function contract_name2(T1_,name1,T2_,name2) result(T)
		type(PTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1,name2
		T%SymTensor=contract(T1_%SymTensor,name1,T2_%SymTensor,name2)
		return
	end function

	type(PTensor) function contract_Same_name(T1_,T2_) result(T)
		type(PTensor),intent(in) :: T1_,T2_
		T%SymTensor=contract(T1_%SymTensor,T2_%SymTensor)
		return
	end function
	type(PTensor) function contract_name_int(T1_,name1,T2_,i2,len_of_contract) result(T)
		type(PTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:)
		integer,intent(in)::i2(:)
		integer,optional,intent(in)::len_of_contract(2)
		T%SymTensor=contract(T1_%SymTensor,name1,T2_%SymTensor,i2,len_of_contract)
		return
	end function
	type(PTensor) function contract_name_int2(T1_,name1,T2_,i2) result(T)
		type(PTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1
		integer,intent(in)::i2
	   T%SymTensor=contract(T1_%SymTensor,name1,T2_%SymTensor,i2)
		return
	end function

	type(PTensor) function contract_int_name(T1_,i1,T2_,name2,len_of_contract) result(T)
		type(PTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name2(:)
		integer,intent(in)::i1(:)
		integer,optional,intent(in)::len_of_contract(2)
		T%SymTensor=contract(T1_%SymTensor,i1,T2_%SymTensor,name2,len_of_contract)
		return
	end function
	type(PTensor) function contract_int_name2(T1_,i1,T2_,name2) result(T)
		type(PTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name2
		integer,intent(in)::i1
		T%SymTensor=contract(T1_%SymTensor,i1,T2_%SymTensor,name2)
		return
	end function
 
	type(PTensor) function Htranspose(T)
		type(PTensor),intent(in) :: T
		Htranspose%SymTensor=.H.T%SymTensor
		return
	end function
	type(PTensor) function Htranspose2(T)
		type(PTensor),intent(in) :: T
		Htranspose2%SymTensor=.Hn.T%SymTensor
		return
	end function
	
	
	type(PTensor) function conjugate(T)
		type(PTensor),intent(in) :: T
		integer :: i
		conjugate%SymTensor=.con.T%SymTensor
		return
	end function
	
	
	subroutine eye_Tensor3(T,m,n,classtype)
		class(PTensor),intent(inout)::T
		integer,intent(in)::m,n
		character(len=*),intent(in)::classtype
		integer::i
		call T%empty()
		call T%allocate((/m,n/),classtype)
		call T%SymTensor%eye()
		return
	end subroutine
	subroutine eye_Tensor4(T,m,n,classtype)
		class(PTensor),intent(inout)::T
		integer,intent(in)::m,n,classtype
		integer::i
		call T%empty()
		call T%allocate((/m,n/),classtype)
		call T%SymTensor%eye()
		return
	end subroutine
	subroutine eye_Tensor5(T,dimen,classtype)
		class(PTensor),intent(inout)::T
		character(len=*),intent(in)::classtype
		integer,intent(in)::dimen(2)
		call T%empty()
		call T%allocate(dimen,classtype)
		call T%SymTensor%eye()
		return
	end subroutine
	
	subroutine eye_Tensor6(T,dimen,classtype)
		class(PTensor),intent(inout)::T
		character(len=*),intent(in)::classtype
		type(Symdimension),intent(in)::dimen
		call T%SymTensor%eye(dimen,classtype)
		return
	end subroutine
	
	
	type(PTensor) function expmPTensor(T)
		type(PTensor),intent(in) :: T
		expmPTensor%SymTensor=expm(T%SymTensor)
		return
	end function
	
	type(PTensor) function ParityFusePTensor(B,outDim,Order,row_)
		class(PTensor),intent(in)::B
		type(SymDimension),intent(inout)::outDim
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
			NewQuanNum=fuseOrder(order,B%quantumnumber(1),B%quantumnumber(2))
			ParityFusePTensor%SymTensor=B%QuanFuse(NewQuanNum,Order,.true.)
		else
			rank=B%getRank()
			i=rank-1
			outDim=B%SymDimension.subdim.[i,rank]
			NewQuanNum=fuseOrder(order,B%quantumnumber(i),B%quantumnumber(rank))
			ParityFusePTensor%SymTensor=B%QuanFuse(NewQuanNum,Order,.false.)
		end if
		return
	end function
	
	
	
	type(PTensor) function ParityFusePTensor2(B,row_)result(ParityFusePTensor)
		class(PTensor),intent(in)::B
		type(Tensor)::Order
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
			NewQuanNum=fuseOrder(order,B%quantumnumber(1),B%quantumnumber(2))
			ParityFusePTensor%SymTensor=B%QuanFuse(NewQuanNum,Order,.true.)
		else
			rank=B%getRank()
			i=rank-1
			NewQuanNum=fuseOrder(order,B%quantumnumber(i),B%quantumnumber(rank))
			ParityFusePTensor%SymTensor=B%QuanFuse(NewQuanNum,Order,.false.)
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
	type(PTensor) function ParityFusePTensor3(B,ith,outDimen,outOrder,orderinfo,row_)result(T)
		class(PTensor),intent(in)::B
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
				T=ParityFusePTensor(B,dimen,order,.true.)
				outDimen=dimen
				outOrder=order
				if(order%getRank().eq.0)then
					call writemess('ERROR in ParityFusePTensor3',-1)
					call error_stop
				end if
				if(order%getRank().eq.1)then
					call orderinfo%setValue(ith-1,1)
				else
					call orderinfo%setValue(ith-1,order%dim(1))
				end if
				do i=2,ith-1
					T=ParityFusePTensor(T,dimen,order,.true.)
					outDimen=outDimen+dimen
					if(order%getRank().eq.1)then
						call orderinfo%setValue(ith-i,1)
					else
						call orderinfo%setValue(ith-i,order%dim(1))
					end if
					outOrder=order.Rpaste.outOrder
				end do
			else
				T=B
			end if
		else
			rank=B%getRank()
			rankV=rank-ith+1
			if(rankV.gt.1)then
				call orderinfo%allocate((/rankV-1/),'integer')
				T=ParityFusePTensor(B,dimen,order,.false.)
				outDimen=dimen
				outOrder=order
				if(order%getRank().eq.0)then
					call writemess('ERROR in ParityFusePTensor3',-1)
					call error_stop
				end if
				if(order%getRank().eq.1)then
					call orderinfo%setValue(rankV-1,1)
				else
					call orderinfo%setValue(rankV-1,order%dim(1))
				end if
				do i=2,rankV-1
					T=ParityFusePTensor(T,dimen,order,.false.)
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
	
	type(PTensor) function ParityFusePTensor4(B,ith,row_)result(T)
		class(PTensor),intent(in)::B
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
				T=ParityFusePTensor2(B,.true.)
				do i=2,ith-1
					T=ParityFusePTensor2(T,.true.)
				end do
			else
				T=B
			end if
		else
			rank=B%getRank()
			rankV=rank-ith+1
			if(rankV.gt.1)then
				T=ParityFusePTensor2(B,.false.)
				do i=2,rankV-1
					T=ParityFusePTensor2(T,.false.)
				end do
			else
				T=B
			end if
		end if
		
		return
	end function
	
	
	
	type(PTensor) function ParitySplitPTensor(B,NewQuanDimen,Order,row_)
		class(PTensor),intent(in)::B
		type(SymDimension),intent(in)::NewQuanDimen
		type(QuanNum)::NewQ1,NewQ2
		type(Tensor),intent(in)::Order
		logical,optional,intent(in)::row_
		logical::row
		if(present(row_))then
			row=row_
		else
			row=.true.
		end if
		NewQ1=NewQuanDimen%QuantumNumber(1)
		NewQ2=NewQuanDimen%QuantumNumber(2)
		if(row)then
			ParitySplitPTensor%SymTensor=B%QuanSplit(NewQ1,NewQ2,Order,.true.)
		else
			ParitySplitPTensor%SymTensor=B%QuanSplit(NewQ1,NewQ2,Order,.false.)
		end if
		if(NewQuanDimen%outNameFlag().eq.1)then
			if(row)then
				call ParitySplitPTensor%setName(1,NewQuanDimen%outName(1))
				call ParitySplitPTensor%setName(2,NewQuanDimen%outName(2))
			else
				call ParitySplitPTensor%setName(ParitySplitPTensor%getRank()-1,NewQuanDimen%outName(1))
				call ParitySplitPTensor%setName(ParitySplitPTensor%getRank(),NewQuanDimen%outName(2))
			end if
		end if
		return
	end function
	type(PTensor) function ParitySplitPTensor2(B,inDimen,inOrder,orderinfo,row)result(T)
		class(PTensor),intent(in)::B
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
			call writemess('ERROR in Split PTensor, input parameter error')
			call error_stop()
		end if
		orderindex=1
		dimen=inDimen.subdim.[rank-1,rank]
		istart=1
		iend=orderinfo%ii(orderindex)
		order=inOrder%subTensor((/-1,istart,iend/))
		istart=iend+1
		T=ParitySplitPTensor(B,dimen,Order,row)
		do i=rank-2,2,-2
			j=i-1
			orderindex=orderindex+1
			dimen=inDimen.subdim.[i-1,i]
			iend=istart+orderinfo%ii(orderindex)-1
			order=inOrder%subTensor((/-1,istart,iend/))
			istart=iend+1
			T=ParitySplitPTensor(T,dimen,order,row)
		end do
		return
	end function
	
	
	function SVD(T,QuanNumCut)
		type(PTensor),allocatable::SVD(:)
		class(PTensor),intent(in)::T
		type(QuanNum),optional,intent(in)::QuanNumCut
		allocate(SVD(3))
		call T%SVDSymTensorRoutine(SVD(1)%SymTensor,SVD(2)%SymTensor,SVD(3)%SymTensor,QuanNumCut)
		return
	end function
	
	subroutine SVDRoutine(T,U,S,V,QuanNumCut)
		class(PTensor),intent(in)::T
		type(PTensor),intent(inout)::U,S,V
		type(QuanNum),optional,intent(in)::QuanNumCut
		call T%SVDSymTensorRoutine(U%SymTensor,S%SymTensor,V%SymTensor,QuanNumCut)
		return
	end subroutine
	function SVDName(Tin,nameU,nameV,QuanNumCut)
		type(PTensor),allocatable::SVDName(:)
		class(PTensor),intent(in)::Tin
		type(QuanNum),optional,intent(in)::QuanNumCut
		character(len=*)::nameU,nameV
		type(PTensor)::T
		integer::rankU,rankV,i,j,rank
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
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
			call writemess("ERROR in SVDcutoff_name")
			call writemess(rankU+','+rankV+','+rank)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name")
			call writemess(nameU)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name")
			call writemess(nameV)
			call error_stop()
		end if
		T=T%ParityFuse(rankU,dimen(1),order(1),orderinfo(1),.true.)
		T=T%ParityFuse(2,dimen(2),order(2),orderinfo(2),.false.)
		allocate(SVDName(3))
		SVDName=SVD(T,QuanNumCut)
		if(Tin%getNameFlag().ne.0)then
			call SVDName(1)%setName(2,'SVD.U')
			call SVDName(2)%setName(1,'SVD.s1')
			call SVDName(2)%setName(2,'SVD.s2')
			call SVDName(2)%setName(1,'SVD.V')
		end if
		SVDName(1)=SVDName(1)%ParitySplit(dimen(1),order(1),orderinfo(1),.true.)
		SVDName(3)=SVDName(3)%ParitySplit(dimen(2),order(2),orderinfo(2),.false.)
		return
	end function 
	function SVDName_back(Tin,nameU,nameV,QuanNumCut)result(SVDName)
		type(PTensor),allocatable::SVDName(:)
		class(PTensor),intent(in)::Tin
		type(QuanNum),optional,intent(in)::QuanNumCut
		character(len=*)::nameU,nameV
		type(PTensor)::T
		integer::rankU,rankV,i,j,rank
		type(QuanNum),allocatable::QLeft(:),QRight(:)
		type(QuanNum),allocatable::QNewLeft(:),QNewRight(:)
		character(len=50),allocatable::Lname(:),Rname(:)
		type(Tensor),allocatable::Lorder(:),Rorder(:)
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
			call writemess("ERROR in SVDcutoff_name")
			call writemess(rankU+','+rankV+','+rank)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name")
			call writemess(nameU)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name")
			call writemess(nameV)
			call error_stop()
		end if
		allocate(QLeft(rankU))
		allocate(Lname(rankU))
		allocate(QRight(rankV))
		allocate(Rname(rankV))
		if(rankU.gt.1)then
			allocate(Lorder(rankU-1))
			allocate(QNewLeft(rankU-1))
		end if
		if(rankV.gt.1)then
			allocate(Rorder(rankV-1))
			allocate(QNewRight(rankV-1))
		end if
		do i=1,rankU
			QLeft(i)=T%quantumnumber(i)
			Lname(i)=T%outName(i)
		end do
		do i=1,rankV
			QRight(i)=T%quantumnumber(rankU+i)
			Rname(i)=T%outName(rankU+i)
		end do
		
		if(rankU.gt.1)then
			QNewLeft(1)=fuseOrder(Lorder(1),QLeft(1),QLeft(2))
			T=T%QuanFuse(QNewLeft(1),Lorder(1),.true.)
			do i=2,rankU-1
				QNewLeft(i)=fuseOrder(Lorder(i),QNewLeft(i-1),QLeft(i+1))
				T=T%QuanFuse(QNewLeft(i),Lorder(i),.true.)
			end do
		end if
		
		if(rankV.gt.1)then
			i=1
			QNewRight(i)=fuseOrder(Rorder(i),QRight(rankV-1),QRight(rankV))
			T=T%QuanFuse(QNewRight(i),Rorder(i),.false.)
			j=rankV-2
			do i=2,rankV-1
				QNewRight(i)=fuseOrder(Rorder(i),QRight(j),QNewRight(i-1))
				T=T%QuanFuse(QNewRight(i),Rorder(i),.false.)
				j=j-1
			end do
		end if
		
		allocate(SVDName(3))
		SVDName=SVD(T,QuanNumCut)
		if(rankU.gt.1)then
			do i=rankU-2,1,-1
				SVDName(1)=SVDName(1)%QuanSplit(QNewLeft(i),QLeft(i+2),Lorder(i+1),.true.)
			end do
			SVDName(1)=SVDName(1)%QuanSplit(QLeft(1),QLeft(2),Lorder(1),.true.)
		end if
		do i=1,rankU
			call SVDName(1)%setName(i,Lname(i))
		end do
		if(rankV.gt.1)then
			j=1
			do i=rankV-2,1,-1
				SVDName(3)=SVDName(3)%QuanSplit(Qright(j),QNewright(i),Rorder(i+1),.false.)
				j=j+1
			end do
			SVDName(3)=SVDName(3)%QuanSplit(Qright(rankV-1),Qright(rankV),Rorder(1),.false.)
		end if
		do i=1,rankV
			call SVDName(3)%setName(i+1,Rname(i))
		end do
		return
	end function
	
	
	subroutine QRdecompositionPTensor(T,Q,R) 
		class(PTensor),intent(in)::T
		type(PTensor),intent(inout)::Q,R
		call T%SymTensor%SymQRRoutine(Q%SymTensor,R%SymTensor)
		return
	end subroutine
	function QRPTensor_noname(T)result(Res)
		type(PTensor),allocatable::Res(:)
		class(PTensor),intent(in)::T
		allocate(Res(2))
		call QRdecompositionPTensor(T,Res(1),Res(2)) 
		return
	end function
	
	
	
	function QRPTensor_name(Tin,nameU,nameV)result(Res)
		type(PTensor),allocatable::Res(:)
		character(len=*),intent(in)::nameU,nameV
		class(PTensor),intent(in)::Tin
		integer::rankU,rankV,i,j,rank
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
		Type(PTensor)::T
		rank=Tin%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if(Tin%outTensorName(i).equ.nameU) rankU=rankU+1
			if(Tin%outTensorName(i).equ.nameV) rankV=rankV+1
		end do
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in QRPTensor_name")
			call writemess(rankU+','+rankV+','+rank)
			call Tin%diminfo()
			call writemess(nameU+','+nameV)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			T=Tin.pf.nameU
		else
			T=Tin.pb.nameV
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in QRPTensor_name,no such name")
			call writemess(nameU)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in QRPTensor_name,no such name")
			call writemess(nameV)
			call error_stop()
		end if
		T=T%ParityFuse(rankU,dimen(1),order(1),orderinfo(1),.true.)
		T=T%ParityFuse(2,dimen(2),order(2),orderinfo(2),.false.)
		allocate(Res(2))
		call QRdecompositionPTensor(T,Res(1),Res(2)) 
		if(Tin%getNameFlag().ne.0)then
			call Res(1)%setName(2,'QR.Q')
			call Res(2)%setName(1,'QR.R')
		end if
		Res(1)=Res(1)%ParitySplit(dimen(1),order(1),orderinfo(1),.true.)
		Res(2)=Res(2)%ParitySplit(dimen(2),order(2),orderinfo(2),.false.)
	end function
	
	function QRPTensor_name_back(Tin,nameU,nameV)result(Res)
		type(PTensor),allocatable::Res(:)
		character(len=*),intent(in)::nameU,nameV
		class(PTensor),intent(in)::Tin
		integer::rankU,rankV,i,j,rank
		type(QuanNum),allocatable::QLeft(:),QRight(:)
		type(QuanNum),allocatable::QNewLeft(:),QNewRight(:)
		character(len=50),allocatable::Lname(:),Rname(:)
		type(Tensor),allocatable::Lorder(:),Rorder(:)
		Type(PTensor)::T
		rank=Tin%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if(Tin%outTensorName(i).equ.nameU) rankU=rankU+1
			if(Tin%outTensorName(i).equ.nameV) rankV=rankV+1
		end do
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in QRPTensor_name")
			call writemess(rankU+','+rankV+','+rank)
			call Tin%diminfo()
			call writemess(nameU+','+nameV)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			T=Tin.pf.nameU
		else
			T=Tin.pb.nameV
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in QRPTensor_name,no such name")
			call writemess(nameU)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in QRPTensor_name,no such name")
			call writemess(nameV)
			call error_stop()
		end if
		allocate(QLeft(rankU))
		allocate(Lname(rankU))
		allocate(QRight(rankV))
		allocate(Rname(rankV))
		if(rankU.gt.1)then
			allocate(Lorder(rankU-1))
			allocate(QNewLeft(rankU-1))
		end if
		if(rankV.gt.1)then
			allocate(Rorder(rankV-1))
			allocate(QNewRight(rankV-1))
		end if
		do i=1,rankU
			QLeft(i)=T%quantumnumber(i)
			Lname(i)=T%outName(i)
		end do
		do i=1,rankV
			QRight(i)=T%quantumnumber(rankU+i)
			Rname(i)=T%outName(rankU+i)
		end do
		
		if(rankU.gt.1)then
			QNewLeft(1)=fuseOrder(Lorder(1),QLeft(1),QLeft(2))
			T=T%QuanFuse(QNewLeft(1),Lorder(1),.true.)
			do i=2,rankU-1
				QNewLeft(i)=fuseOrder(Lorder(i),QNewLeft(i-1),QLeft(i+1))
				T=T%QuanFuse(QNewLeft(i),Lorder(i),.true.)
			end do
		end if
		
		if(rankV.gt.1)then
			i=1
			QNewRight(i)=fuseOrder(Rorder(i),QRight(rankV-1),QRight(rankV))
			T=T%QuanFuse(QNewRight(i),Rorder(i),.false.)
			j=rankV-2
			do i=2,rankV-1
				QNewRight(i)=fuseOrder(Rorder(i),QRight(j),QNewRight(i-1))
				T=T%QuanFuse(QNewRight(i),Rorder(i),.false.)
				j=j-1
			end do
		end if
		
		
		allocate(Res(2))
		call QRdecompositionPTensor(T,Res(1),Res(2)) 
		if(rankU.gt.1)then
			do i=rankU-2,1,-1
				Res(1)=Res(1)%QuanSplit(QNewLeft(i),QLeft(i+2),Lorder(i+1),.true.)
			end do
			Res(1)=Res(1)%QuanSplit(QLeft(1),QLeft(2),Lorder(1),.true.)
		end if
		do i=1,rankU
			call Res(1)%setName(i,Lname(i))
		end do
		if(rankV.gt.1)then
			j=1
			do i=rankV-2,1,-1
				Res(2)=Res(2)%QuanSplit(Qright(j),QNewright(i),Rorder(i+1),.false.)
				j=j+1
			end do
			Res(2)=Res(2)%QuanSplit(Qright(rankV-1),Qright(rankV),Rorder(1),.false.)
		end if
		do i=1,rankV
			call Res(2)%setName(i+1,Rname(i))
		end do
		
		return
	end function
	
	subroutine LQdecompositionPTensor(T,L,Q) 
		class(PTensor),intent(in)::T
		type(PTensor),intent(inout)::L,Q
		call T%SymTensor%SymLQRoutine(L%SymTensor,Q%SymTensor)
		return
	end subroutine
	function LQPTensor_noname(T)result(Res)
		type(PTensor),allocatable::Res(:)
		class(PTensor),intent(in)::T
		allocate(Res(2))
		call LQdecompositionPTensor(T,Res(1),Res(2)) 
		return
	end function
	function LQPTensor_name(Tin,nameU,nameV)result(Res)
		type(PTensor),allocatable::Res(:)
		character(len=*),intent(in)::nameU,nameV
		class(PTensor),intent(in)::Tin
		integer::rankU,rankV,i,j,rank
		Type(PTensor)::T
		type(Tensor)::order(2),orderinfo(2)
		type(SymDimension)::dimen(2)
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
			call writemess("ERROR in LQPTensor_name")
			call writemess(rankU+','+rankV+','+rank)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in LQPTensor_name,no such name")
			call writemess(nameU)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in LQPTensor_name,no such name")
			call writemess(nameV)
			call error_stop()
		end if
		T=T%ParityFuse(rankU,dimen(1),order(1),orderinfo(1),.true.)
		T=T%ParityFuse(2,dimen(2),order(2),orderinfo(2),.false.)
		allocate(Res(2))
		call LQdecompositionPTensor(T,Res(1),Res(2)) 
		if(Tin%getNameFlag().ne.0)then
			call Res(1)%setName(2,'LQ.L')
			call Res(2)%setName(1,'LQ.Q')
		end if
		Res(1)=Res(1)%ParitySplit(dimen(1),order(1),orderinfo(1),.true.)
		Res(2)=Res(2)%ParitySplit(dimen(2),order(2),orderinfo(2),.false.)
		return
	end function
	
	function LQPTensor_name_back(Tin,nameU,nameV)result(Res)
		type(PTensor),allocatable::Res(:)
		character(len=*),intent(in)::nameU,nameV
		class(PTensor),intent(in)::Tin
		integer::rankU,rankV,i,j,rank
		type(QuanNum),allocatable::QLeft(:),QRight(:)
		type(QuanNum),allocatable::QNewLeft(:),QNewRight(:)
		character(len=50),allocatable::Lname(:),Rname(:)
		type(Tensor),allocatable::Lorder(:),Rorder(:)
		Type(PTensor)::T
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
			call writemess("ERROR in LQPTensor_name")
			call writemess(rankU+','+rankV+','+rank)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in LQPTensor_name,no such name")
			call writemess(nameU)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in LQPTensor_name,no such name")
			call writemess(nameV)
			call error_stop()
		end if
		allocate(QLeft(rankU))
		allocate(Lname(rankU))
		allocate(QRight(rankV))
		allocate(Rname(rankV))
		if(rankU.gt.1)then
			allocate(Lorder(rankU-1))
			allocate(QNewLeft(rankU-1))
		end if
		if(rankV.gt.1)then
			allocate(Rorder(rankV-1))
			allocate(QNewRight(rankV-1))
		end if
		do i=1,rankU
			QLeft(i)=T%quantumnumber(i)
			Lname(i)=T%outName(i)
		end do
		do i=1,rankV
			QRight(i)=T%quantumnumber(rankU+i)
			Rname(i)=T%outName(rankU+i)
		end do
		
		if(rankU.gt.1)then
			QNewLeft(1)=fuseOrder(Lorder(1),QLeft(1),QLeft(2))
			T=T%QuanFuse(QNewLeft(1),Lorder(1),.true.)
			do i=2,rankU-1
				QNewLeft(i)=fuseOrder(Lorder(i),QNewLeft(i-1),QLeft(i+1))
				T=T%QuanFuse(QNewLeft(i),Lorder(i),.true.)
			end do
		end if
		
		if(rankV.gt.1)then
			i=1
			QNewRight(i)=fuseOrder(Rorder(i),QRight(rankV-1),QRight(rankV))
			T=T%QuanFuse(QNewRight(i),Rorder(i),.false.)
			j=rankV-2
			do i=2,rankV-1
				QNewRight(i)=fuseOrder(Rorder(i),QRight(j),QNewRight(i-1))
				T=T%QuanFuse(QNewRight(i),Rorder(i),.false.)
				j=j-1
			end do
		end if
		
		
		allocate(Res(2))
		call LQdecompositionPTensor(T,Res(1),Res(2)) 
		if(rankU.gt.1)then
			do i=rankU-2,1,-1
				Res(1)=Res(1)%QuanSplit(QNewLeft(i),QLeft(i+2),Lorder(i+1),.true.)
			end do
			Res(1)=Res(1)%QuanSplit(QLeft(1),QLeft(2),Lorder(1),.true.)
		end if
		do i=1,rankU
			call Res(1)%setName(i,Lname(i))
		end do
		if(rankV.gt.1)then
			j=1
			do i=rankV-2,1,-1
				Res(2)=Res(2)%QuanSplit(Qright(j),QNewright(i),Rorder(i+1),.false.)
				j=j+1
			end do
			Res(2)=Res(2)%QuanSplit(Qright(rankV-1),Qright(rankV),Rorder(1),.false.)
		end if
		do i=1,rankV
			call Res(2)%setName(i+1,Rname(i))
		end do
		
		return
	end function
		
	
	
	type(PTensor) function directProductTensor(T1,T2)result(Res)
		type(PTensor),intent(in) :: T1,T2
		Res%SymTensor=T1%SymTensor.kron.T2%SymTensor
		return
	end function
	
	
	type(Tensor) function TdotTensor(phi1,phi2)result(dotTensor)
		Type(PTensor),intent(in)::phi1,phi2
		dotTensor=phi1%SymTensor.x.phi2%SymTensor
		RETURN
	end function
	
	


! reorder the Tensor, order the leg as the order in allTensorName(:)
	subroutine  DimOrder(inoutT,allTensorName)
		class(PTensor),intent(inout)::inoutT
		character(len=*),intent(in)::allTensorName(:)
		character(len=max_len_of_char_in_TData),allocatable :: TensorName(:)
		character(len=max_len_of_char_in_TData),allocatable :: Names(:)
		integer::rank,i
		if(.not.inoutT%getFlag())then
			call writemess('There is no data in the PTensor')
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
		call writemess('ERROR in checkOrder')
		call error_stop
	end function



!**********************************************************************
!**********************************************************************
!	the code below is for MPI
!**********************************************************************
	subroutine MPI_send_PTensor(Ten1,Ten2,ID1,ID2,ierr,MPIcommon)
		type(PTensor),intent(in)::Ten1
		type(PTensor),intent(inout)::Ten2
		integer,intent(in)::ID1,ID2,ierr
		integer,optional,intent(in)::MPIcommon
		call MPI_send_SymTensor(Ten1%SymTensor,Ten2%SymTensor,ID1,ID2,ierr,MPIcommon)
		return
	end subroutine
	subroutine MPI_BCAST_PTensor(Ten1,ID,ierr,MPIcommon)
		type(PTensor),intent(inout)::Ten1
		integer,intent(in)::ID,ierr
		integer,optional,intent(in)::MPIcommon
		call MPI_BCAST_SymTensor(Ten1%SymTensor,ID,ierr,MPIcommon)
		return
	end subroutine



	subroutine MPI_SUM_PTensor1(inoutTensor,ierr,MPIcommon)
		type(PTensor),intent(inout)::inoutTensor
		integer,intent(in)::ierr
		integer,optional,intent(in)::MPIcommon
		call MPI_SUM_SymTensor(inoutTensor%SymTensor,ierr,MPIcommon)
		return
	end subroutine



end module
