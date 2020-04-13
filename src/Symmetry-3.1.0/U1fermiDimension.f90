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

module U1fermiDimension_Type
	use dimension_typede
	use SymDimension_typede
	use Tools
	use mpi
	implicit none
	private
	
	
	public::U1fDimension
	type,extends (SymDimension) :: U1fDimension
		integer,allocatable :: fermiArrow(:)
		logical::ArrowFlag=.false.
	contains
		procedure,public::getArrowFlag=>FermiArrowFlag
		generic,public::getU1Rule=>getU1DimRule,getU1DimRuleChar,getU1DimAllRule
		generic,public::getFermiArrow=>getFermiArrow,getFermiArrowChar,getFermiArrowAllRule
		procedure::getU1DimRule,getU1DimRuleChar,getU1DimAllRule
		procedure::getFermiArrow,getFermiArrowChar,getFermiArrowAllRule
		procedure::getSymDimRule,getSymDimRuleChar,getSymDimAllRule

	end type SymDimension
	
	public::assignment(=)
	interface assignment(=)
		module procedure U1fDimInitialization!U1fdimension=Symdimension
		module procedure U1fDimInitialization2!U1fdimension=integer(:)
		module procedure U1fDimInitialization3!U1fdimension=dimension
		module procedure U1fDimInitialization4!dimension=U1fdimension
		module procedure U1fDimInitialization5!U1fdimension= (/ QuanNum /)
		module procedure U1fDimInitialization6!U1fdimension=Symdimension
		module procedure U1fQuanNumInitialization2
	end interface
	public::operator(.fuse.)
	interface operator(.fuse.)!overwrite the function in type(Dimension)
		module procedure fuseDimension_val
		module procedure fuseDimension_vec
	end interface
	public::operator(.split.)!overwrite the function in type(Dimension)
	interface operator(.split.)
		module procedure splitDimension2
		module procedure splitDimension3
		module procedure splitDimensionAll
	end interface
	public::operator(.sub.)!overwrite the function in type(Dimension)
	interface operator(.sub.)
		module procedure getSubDim2
		module procedure getSubDim3
		module procedure getSubDim4
		module procedure getSubDim2_name
	end interface
	public::operator(.subdim.)!overwrite the function in type(Dimension)
	interface operator(.subdim.)
		module procedure getSubDim2
		module procedure getSubDim3
		module procedure getSubDim4
		module procedure getSubDim2_name
	end interface
	
	
	public::operator(+)
	interface operator(+)
		module procedure Dimadd
		module procedure DimaddQuanNum
		module procedure QuanNumAddDim
		module procedure QuanNumAddQuanNum
	end interface
	
	public::U1fdimpermute_forwards,U1fDimpermute_backwards,U1fDimpermute_forwards_index,U1fDimpermute_backwards_index
	public::MPI_send_U1fDimension,MPI_BCAST_U1fDimension
contains
	logical function FermiArrowFlag()
		class(U1fDimension),intent(in) :: Dimen
		FermiArrowFlag=Dimen%getArrowFlag
		return
	end function
	
	
	
	
! subroutines and funcitons for SymDimension

	integer function getSymDimRule(Dimen,ith)
		class(SymDimension),intent(in) :: Dimen
		integer,intent(in)::ith
		call writemess("This is the type of U(1) fermi diemsnion")
		call writemess("getU1Rule or getfermiArrow?")
		call error_stop()
	end function
	integer function getSymDimRuleChar(Dimen,w)
		class(SymDimension),intent(in) :: Dimen
		character(len=*),intent(in)::w
		call writemess("This is the type of U(1) fermi diemsnion")
		call writemess("getU1Rule or getfermiArrow?")
		call error_stop()
	end function
	function getSymDimAllRule(Dimen)
		integer,allocatable::getSymDimAllRule(:)
		class(SymDimension),intent(in) :: Dimen
		call writemess("This is the type of U(1) fermi diemsnion")
		call writemess("getU1Rule or getfermiArrow?")
		call error_stop()
	end function
	integer function getU1DimRule(Dimen,ith)
		class(SymDimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(.not.Dimen%getQNflag())then
			call writemess("There is no Quantum Number")
			call error_stop()
		end if
		if(ith.gt.Dimen%outlenDimData())then
			if(Dimen%outlenDimData().eq.0)then
				call writemess('There is no Data in the U1fdimension')
				call error_stop()
			end if
			call writemess("The input index is larger then the length of the dimension")
			call writemess("in getting Rule of the dimension")
			call error_stop()
		end if
		getU1DimRule=Dimen%QN(ith)%getRule()
		return
	end function
	integer function getU1DimRuleChar(Dimen,w)
		class(SymDimension),intent(in) :: Dimen
		character(len=*),intent(in)::w
		integer::ith
		if(.not.Dimen%getQNflag())then
			call writemess("There is no Quantum Number")
			call error_stop()
		end if
		ith=Dimen%NameOrder(w)
		if(ith.gt.Dimen%outlenDimData())then
			if(Dimen%outlenDimData().eq.0)then
				call writemess('There is no Data in the Symdimension')
				call error_stop()
			end if
			call writemess("The input index is larger then the length of the dimension")
			call writemess("in getting Rule of the dimension")
			call error_stop()
		end if
		getU1DimRuleChar=Dimen%QN(ith)%getRule()
		return
	end function
	function getU1DimAllRule(Dimen)
		integer,allocatable::getSymDimAllRule(:)
		class(SymDimension),intent(in) :: Dimen
		integer::i,rank
		if(.not.Dimen%getQNflag())then
			call writemess("There is no Quantum Number")
			call error_stop()
		end if
		rank=Dimen%outlenDimData()
		allocate(getSymDimAllRule(rank))
		do i=1,rank
			getU1DimAllRule(i)=Dimen%QN(i)%getRule()
		end do
		return
	end function
	integer function getFermiArrow(Dimen,ith)
		class(SymDimension),intent(in) :: Dimen
		integer,intent(in)::ith
		if(.not.Dimen%getArrowFlag())then
			call writemess("There is no fermi-arrow")
			call error_stop()
		end if
		if(ith.gt.Dimen%outlenDimData())then
			if(Dimen%outlenDimData().eq.0)then
				call writemess('There is no Data in the U1fdimension')
				call error_stop()
			end if
			call writemess("The input index is larger then the length of the dimension")
			call writemess("in getting Rule of the dimension")
			call error_stop()
		end if
		getFermiArrow=Dimen%fermiArrow(ith)
		return
	end function
	integer function getFermiArrowChar(Dimen,w)
		class(SymDimension),intent(in) :: Dimen
		character(len=*),intent(in)::w
		integer::ith
		if(.not.Dimen%getArrowFlag())then
			call writemess("There is no fermi-arrow")
			call error_stop()
		end if
		ith=Dimen%NameOrder(w)
		if(ith.gt.Dimen%outlenDimData())then
			if(Dimen%outlenDimData().eq.0)then
				call writemess('There is no Data in the Symdimension')
				call error_stop()
			end if
			call writemess("The input index is larger then the length of the dimension")
			call writemess("in getting Rule of the dimension")
			call error_stop()
		end if
		getFermiArrowChar=Dimen%fermiArrow(ith)
		return
	end function
	function getFermiArrowAllRule(Dimen)
		integer,allocatable::getFermiArrowAllRule(:)
		class(SymDimension),intent(in) :: Dimen
		integer::i,rank
		if(.not.Dimen%getArrowFlag())then
			call writemess("There is no fermi-arrow")
			call error_stop()
		end if
		rank=Dimen%outlenDimData()
		allocate(getFermiArrowAllRule(rank))
		do i=1,rank
			getFermiArrowAllRule(i)=Dimen%fermiArrow(i)
		end do
		return
	end function




	subroutine U1fDimInitialization(inoutQN,inQN)
		type(U1fDimension),intent(inout) ::inoutQN
		type(U1fDimension),intent(in) :: inQN
		integer::i,lenDimData
		inoutQN%SymDimension=inQN%SymDimension
		inoutQN%ArrowFlag=inQN%getArrowFlag
		if(inoutQN%getArrowFlag)then
			lenDimData=inQN%outlenDimData()
			call allocateCheck(inoutQN%fermiArrow,lenDimData)
			do i=1,lenDimData
				inoutQN%fermiArrow(i)=inQN%fermiArrow(i)
			end do
		end if
		return
	end subroutine
	subroutine U1fDimInitialization2(inoutQN,dimdata)
		type(U1fDimension),intent(inout) ::inoutQN
		integer,intent(in) :: dimdata(:)
		inoutQN%Symdimension=dimdata
		inoutQN%ArrowFlag=.false.
		return
	end subroutine
	subroutine U1fDimInitialization3(inoutQN,inQN)
		type(U1fDimension),intent(inout) ::inoutQN
		type(Dimension),intent(in) :: inQN
		integer::i,lenDimData
		inoutQN%SymDimension=inQN
		inoutQN%ArrowFlag=.false.
		return
	end subroutine
	subroutine U1fDimInitialization4(inoutQN,inQN)
		type(Dimension),intent(inout) ::inoutQN
		type(U1fDimension),intent(in) :: inQN
		integer::i,lenDimData
		inoutQN=inQN%dimension
		return
	end subroutine
	subroutine U1fDimInitialization5(inoutQN,inQN)
		type(U1fDimension),intent(inout) ::inoutQN
		type(QuanNum),intent(in) :: inQN(:)
		inoutQN%SymDimension=inQN
		inoutQN%ArrowFlag=.false.
		return
	end subroutine
	subroutine U1fDimInitialization6(inoutQN,inQN)
		type(U1fDimension),intent(inout) ::inoutQN
		type(SymDimension),intent(in) :: inQN
		inoutQN%SymDimension=inQN
		inoutQN%ArrowFlag=.false.
		return
	end subroutine
	subroutine U1fQuanNumInitialization2(dimen,inQN)
		type(U1fDimension),intent(inout) ::dimen
		type(QuanNum),intent(in) :: inQN
		inoutQN%SymDimension=inQN
		inoutQN%ArrowFlag=.false.
		return
	end subroutine

	subroutine emptyDimension(Dimen)!overwrite the subroutine in dimension.f90
		class(U1fDimension),intent(inout) ::Dimen
		integer::i
		call Dimen%Symdimension%empty()
		if(allocated(Dimen%fermiArrow))then
			if(deallocate_memory_flag)deallocate(Dimen%fermiArrow)
		end if
		Dimen%ArrowFlag=.false.
		return
	end subroutine
	subroutine cleanDimension(Dimen)!overwrite the subroutine in dimension.f90
		class(U1fDimension),intent(inout) ::Dimen
		integer::i
		call Dimen%Symdimension%deallocate()
		if(allocated(Dimen%fermiArrow))then
			deallocate(Dimen%fermiArrow)
		end if
		Dimen%ArrowFlag=.false.
		return
	end subroutine


	subroutine setSymDimensionRule(inoutQN,ith,rule)
		class(U1fDimension),intent(inout) ::inoutQN
		integer,intent(in)::ith,rule
		call writemess("This is the type of U(1) fermi diemsnion")
		call writemess("setU1Rule or setfermiArrow?")
		call error_stop()
	end subroutine
	subroutine setSymDimensionRule_char(inoutQN,w,rule)
		class(U1fDimension),intent(inout) ::inoutQN
		character(len=*),intent(in)::w
		integer,intent(in)::rule
		call writemess("This is the type of U(1) fermi diemsnion")
		call writemess("setU1Rule or setfermiArrow?")
		call error_stop()
	end subroutine
	subroutine setAllSymDimensionRule(inoutQN,rule)
		class(U1fDimension),intent(inout) ::inoutQN
		integer,intent(in)::rule(:)
		call writemess("This is the type of U(1) fermi diemsnion")
		call writemess("setU1Rule or setfermiArrow?")
		call error_stop()
	end subroutine
	


	subroutine setU1DimensionRule(inoutQN,ith,rule)
		class(U1fDimension),intent(inout) ::inoutQN
		integer,intent(in)::ith,rule
		if(ith.gt.inoutQN%outlenDimData())then
			if(inoutQN%outlenDimData().eq.0)then
				call writemess('There is no Data in the Symdimension')
				call error_stop()
			end if 
			call writemess('ERROR in SETTING Quantum Number')
			call writemess('The index is larger than the length of the dimension')
			call error_stop()
		end if
		if(.not.inoutQN%getArrowFlag())then
			call writemess("One should allocate fermi-arrow before setting the fermi-arrow")
			call error_stop()
		end if
		inoutQN%fermiArrow(ith)=rule
		return
	end subroutine
	subroutine setU1fDimensionRule_char(inoutQN,w,rule)
		class(U1fDimension),intent(inout) ::inoutQN
		character(len=*),intent(in)::w
		integer,intent(in)::rule
		integer::ith
		ith=inoutQN%Nameorder(w)
		if(ith.gt.inoutQN%outlenDimData())then
			if(inoutQN%outlenDimData().eq.0)then
				call writemess('There is no Data in the Symdimension')
				call error_stop()
			end if 
			call writemess('ERROR in SETTING Quantum Number')
			call writemess('The index is larger than the length of the dimension')
			call error_stop()
		end if
		if(.not.inoutQN%getArrowFlag())then
			call writemess("One should allocate fermi-arrow before setting the fermi-arrow")
			call error_stop()
		end if
		inoutQN%fermiArrow(ith)=rule
		return
	end subroutine
	subroutine setAllU1fDimensionRule(inoutQN,rule)
		class(U1fDimension),intent(inout) ::inoutQN
		integer,intent(in)::rule(:)
		integer::i,lenrule
		if(inoutQN%outlenDimData().eq.0)then
			call writemess('There is no Data in the Symdimension, initial first before setting rule')
			call error_stop()
		end if 
		if(size(rule).lt.inoutQN%outlenDimData())then
			call writemess('ERROR in SETTING Quantum rule')
			call writemess('The length of on put rule is smaller than the length of the dimension')
			call error_stop()
		end if
		inoutQN%ArrowFlag=.true.
		lenrule=size(rule)
		call allocateCheck(inoutQN%fermiArrow,lenrule)
		do i=1,lenrule
			call inoutQN%fermiArrow(i)=rule(i)
		end do
		return
	end subroutine







	
	subroutine Dprint(Dimen,uni)!overwrite the subroutine in dimension.f90
		class(U1fDimension),intent(in) ::Dimen
		integer,optional,intent(in)::uni
		integer::i
		character(len=characterlen)::w
		call Dimen%SymDimension%print(uni)
		if(present(uni))then
			if(Dimen%getArrowFlag())then
				write(uni,*)"Below are the fermi-arrow"
				w=''
				do i=1,Dimen%outlenDimData()-1
					w=w+(' '+Dimen%fermiArrow(i))+','
				end do
				w=w+(' '+Dimen%fermiArrow(i))
				write(uni,*)trim(w)
			else
				write(uni,*)"Do not Set the fermi-arrow"
			end if
		else
			if(Dimen%getArrowFlag())then
				call writemess("Below are the fermi-arrow")
				w=''
				do i=1,Dimen%outlenDimData()-1
					w=w+(' '+Dimen%fermiArrow(i))+','
				end do
				w=w+(' '+Dimen%fermiArrow(i))
				call writemess(w)
			else
				call writemess("Do not Set the fermi-arrow")
			end if
		end if
		return
	end subroutine
	
	subroutine Dprint2(Dimen,uni)!overwrite the subroutine in dimension.f90
		class(U1fDimension),intent(in) ::Dimen
		integer,optional,intent(in)::uni
		integer::i
		character(len=characterlen)::w
		call Dimen%Symdimension%info(uni)
		if(present(uni))then
			write(uni,*)"ArrowFlag ",Dimen%ArrowFlag
			if(.not.Dimen%ArrowFlag)return
			write(uni,*)"Below are the fermi-arrow"
			write(uni,*)Dimen%fermiArrow(1:Dimen%outlenDimData())
			write(uni,*)"-----"
		else
			call writemess("ArrowFlag:"+Dimen%ArrowFlag)
			if(.not.Dimen%ArrowFlag)return
			call writemess("Below are the fermi-arrow",-1)
			call writemess("Below are the fermi-arrow")
				w=''
				do i=1,Dimen%outlenDimData()-1
					w=w+(' '+Dimen%fermiArrow(i))+','
				end do
				w=w+(' '+Dimen%fermiArrow(i))
				call writemess(w)
			call writemess("------",-1)
		end if
		return
	end subroutine
	
	subroutine readdimension(dimen,uni)
		class(U1fdimension),intent(inout)::dimen
		integer,intent(in)::uni
		character*50::notused
		integer::i,j,k,total,rule
		integer,allocatable::degeneracy(:)
		real*4,allocatable::QuanNum(:)
		real*4::maxQN
		call dimen%Symdimension%read(uni)
		read(uni,*)notused,Dimen%ArrowFlag
		if(.not.Dimen%ArrowFlag)return
		read(uni,*)notused
		total=Dimen%outlenDimData()
		call allocateCheck(Dimen%fermiArrow,total)
		read(uni,*)(Dimen%fermiArrow(i),i=1,total)
		read(uni,*)notused
		return
	end subroutine
	


!*******  function or subroutine for name   **************	
!	return the inde  dimension	,outpout in a type(dimenison)
	type(U1fdimension) function  getSubDim2(Dimen,inde)
		type(U1fdimension),intent(in) :: Dimen
		integer,intent(in) :: inde
		integer::i,j,k,Dlen,boundary(2)
		getSubDim2%Symdimension=Dimen%Symdimension.sub.inde
		if(.not.Dimen%getArrowFlag())return
		call Dimen%getSubDimboundary(inde,boundary)
		Dlen=boundary(2)-boundary(1)
		if(Dlen.ne.getSubDim2%dimension%outlenDimData())then
			call writemess("ERROR in getSubDim2,U1fDimension.f90")
			call error_stop()
		end if
		call allocateCheck(getSubDim2%fermiArrow,Dlen)
		getSubDim2%ArrowFlag=.true.
		k=1
		do i=boundary(1)+1,boundary(2)
			getSubDim2%fermiArrow(k)=Dimen%fermiArrow(i)
			k=k+1
		end do
		return
	end function
	
	
	type(U1fdimension) function  getSubDim3(Dimen,inde) 
		type(U1fdimension),intent(in) :: Dimen
		integer,intent(in) :: inde(2)
		integer::ith(2),Dlen
		getSubDim3%Symdimension=Dimen%Symdimension.sub.inde
		if(.not.Dimen%getArrowFlag())return
		getSubDim3%ArrowFlag=.true.
		call getSubDim3_index_routine(ith,Dimen%Dimension,inde) 
		Dlen=ith(2)-ith(1)+1
		allocate(getSubDim3%fermiArrow(Dlen))
		getSubDim3%fermiArrow=Dimen%fermiArrow(ith(1):ith(2))
		return
	end function
	type(U1fDimension) function  getSubDim4(Dimen)
		type(U1fDimension),intent(in) :: Dimen
		getSubDim4=Dimen
		return
	end function
	type(U1fDimension) function  getSubDim2_name(Dimen,w)
		type(U1fDimension),intent(in) :: Dimen
		CHARACTER(len=*),intent(in)::w
		integer::inde
		inde=Dimen%Nameorder(w)
		getSubDim2_name=getSubDim2(Dimen,inde)
		return
	end function
		
!**************** fuse   ****************
!		combine two index of the Tensor,which is con_index and con_index+1
	type(U1fDimension) function fuseDimension_val(Dimen,inde)result(fuseDim)!overwrite
		integer,intent(in) :: inde
		type(U1fDimension),intent(in) :: Dimen
		fuseDim%SymDimension=dimen%SymDimension.fuse.(inde)
		fuseDim%ArrowFlag=dimen%ArrowFlag
		if(fuseDim%ArrowFlag)then
			allocate(fuseDim%fermiArrow(fuseDim%outlenDimData()))
			fuseDim%fermiArrow=dimen%fermiArrow
		end if
		return
	end function
	type(U1fDimension) function fuseDimension_vec(dimen,vector)result(fuseDim)
		integer,intent(in) ::vector(2)
		type(U1fDimension),intent(in) :: dimen
		fuseDim%SymDimension=dimen%SymDimension.fuse.vector
		fuseDim%ArrowFlag=dimen%ArrowFlag
		if(fuseDim%ArrowFlag)then
			allocate(fuseDim%fermiArrow(fuseDim%outlenDimData()))
			fuseDim%fermiArrow=dimen%fermiArrow
		end if
		return
	end function		
	type(U1fDimension) function splitDimension2(dimen,vector)result(splitDimension)
		type(U1fDimension),intent(in) :: dimen
		integer,intent(in) ::vector(2)
		splitDimension%SymDimension=dimen%SymDimension.split.vector
		splitDimension%ArrowFlag=dimen%ArrowFlag
		if(splitDimension%ArrowFlag)then
			allocate(splitDimension%fermiArrow(splitDimension%outlenDimData()))
			splitDimension%fermiArrow=dimen%fermiArrow
		end if
		return
	end function
			
	type(U1fDimension) function splitDimension3(dimen,de_index)result(splitDimension)
		type(U1fDimension),intent(in) :: dimen
		integer,intent(in) :: de_index
		splitDimension%SymDimension=dimen%SymDimension.split.de_index
		splitDimension%ArrowFlag=dimen%ArrowFlag
		if(splitDimension%ArrowFlag)then
			allocate(splitDimension%fermiArrow(splitDimension%outlenDimData()))
			splitDimension%fermiArrow=dimen%fermiArrow
		end if
		return
	end function
			
	type(U1fDimension) function splitDimensionAll(Dimen)result(splitDimension)
		type(U1fDimension),intent(in) :: Dimen
		splitDimension%SymDimension=.split.dimen%SymDimension
		splitDimension%ArrowFlag=dimen%ArrowFlag
		if(splitDimension%ArrowFlag)then
			allocate(splitDimension%fermiArrow(splitDimension%outlenDimData()))
			splitDimension%fermiArrow=dimen%fermiArrow
		end if
		return
	end function	
	
	subroutine getSubDimArrow(Dimen,ith,Arrow,outlen)
		type(U1fDimension),intent(in)::Dimen
		integer,intent(in)::ith
		integer,intent(in)::outlen
		integer,intent(inout)::Arrow(outlen)
		integer::boundary(2),i,Dlen,k
		call Dimen%getSubDimboundary(ith,boundary)
		Dlen=boundary(2)-boundary(1)
		if(Dlen.ne.outlen)then
			call writemess("ERROR in getSubDimArrow,U1fDimension.f90")
			call error_stop()
		end if
		k=1
		do i=boundary(1)+1,boundary(2)
			Arrow(k)=Dimen%Arrow(i)
			k=k+1
		end do
		return
	end subroutine
	
	type(U1fDimension) function Dimpermute(dimen,v)
		class(U1fDimension),intent(in) :: Dimen
		integer,intent(in):: v(:)
		integer::i,datalen,subDlen,k
		integer,allocatable :: fermiArrow(:)
		Dimpermute%SymDimension=dimen%SymDimension%Dimpermute(v)
		if(.not.dimen%getArrowFlag()) return
		datalen=dimen%outlenDimData()
		allocate(Dimpermute%fermiArrow(datalen))
		Dimpermute%ArrowFlag=dimen%ArrowFlag
		
		if(dimen%out_sample_flag())then
			do i=1,datalen
				Dimpermute%fermiArrow(i)=dimen%fermiArrow(v(i))
			end do
			return
		end if
		k=1
		do i=1,Dimen%getRank()
			subDlen=dimen%getsubDimlen(v(i))
			call allocateCheck(fermiArrow,subDlen)
			call getSubDimArrow(Dimen,v(i),fermiArrow,subDlen)
			Dimpermute%fermiArrow(k:subDlen+k-1)=fermiArrow(1:subDlen)
			k=k+subDlen
		end do
		return
	end function	

	! order is [ith,1,2,3...]

	subroutine  U1fDimpermute_forwards(outdim,dimen,ith)
		type(U1fDimension),intent(inout)::outdim
		type(U1fDimension),intent(in) :: Dimen
		integer,intent(in):: ith
		integer,allocatable::order(:)
		integer::length,i
		length=Dimen%getRank()
		allocate(order(length))
		order(1)=ith
		do i=2,ith
			order(i)=i-1
		end do
		do i=ith+1,length
			order(i)=i
		end do
		outdim=Dimpermute(Dimen,order)
		return
	end subroutine

! order is [1,2,3...,n,ith]

	subroutine U1fDimpermute_backwards(outdim,dimen,ith)
		type(U1fDimension),intent(inout)::outdim
		type(U1fDimension),intent(in) :: Dimen
		integer,intent(in):: ith
		integer,allocatable::order(:)
		integer::length,i
		length=Dimen%getRank()
		allocate(order(length))
		order(length)=ith
		do i=1,ith-1
			order(i)=i
		end do
		do i=ith,length-1
			order(i)=i+1
		end do
		outdim=Dimpermute(Dimen,order)
		return
	end subroutine


! oeder is [2,3,4,...,ith,1,ith+1,...]

	subroutine U1fDimpermute_forwards_index(outdim,dimen,ith)
		type(U1fDimension),intent(inout)::outdim
		type(U1fDimension),intent(in) :: Dimen
		integer,intent(in):: ith
		integer,allocatable::order(:)
		integer::length,i
		length=Dimen%getRank()
		allocate(order(length))
		order(ith)=1
		do i=1,ith-1
			order(i)=i+1
		end do
		do i=ith+1,length
			order(i)=i
		end do
		outdim=Dimpermute(Dimen,order)
		return
	end subroutine

	! oeder is [1,2,3,4,...,ith,n,ith+1,...,n-1]

	subroutine U1fDimpermute_backwards_index(outdim,dimen,ith)
		type(U1fDimension),intent(inout)::outdim
		type(U1fDimension),intent(in) :: Dimen
		integer,intent(in):: ith
		integer,allocatable::order(:)
		integer::length,i
		length=Dimen%getRank()
		allocate(order(length))
		order(ith)=length
		do i=1,ith-1
			order(i)=i
		end do
		do i=ith+1,length
			order(i)=i-1
		end do
		outdim=Dimpermute(Dimen,order)
		return
	end subroutine

	
	
	type(U1fDimension) function  DimAdd(Dimen,Dimen2)
		type(U1fDimension),intent(in) :: Dimen,Dimen2
		integer::l1,l2
		DimAdd%SymDimension=Dimen%SymDimension+Dimen2%SymDimension
		if(Dimen%getArrowFlag() .and. Dimen2%getArrowFlag())then
			l1=Dimen%outlenDimData()
			l2=Dimen2%outlenDimData()
			allocate(DimAdd%fermiArrow(l1+l2))
			DimAdd%fermiArrow(1:l1)=Dimen%fermiArrow(1:l1)
			DimAdd%fermiArrow(l1+1:)=Dimen2%fermiArrow(1:l2)
			DimAdd%ArrowFlag=.true.
			return
		end if
		if( (.not.Dimen%getArrowFlag()) .and. (.not.Dimen2%getArrowFlag()))then
			DimAdd%ArrowFlag=.false.
			return
		end if
		call writemess('ERRO in U1fDimension + U1fDimension')
		call writemess('There is no fermi-arrow in one or both Symdimension')
		call error_stop()
	end function
	type(U1fDimension) function  DimAddQuanNum(Dimen,Quan)result(DimAdd)
		type(U1fDimension),intent(in) :: Dimen
		type(QuanNum),intent(in)::Quan
		integer::l1
		DimAdd%Symdimension=Dimen%Symdimension+Quan
		if(Dimen%getArrowFlag())then
			l1=Dimen%outlenDimData()
			allocate(DimAdd%fermiArrow(l1+1))
			DimAdd%fermiArrow(1:l1)=Dimen%fermiArrow(1:l1)
			DimAdd%fermiArrow(l1+1)=0
			DimAdd%ArrowFlag=.true.
			return
		end if
		return
	end function
	type(U1fDimension) function  QuanNumAddDim(Quan,Dimen)result(DimAdd)
		type(QuanNum),intent(in)::Quan
		type(U1fDimension),intent(in) :: Dimen
		integer::l1,l2
		DimAdd%Symdimension=Quan+ Dimen%Symdimension
		if(Dimen%getArrowFlag())then
			l2=Dimen%outlenDimData()
			allocate(DimAdd%fermiArrow(l2+1))
			DimAdd%fermiArrow(1)=0
			DimAdd%fermiArrow(2:l2+1)=Dimen%fermiArrow(1:l2)
			DimAdd%ArrowFlag=.true.
			return
		end if
		return
	end function	
	
	

	
	

!**********************************************************************
!**********************************************************************
!	the code below is for MPI
!**********************************************************************
		
	
	subroutine MPI_send_U1fDimension(Dimen1,Dimen2,ID1,ID2,ierr,MPIcommon)
		type(U1fDimension),intent(in)::Dimen1
		type(U1fDimension),intent(inout)::Dimen2
		integer,intent(in)::ID1,ID2
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,istatus(MPI_STATUS_SIZE),mpi_comm,lendata,i
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		if(present(MPIcommon))then
			if((ID1.ge.proNum).or.(ID2.ge.proNum))return
		end if
		
		if(ID1.eq.ID2) return !The same cpu, do nothing
		
		if((proID.ne.ID1).and.(proID.ne.ID2)) return!The proID do not send or recv, return
		
		call MPI_send_SymDimension(Dimen1%SymDimension,Dimen2%SymDimension,ID1,ID2,ierr,MPIcommon)
!**************************flag************************************************				
		if(proID.eq.ID1) then
			call mpi_send(Dimen1%ArrowFlag,1,MPI_logical,ID2,tag,MPI_Comm,ierr)
			if(.not.Dimen1%ArrowFlag)then
				return
			else
				lendata=Dimen1%outlenDimData()
			end if

		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%ArrowFlag,1,MPI_logical,ID1,tag,MPI_Comm,istatus,ierr)
			if(.not.Dimen2%ArrowFlag)then
				return
			else
				lendata=Dimen2%outlenDimData()
				allocate(Dimen2%fermiArrow(lendata))
			end if
		end if	

		if(proID.eq.ID1) then	
			call mpi_send(Dimen1%fermiArrow,lendata,MPI_integer,ID2,tag,MPI_Comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Dimen2%fermiArrow,lendata,MPI_integer,ID1,tag,MPI_Comm,istatus,ierr)
		end if
		return
	end subroutine
	
	subroutine MPI_BCAST_U1fDimension(Dimen,ID,ierr,MPIcommon)
		type(U1fDimension),intent(inout)::Dimen
		integer,intent(in)::ID
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,istatus(MPI_STATUS_SIZE),mpi_comm,lendata,i
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
		
		call MPI_BCAST_SymDimension(Dimen%SymDimension,ID,ierr,MPIcommon)
		
		call MPI_BCAST(Dimen%ArrowFlag,1,MPI_logical,ID,mpi_comm,ierr)		
		if(Dimen%ArrowFlag)then
			lendata=Dimen%outlenDimData()
			if(proId.ne.ID)then
				call allocateCheck(Dimen%fermiArrow,lendata)
			end if
		else
			return
		end if

		call MPI_BCAST(Dimen%fermiArrow,lendata,MPI_integer,ID,mpi_comm,ierr)
		return
	end subroutine
		

end module










