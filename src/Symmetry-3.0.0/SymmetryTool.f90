
module Symmetry_tool
	use SymDimension_typede
	use usefull_function
	use Tensor_complex
	implicit none
	private
	public::SymTool
	type SymTool
		integer::SymmetryType=0! 0 Parity
		                       ! 1 U1
	contains
		procedure,public::outType=>outsymmetryType
		procedure,nopass::PRule1,PRule2
		generic,public::PRule=>PRule1,PRule2
		generic,public::U1Rule=>U1Rule1,U1RuleArray	
		procedure,nopass::U1Rule1
		procedure,nopass::U1RuleArray
		procedure,public::fuseOrder=>ParityU1FuseOrder
		procedure,public::Rule
		procedure,public::NonZeroElement
		generic,public::QN=>U1_1,U1_2,U1_3,U1_4,Parity
		procedure::U1_1,U1_2,U1_3,U1_4,Parity
		generic,public::setSymmetry=>setSymmetry1,setSymmetry2
		procedure::setSymmetry1,setSymmetry2
	end type SymTool

!
!order of codes for do
! outQ=Q1+Q2 
!  first order Q1 and then Q2 at lase outQ
!Parity example:
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
contains
	integer function outsymmetryType(Tool)
		class(SymTool),intent(in)::Tool
		outsymmetryType=Tool%SymmetryType
		return
	end function
	logical function Rule(PTool,i,j,k,QN1,QN2,QN3)
		class(SymTool),intent(in)::PTool
		integer,intent(in)::i,j,k
		type(QuanNum),intent(in)::QN1,QN2,QN3
		real*4::Q1,Q2,Q3
		integer::Urule1,Urule2,Urule3
		Q1=QN1%getQN(i)
		Q2=QN2%getQN(j)
		Q3=QN3%getQN(k)
		select case (PTool%SymmetryType)
			case (0)
				Rule=PRule1(Q1,Q2,Q3)
			case (1)
				Urule1=QN1%getRule()
				Urule2=QN2%getRule()
				Urule3=QN3%getRule()
				Rule=U1Rule1(Q1,Q2,Q3,Urule1,Urule2,Urule3)
		end select
	end function
	subroutine setSymmetry1(PTool,cha)
		class(SymTool),intent(inout)::PTool
		character(len=*),intent(in)::cha
		if(cha.equ.'U1') then
			PTool%SymmetryType=1
			return
		end if
		if(cha.equ.'Parity') then
			PTool%SymmetryType=0
			return
		end if
		call writemess("No such symmetry type")
		call error_stop()
	end subroutine
	subroutine setSymmetry2(PTool,cha)
		class(SymTool),intent(inout)::PTool
		integer,intent(in)::cha
		select case (cha)
			case (0)
				PTool%SymmetryType=0
			case (1)
				PTool%SymmetryType=1
			case default
				call writemess("No such symmetry type")
			call error_stop()
		end select
		return
	end subroutine
	
	
	
	
	
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
	
	type(QuanNum) function Parity(PTool,deg)
		class(SymTool),intent(in)::PTool
		integer,intent(in)::deg(:)
		if(PTool%SymmetryType.ne.0)then
			call writemess("ERROR Symmetry")
			call writemess("Symmetry is")
			select case (PTool%SymmetryType)
				case (0)
					call writemess("Parity")
				case (1)
					call writemess("U(1)")
			end select
			call error_stop()
		end if
		call Parity%setRule(1)
		call Parity%setMaxQN(1.)
		call Parity%setQN((/-1.,1./))
		call Parity%setDeg(Deg)
		return
	end function
	
	type(QuanNum) function U1_1(PTool,maxQN_)result(outQ)
		class(SymTool),intent(in)::PTool
		real*4,intent(in)::maxQN_
		integer::lenQN,i
		real*4,allocatable::outQuanNum(:)
		real*4::temp,maxQN
		logical::goon
		type(QuanNum)::Q1,Q2
		if(PTool%SymmetryType.ne.1)then
			call writemess("ERROR Symmetry")
			call writemess("Symmetry is")
			select case (PTool%SymmetryType)
				case (0)
					call writemess("Parity")
				case (1)
					call writemess("U(1)")
			end select
			call error_stop()
		end if
		maxQN=abs(maxQN_)
		if(abs(maxQN-0.5).le.1e-6)then
			call outQ%setRule(1)
			call outQ%SetQN((/-0.5,0.5/),(/1,1/))
			return
		end if
		call Q1%setRule(1)
		call Q1%SetQN((/-0.5,0.5/),(/1,1/))
		call Q2%setRule(1)
		call Q2%SetQN((/-0.5,0.5/),(/1,1/))
		goon =.true.
		do while (goon)
			if(maxQN.gt.maxQN_)then
				call writemess("ERROR Quantum Number input")
				call error_stop()
			end if
			maxQN=Q1%getMaxQN()+Q2%getMaxQN()
			call ADDQuan(outQ,Q1,Q2,maxQN)
			Q1=outQ
			if(abs(maxQN-maxQN_).le.1e-6) goon=.false.
		end do
		call outQ%setRule(1)
		return
	end function
	type(QuanNum) function U1_2(PTool,maxQN_,rule)result(outQ)
		class(SymTool),intent(in)::PTool
		real*4,intent(in)::maxQN_
		integer,intent(in)::rule
		integer::lenQN,i
		real*4,allocatable::outQuanNum(:)
		real*4::temp,maxQN
		logical::goon
		type(QuanNum)::Q1,Q2
		if(PTool%SymmetryType.ne.1)then
			call writemess("ERROR Symmetry")
			call writemess("Symmetry is")
			select case (PTool%SymmetryType)
				case (0)
					call writemess("Parity")
				case (1)
					call writemess("U(1)")
			end select
			call error_stop()
		end if
		maxQN=abs(maxQN_)
		if(abs(maxQN-0.5).le.1e-6)then
			call outQ%setRule(1)
			call outQ%SetQN((/-0.5,0.5/),(/1,1/))
			return
		end if
		call Q1%setRule(1)
		call Q1%SetQN((/-0.5,0.5/),(/1,1/))
		call Q2%setRule(1)
		call Q2%SetQN((/-0.5,0.5/),(/1,1/))
		goon =.true.
		do while (goon)
			if(maxQN.gt.maxQN_)then
				call writemess("ERROR Quantum Number input")
				call error_stop()
			end if
			maxQN=Q1%getMaxQN()+Q2%getMaxQN()
			call ADDQuan(outQ,Q1,Q2,maxQN)
			Q1=outQ
			if(abs(maxQN-maxQN_).le.1e-6) goon=.false.
		end do
		call outQ%setRule(rule)
		return
	end function
	
	type(QuanNum) function U1_1_old(PTool,maxQN_)result(outQ)
		class(SymTool),intent(in)::PTool
		real*4,intent(in)::maxQN_
		integer::lenQN,i
		real*4,allocatable::outQuanNum(:)
		real*4::temp,maxQN
		if(PTool%SymmetryType.ne.1)then
			call writemess("ERROR Symmetry")
			call writemess("Symmetry is")
			select case (PTool%SymmetryType)
				case (0)
					call writemess("Parity")
				case (1)
					call writemess("U(1)")
			end select
			call error_stop()
		end if
		maxQN=abs(maxQN_)
		call outQ%setRule(1)
		call outQ%setMaxQN(maxQN)
		lenQN=2*maxQN+1
		allocate(outQuanNum(lenQN))
		temp=-maxQN
		do i=1,lenQN
			outQuanNum(i)=temp
			temp=temp+1
		end do
		call outQ%setQN(outQuanNum)
		return
	end function
	type(QuanNum) function U1_2_old(PTool,maxQN_,rule)result(outQ)
		class(SymTool),intent(in)::PTool
		real*4,intent(in)::maxQN_
		integer,intent(in)::rule
		integer::lenQN,i
		real*4,allocatable::outQuanNum(:)
		real*4::temp,maxQN
		if(PTool%SymmetryType.ne.1)then
			call writemess("ERROR Symmetry")
			call writemess("Symmetry is")
			select case (PTool%SymmetryType)
				case (0)
					call writemess("Parity")
				case (1)
					call writemess("U(1)")
			end select
			call error_stop()
		end if
		maxQN=abs(maxQN_)
		call outQ%setRule(rule)
		call outQ%setMaxQN(maxQN)
		lenQN=2*maxQN+1
		allocate(outQuanNum(lenQN))
		temp=-maxQN
		do i=1,lenQN
			outQuanNum(i)=temp
			temp=temp+1
		end do
		return
	end function
!We suppose S is made up of many si=0.5 spin
!
	subroutine ADDQuan(outQ,Q1,Q2,maxQN)
		type(QuanNum),intent(inout)::outQ
		type(QuanNum),intent(in)::Q1,Q2
		real*4,intent(in)::maxQN
		integer::i,j,k,lenQN,rule1,rule2,deg1,deg2,degstart,rule
		real*4::ouQN,QN1,QN2,temp
		real*4,allocatable::outQuanNum(:)
		integer,allocatable::outDegen(:)
		logical::goon
		if(Q1%getmaxQN()+Q2%getmaxQN().ne.maxQN)then
			if(abs(Q1%getmaxQN()-Q2%getmaxQN()).ne.maxQN)then
				call writemess("Can not fuse the quantum number,error in SymmetryTool.f90")
				write(*,*)Q1%getmaxQN(),Q2%getmaxQN(),maxQN
				call error_stop()
			end if
		end if
		call outQ%setMaxQN(maxQN)
		lenQN=2*maxQN+1
		allocate(outQuanNum(lenQN))
		allocate(outDegen(lenQN))
		temp=-maxQN
		do i=1,lenQN
			outQuanNum(i)=temp
			temp=temp+1
		end do
		call outQ%setQN(outQuanNum)
		outDegen=0
		rule1=1
		rule2=1
		rule=-1
		do i=1,lenQN
			do j=1,Q2%getQNlength()
				do k=1,Q1%getQNlength()
					ouQN=outQ%getQN(i)
					QN1=Q1%getQN(k)
					QN2=Q2%getQN(j)
					goon=U1Rule1(ouQN,QN1,QN2,rule,rule1,rule2)
					if(goon)then
						deg1=Q1%getDeg(k)
						deg2=Q2%getDeg(j)
						outDegen(i)=outDegen(i)+deg1*deg2
					end if
				end do
			end do
		end do
		call outQ%setDeg(outDegen)
	end subroutine
	type(QuanNum) function U1_3(PTool,maxQN_,outDegen)result(outQ)
		class(SymTool),intent(in)::PTool
		real*4,intent(in)::maxQN_
		integer,intent(in)::outDegen(:)
		integer::lenQN,i
		real*4,allocatable::outQuanNum(:)
		real*4::temp,maxQN
		if(PTool%SymmetryType.ne.1)then
			call writemess("ERROR Symmetry")
			call writemess("Symmetry is")
			select case (PTool%SymmetryType)
				case (0)
					call writemess("Parity")
				case (1)
					call writemess("U(1)")
			end select
			call error_stop()
		end if
		maxQN=abs(maxQN_)
		call outQ%setRule(1)
		call outQ%setMaxQN(maxQN)
		lenQN=2*maxQN+1
		allocate(outQuanNum(lenQN))
		temp=-maxQN
		do i=1,lenQN
			outQuanNum(i)=temp
			temp=temp+1
		end do
		call outQ%setQN(outQuanNum)
		call outQ%setDeg(outDegen)
		return
	end function
	type(QuanNum) function U1_4(PTool,maxQN_,outDegen,rule)result(outQ)
		class(SymTool),intent(in)::PTool
		real*4,intent(in)::maxQN_
		integer,intent(in)::rule,outDegen(:)
		integer::lenQN,i
		real*4,allocatable::outQuanNum(:)
		real*4::temp,maxQN
		if(PTool%SymmetryType.ne.1)then
			call writemess("ERROR Symmetry")
			call writemess("Symmetry is")
			select case (PTool%SymmetryType)
				case (0)
					call writemess("Parity")
				case (1)
					call writemess("U(1)")
			end select
			call error_stop()
		end if
		maxQN=abs(maxQN_)
		call outQ%setRule(rule)
		call outQ%setMaxQN(maxQN)
		lenQN=2*maxQN+1
		allocate(outQuanNum(lenQN))
		temp=-maxQN
		do i=1,lenQN
			outQuanNum(i)=temp
			temp=temp+1
		end do
		call outQ%setQN(outQuanNum)
		call outQ%setDeg(outDegen)
		return
	end function

!outQ=Q1+Q2
	type(QuanNum) function ParityU1FuseOrder(Tool,order,Q1,Q2,maxQN,rule)
		class(SymTool),intent(in)::Tool
		type(QuanNum),intent(in)::Q1,Q2
		type(Tensor),intent(inout)::order
		real*4,optional,intent(in)::maxQN
		integer,optional,intent(in)::rule
		select case (Tool%SymmetryType)
			case (0)
				call ParityFuseOrder(order,ParityU1FuseOrder,Q1,Q2)
			case (1)
				call U1FuseOrder(order,ParityU1FuseOrder,Q1,Q2,maxQN,rule)
		end select
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
! Q3=Q1+Q2
!                    Q3
!             -                  +
!        (+ -)  (- +)         (- -)  (+ +)
!         3*4    2*5           2*4   3*5 
!
! 
	subroutine ParityFuseOrder(order,outQ,Q1,Q2)
		type(QuanNum),intent(inout)::outQ
		type(QuanNum),intent(in)::Q1,Q2
		type(Tensor),intent(inout)::order
		integer::i,j,k,deg(2),deg1,deg2,degstart
		real*4::ouQN,QN1,QN2
		logical::goon
		call outQ%setRule(1)
		call outQ%setMaxQN(1.)
		call outQ%setQN((/-1.,1./))
		deg=0
		call order%empty()
		do i=1,2
			degstart=0
			do j=1,2
				do k=1,2
					ouQN=outQ%getQN(i)
					QN1=Q1%getQN(k)
					QN2=Q2%getQN(j)
					goon=PRule1(ouQN,QN1,QN2)
					if(goon)then
						deg1=Q1%getDeg(k)
						deg2=Q2%getDeg(j)
						deg(i)=deg(i)+deg1*deg2
						if(order%getFlag())then
							call order%addrow(Tensor( (/i,k,j,deg1,deg2,degstart+1,deg(i)/)))
							degstart=deg(i)
						else
							order=(/i,k,j,deg1,deg2,1,deg(i)/)
							degstart=deg(i)
						end if
					end if
				end do
			end do
		end do
		call outQ%setDeg(Deg)
		return
	end subroutine
	
	
	logical function U1Rule1(Q1,Q2,Q3,rule1,rule2,rule3)
		real*4,intent(in)::Q1,Q2,Q3
		integer,intent(in)::rule1,rule2,rule3
		real*4::temp
		temp=rule1*Q1+rule2*Q2+rule3*Q3
		U1Rule1=abs(temp).le.1e-7
		return
	end function
	logical function U1RuleArray(Q,rule)
		real*4,intent(in)::Q(:)
		integer,intent(in)::rule(:)
		real*4::temp
		integer::I
		temp=0.
		do i=1,size(Q)
			temp=temp+Q(i)*rule(i)
		end do
		U1RuleArray=abs(temp).le.1e-7
		return
	end function
	

!
!order of codes for do
! outQ=Q1+Q2 
!  first order Q1 and then Q2 at lase outQ


	subroutine U1FuseOrder(order,outQ,Q1,Q2,maxQN,rule_)
		type(QuanNum),intent(inout)::outQ
		type(QuanNum),intent(in)::Q1,Q2
		type(Tensor),intent(inout)::order
		integer,intent(in)::rule_
		real*4,intent(in)::maxQN
		integer::i,j,k,lenQN,rule1,rule2,deg1,deg2,degstart,rule
		real*4::ouQN,QN1,QN2,temp
		real*4,allocatable::outQuanNum(:)
		integer,allocatable::outDegen(:)
		logical::goon
		if(Q1%getmaxQN()+Q2%getmaxQN().ne.maxQN)then
			if(abs(Q1%getmaxQN()-Q2%getmaxQN()).ne.maxQN)then
				call writemess("Can not fuse the quantum number,error in SymmetryTool.f90")
				write(*,*)Q1%getmaxQN(),Q2%getmaxQN(),maxQN
				call error_stop()
			end if
		end if
		rule=-1*rule_
		call outQ%setMaxQN(maxQN)
		lenQN=2*maxQN+1
		allocate(outQuanNum(lenQN))
		allocate(outDegen(lenQN))
		temp=-maxQN
		do i=1,lenQN
			outQuanNum(i)=temp
			temp=temp+1
		end do
		call outQ%setQN(outQuanNum)
		outDegen=0
		call order%empty()
		do i=1,lenQN
			degstart=0
			do j=1,Q2%getQNlength()
				do k=1,Q1%getQNlength()
					ouQN=outQ%getQN(i)
					QN1=Q1%getQN(k)
					QN2=Q2%getQN(j)
					rule1=Q1%getRule()
					rule2=Q2%getRule()
					goon=U1Rule1(ouQN,QN1,QN2,rule,rule1,rule2)
					if(goon)then
						deg1=Q1%getDeg(k)
						deg2=Q2%getDeg(j)
						outDegen(i)=outDegen(i)+deg1*deg2
						if(order%getFlag())then
							call order%addrow(Tensor( (/i,k,j,deg1,deg2,degstart+1,outDegen(i)/)) )
							degstart=outDegen(i)
						else
							order=(/i,k,j,deg1,deg2,1,outDegen(i)/)
							degstart=outDegen(i)
						end if
					end if
				end do
			end do
		end do
		call outQ%setDeg(outDegen)
		call outQ%setRule(rule_)
	end subroutine
	
	
	type(Tensor) function NonZeroElement(Tool,dimen)
		class(SymTool),intent(in)::Tool
		class(SymDimension),intent(in)::dimen
		select case (Tool%SymmetryType)
			case (0)
				call NonZeroElementParity(NonZeroElement,dimen)
			case (1)
				call NonZeroElementU1(NonZeroElement,dimen)
		end select
		return
	end function
	
	subroutine NonZeroElementU1(Res,dimen)
		type(Tensor),intent(inout)::Res
		class(SymDimension),intent(in)::dimen
		integer,allocatable::indices(:),maxdim(:),mindim(:),rule(:)
		integer::Rank
		logical::goon,rulefit
		real*4,allocatable::QN(:)
		rank=dimen%getRank()
		allocate(indices(rank))
		allocate(maxdim(rank))
		allocate(mindim(rank))
		allocate(QN(rank))
		allocate(rule(rank))
		mindim=1
		indices=mindim
		goon=.true.
		rule=dimen%getRule()
		maxdim=dimen%dim()
		call Res%empty()
		do while(goon)
			QN=dimen%QNDim(indices)
			rulefit=U1ruleArray(QN,rule)
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

