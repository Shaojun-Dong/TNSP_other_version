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

module U1Tool
	use QuantumNumber_Type
	use SymDimension_typede
	use SymTensor_type
	use Tools
	use mpi
	use Tensor_type
	implicit none
	private
	public::NonZeroElement,rule,fuseOrder,quantumnumber,checkSymmetryRule,reverseSymmetryRule
	interface rule
		module procedure U1Rule1
		module procedure U1Rule2
		module procedure U1Rule3
	end interface
	interface fuseOrder
		module procedure QuanFuseOrderFunction
	end interface
	interface quantumnumber
		module procedure U1QuantumNumber0
		module procedure U1QuantumNumber0_
		module procedure U1QuantumNumber0_2
		module procedure U1QuantumNumber1
		module procedure U1QuantumNumber2
		module procedure U1QuantumNumber3
		module procedure U1QuantumNumber1_
		module procedure U1QuantumNumber2_
		module procedure U1QuantumNumber3_
	end interface

	public::ifParity
	interface ifParity
		module procedure ifParityIndex1
		module procedure ifParityIndex2
		module procedure ifParityIndex3
	end interface

	public::QuanNumParity
	interface QuanNumParity
		module procedure Parity1
		module procedure Parity2
	end interface
	
	
	
contains
	
	
	logical function U1Rule1(Q1,Q2,Q3,R1,R2,R3)
		real*4,intent(in)::Q1,Q2,Q3
		integer,intent(in)::R1,R2,R3
		U1Rule1=((Q1*R1)+(Q2*R2)+(Q3*R3)).equ.0e0
		return
	end function
	logical function U1Rule2(Q,R)
		real*4,intent(in)::Q(:)
		integer,intent(in)::R(:)
		real*4::temp
		integer::i
		temp=Q(1)*R(1)
		do i=2,size(Q)
			temp=temp+(Q(i)*R(i))
		end do
		U1Rule2=temp.equ.0e0
		return
	end function
	
	logical function U1Rule3(dimen,indices)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::indices(:)
		real*4::Q(size(indices))
		integer::R(size(indices))
		integer::i
		real*4::temp
		do i=1,size(indices)
			Q(i)=dimen%getQN(i,indices(i))
			R(i)=dimen%getRule(i)
		end do
		temp=Q(1)*R(1)
		do i=2,size(Q)
			temp=temp+(Q(i)*R(i))
		end do
		U1Rule3=temp.equ.0e0
		return
	end function


	subroutine QNcheck(QN)
		real*4,intent(in)::QN(:)
		integer::lenQN,lenQN2
		lenQN=size(QN)
		lenQN2=anint(QN(lenQN)-QN(1))+1
		if(lenQN.ne.lenQN2)then
			call writemess('ERROR in the quantum number of U(1) ')
			call writemess('the input quantum number shoulb be look like:')
			call writemess(' i,i+1,i+3,i+4,....')
			call error_stop
		end if
		return
	end subroutine

	type(QuanNum) function U1QuantumNumber0(QN)
		real*4,intent(in)::QN(:)
		call QNcheck(QN)
		call U1QuantumNumber0%setQN(QN)
		return
	end function

	type(QuanNum) function U1QuantumNumber0_(QN1,QN2)
		real*4,intent(in)::QN1,QN2
		real*4,allocatable::QN(:)
		integer::lenQN,i
		if(QN1.gt.QN2)then
			call writemess('ERROR in input QN1,QN2, QN1 > QN2')
			call error_stop
		end if
		lenQN=int(anint(QN2))-int(anint(QN1))+1
		allocate(QN(lenQN))
		do i=1,lenQN
			QN(i)=QN1+(i-1)
		end do		
		call U1QuantumNumber0_%setQN(QN)
		return
	end function
	type(QuanNum) function U1QuantumNumber0_2(QN1)
		real*4,intent(in)::QN1
		call U1QuantumNumber0_2%setQN([QN1])
		return
	end function
	
	type(QuanNum) function U1QuantumNumber1(QN,deg)
		integer,intent(in)::deg(:)
		real*4,intent(in)::QN(:)
		call QNcheck(QN)
		call U1QuantumNumber1%setQN(QN)
		call U1QuantumNumber1%setDeg(Deg)
		return
	end function
	type(QuanNum) function U1QuantumNumber1_(QN1,QN2,deg)result(U1QuantumNumber1)
		integer,intent(in)::deg(:)
		real*4,intent(in)::QN1,QN2
		real*4,allocatable::QN(:)
		integer::lenQN,i
		if(QN1.gt.QN2)then
			call writemess('ERROR in input QN1,QN2, QN1 > QN2')
			call error_stop
		end if
		lenQN=int(anint(QN2))-int(anint(QN1))+1
		if(lenQN.ne.size(deg))then
			call writemess('ERROR in input deg, size(deg) .ne. lenQN')
			call error_stop
		end if
		allocate(QN(lenQN))
		do i=1,lenQN
			QN(i)=QN1+(i-1)
		end do
		call U1QuantumNumber1%setQN(QN)
		call U1QuantumNumber1%setDeg(Deg)
		return
	end function
	type(QuanNum) function U1QuantumNumber2(QN,deg,rule)
		integer,intent(in)::deg(:),rule
		real*4,intent(in)::QN(:)
		call QNcheck(QN)
		call U1QuantumNumber2%setQN(QN)
		call U1QuantumNumber2%setDeg(Deg)
		call U1QuantumNumber2%setRule(rule)
		return
	end function
	type(QuanNum) function U1QuantumNumber2_(QN1,QN2,deg,rule)result(U1QuantumNumber2)
		integer,intent(in)::deg(:),rule
		real*4,intent(in)::QN1,QN2
		real*4,allocatable::QN(:)
		integer::lenQN,i
		if(QN1.gt.QN2)then
			call writemess('ERROR in input QN1,QN2, QN1 > QN2')
			call error_stop
		end if
		lenQN=int(anint(QN2))-int(anint(QN1))+1
		if(lenQN.ne.size(deg))then
			call writemess('ERROR in input deg, size(deg) .ne. lenQN')
			call error_stop
		end if
		allocate(QN(lenQN))
		do i=1,lenQN
			QN(i)=QN1+(i-1)
		end do
		call U1QuantumNumber2%setQN(QN)
		call U1QuantumNumber2%setDeg(Deg)
		call U1QuantumNumber2%setRule(rule)
		return
	end function
	type(QuanNum) function U1QuantumNumber3(QN,deg,rule,arrow)
		integer,intent(in)::deg(:),rule,arrow
		real*4,intent(in)::QN(:)
		call QNcheck(QN)
		call U1QuantumNumber3%setQN(QN)
		call U1QuantumNumber3%setDeg(Deg)
		call U1QuantumNumber3%setRule(rule)
		call U1QuantumNumber3%setFermiArrow(arrow)
		return
	end function
	type(QuanNum) function U1QuantumNumber3_(QN1,QN2,deg,rule,arrow)result(U1QuantumNumber3)
		integer,intent(in)::deg(:),rule,arrow
		real*4,intent(in)::QN1,QN2
		real*4,allocatable::QN(:)
		integer::lenQN,i
		if(QN1.gt.QN2)then
			call writemess('ERROR in input QN1,QN2, QN1 > QN2')
			call error_stop
		end if
		lenQN=int(anint(QN2))-int(anint(QN1))+1
		if(lenQN.ne.size(deg))then
			call writemess('ERROR in input deg, size(deg) .ne. lenQN')
			call error_stop
		end if
		allocate(QN(lenQN))
		do i=1,lenQN
			QN(i)=QN1+(i-1)
		end do
		call U1QuantumNumber3%setQN(QN)
		call U1QuantumNumber3%setDeg(Deg)
		call U1QuantumNumber3%setRule(rule)
		call U1QuantumNumber3%setFermiArrow(arrow)
		return
	end function
	
!outQ=Q1+Q2
	type(QuanNum) function QuanFuseOrderFunction(order,Q1,Q2,newRule_)
		type(QuanNum),target,intent(in)::Q1,Q2
		type(Tensor),intent(inout)::order
		integer,optional,intent(in)::newRule_
		integer::newRule
		real*4::meanQN1,meanQN2
		type(QuanNum),pointer::Qp1,Qp2
		type(QuanNum),target::nonZeroQ1,nonZeroQ2
		if(Q1%ifzeroDeg())then
			nonZeroQ1=Q1%NonZeroDeg()
			Qp1=>nonZeroQ1
		else
			Qp1=>Q1
		end if
		if(Q2%ifzeroDeg())then
			nonZeroQ2=Q2%NonZeroDeg()
			Qp2=>nonZeroQ2
		else
			Qp2=>Q2
		end if
		if(present(newRule_))then
			newRule=newRule_
		else
			meanQN1=real(Qp1%getmaxQN()+Qp1%getminQN())/2.
			meanQN2=real(Qp2%getmaxQN()+Qp2%getminQN())/2.
			if(meanQN1.ge.meanQN2)then
				newRule=Qp1%GetRule()
			else
				newRule=Qp2%GetRule()
			end if
		end if
		call QuanFuseOrder(order,QuanFuseOrderFunction,Qp1,Qp2,newRule)
		Qp1=>null()
		Qp2=>null()
		return
	end function

	subroutine QuanFuseOrder(order,outQ,Q1,Q2,newRule)
		type(QuanNum),intent(inout)::outQ
		type(QuanNum),intent(in)::Q1,Q2
		integer,intent(in)::newRule
		type(Tensor),intent(inout)::order
		integer::i,j,k,deg1,deg2,degstart,maxi,maxj,maxk
		integer::rule1,rule2
		real*4::outQN,QN1,QN2,minQN,maxQN
		logical::goon,deg_not_zero_flag
		real*4,allocatable::newQN(:)
		integer,allocatable::deg(:)
		real*4::tempQN(4)
		maxj=Q2%getQNlength()
		maxk=Q1%getQNlength()
		rule1=-1*Q1%GetRule()
		rule2=-1*Q2%GetRule()
		tempQN(1)=Q1%GetQN(1)*rule1+Q2%GetQN(1)*rule2
		tempQN(2)=Q1%GetQN(1)*rule1+Q2%GetQN(maxj)*rule2
		tempQN(3)=Q1%GetQN(maxk)*rule1+Q2%GetQN(1)*rule2
		tempQN(4)=Q1%GetQN(maxk)*rule1+Q2%GetQN(maxj)*rule2
		if(newRule.gt.0)then
			tempQN=-1*tempQN
		endif
		minQN=minval(tempQN)
		maxQN=maxval(tempQN)


		if(maxQN.lt.minQN)then
			call writemess('ERROR maxQN and minQN ,rule of newQN='+newRule,-1)
			call Q1%print()
			call Q2%print()
			call error_stop()
		end if

		maxi=int(anint(maxQN))-int(anint(minQN))+1
		allocate(newQN(maxi))
		allocate(deg(maxi))
		do i=1,maxi
			newQN(i)=minQN+(i-1)
		end do
		if(newQN(maxi).nequ.maxQN)then
			call writemess('ERROR in QuanFuseOrder,1')
			call error_stop
		end if

		if(maxi.eq.1)then
			call QuanFuseOrder_one_QN(order,outQ,Q1,Q2,newRule)
			return
		end if

		call outQ%setQN(newQN)
		call outQ%setRule(newRule)
		
		deg=0
		call order%empty()
		do i=1,maxi
			degstart=0
			do j=1,maxj
				do k=1,maxk
					outQN=outQ%getQN(i)
					QN1=Q1%getQN(k)
					QN2=Q2%getQN(j)
					goon=U1Rule1(outQN,QN1,QN2,newRule,rule1,rule2)
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


	subroutine QuanFuseOrder_one_QN(order,outQ,Q1,Q2,newRule)
		type(QuanNum),intent(inout)::outQ
		type(QuanNum),intent(in)::Q1,Q2
		type(Tensor),intent(inout)::order
		integer,intent(in)::newRule
		integer::i,j,k,deg,deg1,deg2
		real*4::outQN,QN1,QN2
		logical::deg_not_zero_flag
		integer::rule1,rule2
		QN1=Q1%getQN(1)
		QN2=Q2%getQN(1)
		rule1=-Q1%GetRule()
		rule2=-Q2%GetRule()
		outQN=Q1%GetQN(1)*rule1+Q2%GetQN(1)*rule2
		if(newRule.gt.0)then
			outQN=-outQN
		end if
		call outQ%setRule(newRule)
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
		call NonZeroElementU1(NonZeroElement,dimen)
		if(NonZeroElement%getRank().eq.1)then
			call NonZeroElement%resetDim((/1,NonZeroElement%getTotalData()/))
		end if
		return
	end function
	subroutine NonZeroElementU1(Res,dimen)
		type(Tensor),intent(inout)::Res
		class(SymDimension),intent(in)::dimen
		integer,allocatable::indices(:),maxdim(:),mindim(:),allrule(:)
		integer::Rank
		logical::goon,rulefit
		real*4,allocatable::QN(:)
		rank=dimen%getRank()
		allocate(indices(rank))
		allocate(maxdim(rank))
		allocate(mindim(rank))
		allocate(QN(rank))
		allocate(allrule(rank))
		mindim=1
		indices=mindim
		maxdim=dimen%dim()
		goon=.true.
		call Res%empty()
		allrule=Dimen%GetRule()
		do while(goon)
			QN=dimen%QNDim(indices)
			rulefit=U1Rule2(QN,allrule)
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
	

	integer function Parity1(dimen,ith,jth)
		class(SymDimension),intent(in)::dimen
		integer,intent(in)::ith,jth
		integer::QN
		QN=int(anint(Dimen%getQN(ith,jth)))
		if(mod(QN,2).eq.0)then
			Parity1=1
		else
			Parity1=-1
		end if
		return
	end function

	logical function ifParityIndex1(dimen,ith,jth)
		class(SymDimension),intent(in)::dimen
		integer,intent(in)::ith,jth
		integer::QN
		QN=int(anint(Dimen%getQN(ith,jth)))
		if(mod(QN,2).eq.0)then
			ifParityIndex1=.true.
		else
			ifParityIndex1=.false.
		end if
		return
	end function

	logical function ifParityIndex2(dimen,vec,ith,jth)
		class(SymDimension),intent(in)::dimen
		integer,intent(in)::vec(:)
		integer,intent(in)::ith,jth
		integer,allocatable::QN(:)
		integer::rank
		rank=dimen%getRank()
		if(rank.ne.size(vec))then
			call writemess("ERROR in calculate the parity",-1)
			call error_stop
		end if
		allocate(QN(rank))
		QN=int(anint(Dimen%getQN(vec)))
		if(mod(sum(QN(ith:jth)),2).eq.0)then
			ifParityIndex2=.true.
		else
			ifParityIndex2=.false.
		end if
		return
	end function

	logical function ifParityIndex3(dimen,vec)
		class(SymDimension),intent(in)::dimen
		integer,intent(in)::vec(:)
		integer,allocatable::QN(:)
		integer::rank
		rank=dimen%getRank()
		if(rank.ne.size(vec))then
			call writemess("ERROR in calculate the parity",-1)
			call error_stop
		end if
		allocate(QN(rank))
		QN=int(anint(Dimen%getQN(vec)))
		if(mod(sum(QN),2).eq.0)then
			ifParityIndex3=.true.
		else
			ifParityIndex3=.false.
		end if
		return
	end function

	function Parity2(dimen,vec)
		integer,allocatable::Parity2(:)
		class(SymDimension),intent(in)::dimen
		integer,intent(in)::vec(:)
		integer,allocatable::QN(:)
		integer::rank,i
		rank=dimen%getRank()
		if(rank.ne.size(vec))then
			call writemess("ERROR in calculate the parity",-1)
			call error_stop
		end if
		allocate(QN(rank))
		allocate(Parity2(rank))
		QN=int(anint(Dimen%getQN(vec)))
		do i=1,rank
			if(mod(QN(i),2).eq.0)then
				Parity2(i)=1
			else
				Parity2(i)=-1
			end if
		end do
		return
	end function

	subroutine checkSymmetryRule(Rule1,Rule2,lenname1,lenname2)
		integer,intent(in)::Rule1,Rule2
		character(len=*),optional,intent(in)::lenname1,lenname2
		if(rule1*rule2.ge.0)then
			call writemess('ERROR in U(1) Symmetry Rule',-1)
			call writemess('Ruel1='+rule1+',Rule2='+rule2)
			if(present(lenname1))then
				call writemess('lenname1='+lenname1+',lenname2='+lenname2)
			end if
			call error_stop
		end if
		return
	end subroutine

	subroutine reverseSymmetryRule(T)
		Type(SymTensor),intent(inout)::T
		type(Tensor),allocatable::TData(:,:)
		type(Tensor),pointer::p(:,:)
		integer::i,j,m,n
		integer::rule,arrow
		real*4,allocatable::QN(:)
		integer,allocatable::deg(:)
		type(QuanNum)::qunNum
		type(SymDimension)::dimen
		if(T%getRank().ne.2)then
			call writemess('The input should be a matrix in reverseSymmetryRule ')
			call error_stop
		end if
		if(.not.T%out_sample_flag())then
			call writemess('ERROR in  reverseSymmetryRule')
			call error_stop 
		end if
		call T%pointer(p)
		m=T%dim(1)
		n=T%dim(2)
		allocate(TData(m,n))
		do j=1,n
			do i=1,m
				TData(i,j)=p(i,n-j+1)
			end do
		end do
		n=T%getQNlength(2)
		allocate(QN(n))
		allocate(deg(n))
		do i=1,n
			QN(i)=T%getQN(2,n-i+1)
			if((  QN(i).nequ.(0e0) ) )QN(i)=-1*QN(i)
			deg(i)=T%getdeg(2,n-i+1)
		end do
		rule=-1*T%getRule(2)
		Arrow=T%getFermiArrow(2)
		qunNum=quantumnumber(QN,deg,rule,Arrow)
		dimen=[T%Quantumnumber(1),qunNum]
		if(T%getnameFlag().ne.0)then
			do i=1,T%getRank()
				call dimen%setName(i,T%getName(i))
			end do
		end if
		call T%resetdim(dimen)
		do j=1,n
			do i=1,m
				p(i,j)=TData(i,j)
			end do
		end do
		p=>null()
		return
	end subroutine
	

end module

