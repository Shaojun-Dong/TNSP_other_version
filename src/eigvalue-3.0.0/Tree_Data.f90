!This file test OK

!TreeNode: [T]
!DataNode: [D]

!                      [T] <---------------------------------------------DataTree%head point here
!               ________|___________
!             /                     \
!            |                      |
!           [T]                    [T]
!        ____|____              ____|____
!       /         \            /         \
!      |           |          |          |
!     [T]         [T]        [T]        [T]
!      |           |          |          |
!    /   \       /   \      /  \       /   \
!   |     |     |     |    |    |     |     |     
!  [T1]  [T1]  [T3] [T4]  [T5] [T6]  [T7] [T8] 
!    
!T1,T2,T3... TreeNode%DataNode point to DataNode
!            if  TreeNode%DataNode point to null(), it means that there is no data
!
!
!    
!   [D1]<-->[D2]<-->[D3]<-->[D4]<-->[D5]<-->[D6]  <--------------------- DataTree%Dataend point here
!    ^
!    |
!    `------------------------------------------------------------------ DataTree%DataHead point here
!
!  The data Tree use to find the DataNode, For example, lengthIndices=2
!
!  find the data of [1,2,1], use the tree, find the data in T3
!  if T3 point to DataNode, read the data in DataNode. If not, add data to the link of DataNode
!
module Tree_Data
	use Tensor_complex
	use SymDimension_typede
	use Dimension_typede
	use SymTensor_type
	use fermiParityTool
	use fermiOperator
	use fermiTensor
	implicit none
	
	type DataTree
		type(TreeNode),pointer  :: head=>null()
		type(leaf),pointer  :: LeafHead=>null()
		type(leaf),pointer  :: LeafEnd=>null()
		integer::lengthleaf=0
		integer::lengthTreeNode=0
		integer::lengthChild=0
		logical::flag=.false.
	contains
			procedure,public::getFlag	=>TreeGetFlag
			procedure,public::getlengthleaf
			procedure,public::getTreeDeep=>GetlengthTreeNode
			procedure,public::getChildNum=>GetTreechildlength
			procedure,public::initial=>initialTree
			procedure,public::push_back=>Tree_push_back
			generic,public::Leafpointer=>Tree_leaf,Tree_leaf2
			procedure::Tree_leaf,Tree_leaf2
			procedure,public::Leaf=>Tree_leaf3
			procedure,public::deallocate=>cleanTree
			generic,public::delete=>deleteleaf2,deleteleaf1
			procedure::deleteleaf2,deleteleaf1
			procedure,public::leafindex
	end type DataTree 
	
	type  Branch!Use to pointer the next level,it is the data link
		type(BranchNode),pointer::Head=>null()
		type(BranchNode),pointer::End=>null()
	contains
			procedure,public::allocate	=>AllocateBranch
			procedure,public::Branchpointer	=>outBranchpointer
			procedure,public::deallocate=>cleanBranch
	end type Branch
	
	type ,extends (TreeNode) ::BranchNode!The element of Branch
	                                     !NextTreeNode and NextLeaf is the pointer to point to next level
		type(BranchNode),pointer::Next=>null()
		type(BranchNode),pointer::prior=>null()
		type(TreeNode),pointer::NextTreeNode=>null()!pointe to TreeNode
		type(leaf),pointer::NextLeaf=>null()!point to Leaf
		logical::NextTreeNodeFlag=.false.
	end type BranchNode
	
	type TreeNode
		type(TreeNode),pointer::parent=>null()!pointe to prior level
		type(Branch)::child!pointe to next level
		logical::childFlag=.false.!If True, child pointer to some data
		integer::childlength=0!How many TreeNode connect to this level
		integer::index
	contains
			procedure,public::getFlag	=>TreeNodeGetFlag	
			generic,public::push_back=>TreeNode_push_back,TreeNode_push_back2
			procedure::TreeNode_push_back,TreeNode_push_back2
			generic,public::ifNextLevel=>ifNextLevelLeaf,ifNextLevelNode
			procedure::ifNextLevelLeaf,ifNextLevelNode
			procedure,public::Branch=>NodeoutBranchpointer
			procedure,public::getLength=>Getchildlength
			procedure,public::allocate=>AllocateNode
			procedure,public::deallocate=>cleanTreeNode
	end type TreeNode
	
	
	
	type,extends (TreeNode) :: leaf!Use to store data
		type(leaf),pointer::next=>null()
		type(leaf),pointer::prior=>null()
		!*************************************
		!    The data to be stored
		!************************************
		real*8::Ws
	end type leaf
	
	public::assignment(=)
	interface assignment(=)
		module procedure copyLeaf
	end interface
contains
	subroutine copyLeaf(inoutLeaf,inleaf)
		type(leaf),intent(inout)::inoutLeaf
		type(leaf),intent(in)::inleaf
		inoutLeaf%Ws=inleaf%Ws
		inoutLeaf%index=inleaf%index
		return
	end subroutine
	logical function TreeGetFlag(Tree)
		class(DataTree),intent(in)::Tree
		TreeGetFlag=Tree%flag
		return
	end function
	logical function TreeNodeGetFlag(Node)
		class(TreeNode),intent(in)::Node
		select type (Node)
			type is (TreeNode)
				TreeNodeGetFlag=Node%childFlag
			type is (BranchNode)
				TreeNodeGetFlag=Node%NextTreeNodeFlag
		end select
		return
	end function
	integer function Getchildlength(Node)
		class(TreeNode),intent(in)::Node
		Getchildlength=Node%childlength
		return
	end function
	integer function GetTreechildlength(Tree)
		class(DataTree),intent(in)::Tree
		GetTreechildlength=Tree%lengthChild
		return
	end function
	subroutine initialTree(Tree,lengthChild,lengthTreeNode)
		class(DataTree),intent(inout)::Tree
		integer,intent(in)::lengthChild,lengthTreeNode
		Tree%lengthChild=lengthChild
		Tree%lengthTreeNode=lengthTreeNode
		Tree%flag=.true.
		return
	end subroutine
	
	
	integer function getlengthleaf(Tree)
		class(DataTree),intent(in)::Tree
		getlengthleaf=Tree%lengthleaf
		return
	end function
	integer function GetlengthTreeNode(Tree)
		class(DataTree),intent(in)::Tree
		GetlengthTreeNode=Tree%lengthTreeNode
		return
	end function
	
	subroutine AllocateBranch(Br,lenindex)
		class(Branch),intent(inout)::Br
		integer,intent(in)::lenindex
		type(BranchNode),pointer::p,p2
		integer::i
		nullify(p)
		allocate(p)
		Br%head=>p
		p2=>p
		do i=2,lenindex
			nullify(p)
			allocate(p)
			p2%next=>p
			p%prior=>p2
			p2=>p
		end do
		Br%end=>p2
		nullify(p)
		nullify(p2)
		return
	end subroutine
	subroutine AllocateNode(Node,lenindex)
		class(TreeNode),intent(inout)::Node
		integer,intent(in)::lenindex
		call Node%child%allocate(lenindex)
		Node%childFlag=.true.
		return
	end subroutine
	
	subroutine cleanBranch(Br)
		class(Branch),intent(inout)::Br
		type(BranchNode),pointer::p,p1
		p=>Br%head
		do while(associated(p))
			p1=>p%next
			nullify(p%Next)
			nullify(p%prior)
			nullify(p%NextTreeNode)
			nullify(p%NextLeaf)
			p%NextTreeNodeFlag=.false.
			deallocate(p)
			p=>p1
		end do
		nullify(p)
		nullify(p1)
		return
	end subroutine
	
	subroutine cleanTreeNode(Node)
		class(TreeNode),intent(inout)::node
		select type (Node)
			type is (TreeNode)
				nullify(node%parent)
				call node%child%deallocate()
				node%childFlag=.false.
				node%childlength=0
				node%index=0
			type is (leaf)
				nullify(node%parent)
				nullify(node%next)
				nullify(node%prior)
			type is (BranchNode)
				nullify(node%Next)
				nullify(node%prior)
				nullify(node%NextTreeNode)
				nullify(node%NextLeaf)
				node%NextTreeNodeFlag=.false.
		end select
		return
	end subroutine
	
	subroutine cleanTree(Tree)
		class(DataTree),intent(inout)::Tree
		type(leaf),pointer::lp1,lp2
		type(TreeNode),pointer::p1,p2
		integer::i
		logical::goon
		lp1=>Tree%LeafHead
		do while(associated(lp1))
			p1=>lp1%parent
			lp2=>lp1%next
			call lp1%deallocate()
			deallocate(lp1)
			lp1=>lp2
			p1%childlength=p1%childlength-1
			goon=.false.
			if(associated(p1))then
				if(p1%childlength.eq.0)then
					goon=.true.
				end if
			end if
			do while(goon)
				p2=>p1%parent
				call p1%deallocate()
				deallocate(p1)
				p1=>p2
				goon=.false.
				if(associated(p1))then
					p1%childlength=p1%childlength-1
					if(p1%childlength.eq.0)then
						goon=.true.
					end if
				end if
			end do
		end do
		Tree%head=>null()
		Tree%LeafHead=>null()
		Tree%LeafEnd=>null()
		Tree%lengthleaf=0
		Tree%lengthTreeNode=0
		Tree%lengthChild=0
		Tree%flag=.false.
		return
	end subroutine
	
		
	
	
	subroutine outBranchpointer(Br,outpointer,ith)
		class(Branch),intent(in)::Br
		type(BranchNode),pointer,intent(inout)::outpointer
		integer,intent(in)::ith
		integer::i
		outpointer=>Br%head
		do i=2,ith
		!	if(associated(outpointer%next))then
				outpointer=>outpointer%next
		!	else
		!		call writemess('ERROR in Finding the pointer')
		!		call error_stop()
		!	end if
		end do
		return
	end subroutine
	
	
	subroutine NodeoutBranchpointer(Node,outpointer,ith)
		class(TreeNode),intent(in)::Node
		type(BranchNode),pointer,intent(inout)::outpointer
		integer,intent(in)::ith
		integer::i
		!if(.not.Node%childFlag)then
		!	call writemess('ERROR. Do no allocate data')
		!	call error_stop()
		!end if
		call Node%child%Branchpointer(outpointer,ith)
		return
	end subroutine


!return 0, do not allocate
!return -1, point to null
!return 1 pointer to next level	
	integer function ifNextLevelNode(Node,outpointer,ith)result(ifNextLevel)
		class(TreeNode),intent(in)::Node
		type(TreeNode),pointer,intent(inout)::outpointer
		integer,intent(in)::ith
		type(BranchNode),pointer::p
		if(.not.Node%childFlag)then
			ifNextLevel=0
			nullify(outpointer)
			return
		end if
		call Node%child%Branchpointer(p,ith)
		if(associated(p%NextTreeNode))then
			outpointer=>p%NextTreeNode
			ifNextLevel=1
		else
			ifNextLevel=-1
			nullify(outpointer)
		end if
		nullify(p)
		return
	end function
	integer function ifNextLevelLeaf(Node,outpointer,ith)result(ifNextLevel)
		class(TreeNode),intent(in)::Node
		type(leaf),pointer,intent(inout)::outpointer
		integer,intent(in)::ith
		type(BranchNode),pointer::p
		if(.not.Node%childFlag)then
			ifNextLevel=0
			nullify(outpointer)
			return
		end if
		call Node%child%Branchpointer(p,ith)
		if(associated(p%NextLeaf))then
			outpointer=>p%NextLeaf
			ifNextLevel=1
		else
			ifNextLevel=-1
			nullify(outpointer)
		end if
		nullify(p)
		return
	end function
	
	subroutine TreeNode_push_back(Node,ith,Node2)
		class(TreeNode),intent(inout)::Node
		integer,intent(in)::ith
		type(TreeNode),pointer,intent(in)::Node2
		type(BranchNode),pointer::p
		if(.not.Node%getFlag())then
			call writemess('ERROR in push back TreeNode into TreeNode0, the TreeNode0 should be allocate ')
			call error_stop()
		end if
		call Node%Child%Branchpointer(p,ith)
		p%NextTreeNode=>Node2
		p%NextTreeNodeFlag=.true.
		Node%ChildLength=Node%ChildLength+1
		select type (Node)
			type is (TreeNode)
				Node2%parent=>Node
		end select
		nullify(p)
		return
	end subroutine
	subroutine TreeNode_push_back2(Node,ith,Node2)
		class(TreeNode),intent(inout)::Node
		integer,intent(in)::ith
		type(leaf),pointer,intent(in)::Node2
		type(BranchNode),pointer::p
		if(.not.Node%getFlag())then
			call writemess('ERROR in push back TreeNode into TreeNode0, the TreeNode0 should be allocate ')
			call error_stop()
		end if
		call Node%Child%Branchpointer(p,ith)
		p%NextLeaf=>Node2
		p%NextTreeNodeFlag=.true.
		Node%ChildLength=Node%ChildLength+1
		select type (Node)
			type is (TreeNode)
				Node2%parent=>Node
		end select
		nullify(p)
		return
	end subroutine
	
	subroutine Tree_push_back(Tree,Ws,indices)
		class(DataTree),intent(inout)::Tree
		real*8,intent(in)::Ws
		integer,intent(in)::indices(:)
		integer::i,deep,numChild
		type(leaf),pointer::leafp,lp2
		type(TreeNode),pointer::Nodep,Np
		integer::goon
		if(.not.Tree%getFlag())then
			call writemess('ERROR in push back leaf into Tree, the Tree should be initial ')
			call error_stop()
		end if
		deep=Tree%getTreeDeep()
		numChild=Tree%getChildNum()
		if(maxval(indices).gt.numChild)then
			call writemess('ERROR in push back leaf into Tree, Indices overflow')
			call error_stop()
		end if
		nullify(leafp)
		allocate(leafp)
		leafp%Ws=Ws
		leafp%index=indices(deep)
		if(Tree%getlengthleaf().eq.0)then
			Tree%leafHead=>leafp
			Tree%leafEnd=>leafp
			tree%lengthleaf=tree%lengthleaf+1
			allocate(Nodep)
			Tree%Head=>Nodep
			nullify(Nodep)
		else
			Tree%leafEnd%next=>leafp
			leafp%prior=>Tree%leafEnd
			Tree%leafEnd=>leafp
			tree%lengthleaf=tree%lengthleaf+1
		end if
		
		Nodep=>Tree%Head
		do i=1,deep-1
			goon= Nodep%ifNextLevel(Np,indices(i))
			if (goon.eq.0) then
				nullify(Np)
				allocate(Np)
				Np%index=indices(i)
				call Nodep%allocate(numChild)
				call Nodep%push_back(indices(i),Np)
			else if(goon.eq.-1)then
				nullify(Np)
				allocate(Np)
				Np%index=indices(i)
				call Nodep%push_back(indices(i),Np)
			end if
			Nodep=>Np
		end do
		
		goon= Nodep%ifNextLevel(lp2,indices(deep))
		if (goon.eq.0) then
			call Nodep%allocate(numChild)
			call Nodep%push_back(indices(deep),leafp)
		else if(goon.eq.-1)then
			call Nodep%push_back(indices(deep),leafp)
		else
			call writemess('There already data in the leaf')
			write(*,*)goon
			call error_stop()
		end if
		
		return
	end subroutine

!output the pointer that pointer to indices=[n1,n2,n3,n4]
! if there is the pointer, output .true., and nodepointer point t
! if there is no data, output .false., and nodepointer=>null()
	logical function Tree_leaf(Tree,leafp,indices)
		class(DataTree),intent(in)::Tree
		integer,intent(in)::indices(:)
		type(leaf),pointer,intent(inout)::leafp
		type(TreeNode),pointer::p,p2
		integer::i,goon
		if(.not.Tree%getFlag())then
			call writemess('ERROR. the Tree should be initial ')
			call error_stop()
		end if
		!if(maxval(indices).gt.Tree%getChildNum())then
		!	call writemess('ERROR in getting leaf out of the Tree, Indices overflow')
		!	call error_stop()
		!end if
		Tree_leaf=.true.
		p=>Tree%Head
		if(.not.associated(p))then
			Tree_leaf=.false.
			nullify(leafp)
			return
		end if
		do i=1,Tree%getTreeDeep()-1
			goon=p%ifNextLevel(p2,indices(i))
			if(goon.ne.1)then
				nullify(leafp)
				Tree_leaf=.false.
				return
			end if
			p=>p2
		end do
		goon=p%ifNextLevel(leafp,indices(i))
		if(goon.ne.1)then
			leafp=>null()
			Tree_leaf=.false.
		end if
		return
	end function
	logical function Tree_leaf2(Tree,leafp,ith)
		class(DataTree),intent(in)::Tree
		integer,intent(in)::ith
		type(leaf),pointer,intent(inout)::leafp
		integer::i
		if(.not.Tree%getFlag())then
			call writemess('ERROR. the Tree should be initial ')
			call error_stop()
		end if
		if(ith.gt.Tree%getlengthleaf())then
			Tree_leaf2=.false.
			nullify(leafp)
			return
		end if
		leafp=>Tree%LeafHead
		Tree_leaf2=.true.
		do i=2,ith
			if(associated(leafp%Next))then
				leafp=>leafp%Next
			else
				Tree_leaf2=.false.
				leafp=>null()
				return
			end if
		end do
		return
	end function
	
	type(leaf) function Tree_leaf3(Tree,ith)
		class(DataTree),intent(in)::Tree
		integer,intent(in)::ith
		type(leaf),pointer::leafp
		integer::i
		if(.not.Tree%getFlag())then
			call writemess('ERROR. the Tree should be initial ')
			call error_stop()
		end if
		if(ith.gt.Tree%getlengthleaf())then
			call writemess('ERROR in getting leaf out of the Tree, ith overflow')
			call error_stop()
		end if
		leafp=>Tree%LeafHead
		do i=2,ith
			if(associated(leafp%Next))then
				leafp=>leafp%Next
			else
				call writemess('ERROR in getting leaf out of the Tree, NO.1')
				call error_stop()
			end if
		end do
		Tree_leaf3=leafp
		return
	end function
	function leafindex(Tree,ith)
		integer,allocatable::leafindex(:)
		class(DataTree),intent(in)::Tree
		integer,intent(in)::ith
		logical::goon
		integer::i,leng
		type(leaf),pointer::lp
		type(TreeNode),pointer::Np
		leng=Tree%getTreeDeep()
		allocate(leafindex(leng))
		goon=Tree%LeafPointer(lp,ith)
		if(.not.goon)then
			call writemess('ERROR in getting leafindex, ith may larger then length')
			call error_stop()
		end if
		leafindex(leng)=lp%index
		Np=>lp%parent
		do i=2,leng
			leafindex(leng-i+1)=Np%index
			Np=>Np%parent
		end do
		return
	end function
	
	
	subroutine deleteleaf1(Tree,ith)
		class(DataTree),intent(inout)::Tree
		integer,intent(in)::ith
		type(leaf),pointer::leafp,lp1,lp2
		type(TreeNode),pointer::p1,p2
		integer::i
		logical::goon
		if(.not.Tree%getFlag())then
			call writemess('ERROR. the Tree should be initial ')
			call error_stop()
		end if
		if(ith.gt.Tree%getlengthleaf())then
			call writemess('ERROR in getting leaf out of the Tree, ith overflow')
			write(*,*)ith,Tree%getlengthleaf()
			call error_stop()
		end if
		leafp=>Tree%LeafHead
		do i=2,ith
			if(associated(leafp%Next))then
				leafp=>leafp%Next
			else
				call writemess('ERROR in getting leaf out of the Tree, NO.1')
				call error_stop()
			end if
		end do
		
		p1=>leafp%parent
		lp1=>leafp%prior
		lp2=>leafp%next
		call leafp%deallocate()
		deallocate(leafp)
		nullify(leafp)
		if(associated(lp1))then
			lp1%next=>lp2
		else!The first leaf
			Tree%LeafHead=>lp2
		end if
		if(associated(lp2))then
			lp2%prior=>lp1
		else!The last leaf
			Tree%LeafEnd=>lp1
		end if
		tree%lengthleaf=tree%lengthleaf-1
		
		p1%childlength=p1%childlength-1
		goon=.false.
		if(p1%childlength.eq.0)goon=.true.
			
		do while(goon)
			p2=>p1%parent
			call p1%deallocate()
			deallocate(p1)
			p1=>p2
			goon=.false.
			if(associated(p1))then
				p1%childlength=p1%childlength-1
				if(p1%childlength.eq.0)then
					goon=.true.
				end if
			end if
		end do
		if(tree%lengthleaf.eq.0)nullify(Tree%head)
		return
	end subroutine
	
	subroutine deleteleaf2(Tree,indices)
		class(DataTree),intent(inout)::Tree
		integer,intent(in)::indices(:)
		type(leaf),pointer::leafp,lp1,lp2
		type(TreeNode),pointer::p1,p2
		integer::i
		logical::goon
		if(.not.Tree%getFlag())then
			call writemess('ERROR. the Tree should be initial ')
			call error_stop()
		end if
		goon=Tree%LeafPointer(leafp,indices)
		if(.not.goon)return
		
		p1=>leafp%parent
		lp1=>leafp%prior
		lp2=>leafp%next
		call leafp%deallocate()
		deallocate(leafp)
		nullify(leafp)
		if(associated(lp1))then
			lp1%next=>lp2
		else!The first leaf
			Tree%LeafHead=>lp2
		end if
		if(associated(lp2))then
			lp2%prior=>lp1
		else!The last leaf
			Tree%LeafEnd=>lp1
		end if
		tree%lengthleaf=tree%lengthleaf-1
		
		p1%childlength=p1%childlength-1
		goon=.false.
		if(p1%childlength.eq.0)goon=.true.
			
		do while(goon)
			p2=>p1%parent
			call p1%deallocate()
			deallocate(p1)
			p1=>p2
			goon=.false.
			if(associated(p1))then
				p1%childlength=p1%childlength-1
				if(p1%childlength.eq.0)then
					goon=.true.
				end if
			end if
		end do
		if(tree%lengthleaf.eq.0)nullify(Tree%head)
		return
	end subroutine

















end module



!  May be the link will be faster than Tree, see the example below:
!
!	type(DataTree)::Tree
!	type(leaf),pointer::lea
!	type(leaf)::lea2
!	logical::goon
!	integer::ide(100),ide2(100)
!	call Tree%initial(4,100)
!	ide=1
!	call Tree%push_back(0d0,ide)
!	ide=2
!	call Tree%push_back(2d0,ide)
!	ide=3
!	call Tree%push_back(3d0,ide)
!	call cpu_time(time1)
!	do i=1,999999
!		ide=1
!		goon=Tree%LeafPointer(lea,ide)
!		ide=2
!		goon=Tree%LeafPointer(lea,ide)
!		ide=4
!		goon=Tree%LeafPointer(lea,ide)
!	end do
!	call cpu_time(time2)
!	write(*,*)time2-time1
!	ide2=1
!	ide=1
!	ide2(59)=2
	
!	call cpu_time(time1)
!	do i=1,999999
!		goon=.true.
!		do j=1,100
!			goon=goon.and.(ide(j).eq.ide2(j))
!			if(.not.goon)exit
!		end do
!	end do
!	call cpu_time(time2)
!	write(*,*)time2-time1
