!to be modify function fuse
!The permutation can be more fast,permute the block and at the same time permute the data
!in the block
! .equ. may not right(SymDimenion)
!It is sure that there is something wrong in symfuseTensor,when QN is large, it may gose out of memory.
!But it is ok for QN<6.
module SymTensor_complex
	use Tensor_complex
	use SymDimension_typede
	implicit none
!index of the Tensor is 1,2,3...N
!index of the physis Tensor should be related to the quantum number(QN)
!If the index of the Tensor is(QN=0.5)
!index of Tensor		[1]
!index of SymTen		[0]
!
!1/2 spin:
!index of Tensor		 [1,2]  
!index of SymTen		[-0.5,0.5]
!
!1 spin:
!index of Tensor		 [1,2,3]
!index of SymTen		[-1,0,1]

!3/2 spin:
!index of Tensor		 [1,2,3,4]
!index of SymTen		[-3/2,-1/2,1/2,3/2]
!
!index of Tensor		 [1,2,...,N]
!index of SymTen		[-(N-1)*Qnum_,-(N-1)*Qnum_+Delta_Q_,..,(N-1)*Qnum_-Delta_Q_,(N-1)*Qnum_]
	type SymTensor
		integer,private :: totalblock=0!The number of the total partitioned Tensors,that is the non-zero block
		integer,private :: totalData=0!size(block)
		type(Tensor),allocatable::block(:)!the partitioned Tensor
		type(SymDimension),private::TenDim
		integer,private :: rank=0!length of Dimension.TenDim%Dimsize
		logical,private :: flag=.false. !if flag=true ,it means there are data in SymTensor
	end type SymTensor
	!the rule for partition,in U(1),the rule is \sum S_in=\sum S_out
	!==>\sum S_in - \sum S_out =0
	!==>rule(1)*(N_1-1)*Qnum_+rule(1)*(N_2-1)*Qnum_+..+rule(n)*(N_n-1)*Qnum_=0
	!N_i is the ith dimension
!
!   s1 -->--[Tensor]--->-- s3    : a [2*3*4] Tensor,the non-zero partitioned Tensor is 1*s1+1*s2-1*s3=0
!              |                rule:[1,1,-1] ,1 mean the QN gose in, and -1 mean gose out
!              ^                In svd ,use for determine the new dimenion of QTensor
!              | s2
!
!example of degeneracy:
!QN:1,0.5,1.5
!   1 -->--[Tensor]--->-- 1.5 
!              |        
!              ^ 
!              | 0.5  
!first index of SymTen[-1,0,1],degeneracy=[1,2,1]( [d1_1,d1_2,d1_3])
!second index of SymTen[-0.5,0.5],degeneracy=[1,1]( [d2_1,d2_2])
!third index of SymTen[-1.5,-0.5,0.5,1.5],degeneracy=[1,2,2,1]( [d3_1,d3_2,d3_3,d3_4])
!The non-zero partitioned Tensors of T_{s1,s2,s3} are
!T_{-1,-0.5,-1.5} a Tesnor of 1*1*1 ,or d1_1*d2_1*d3_1
!T_{-1, 0.5,-0.5} a Tesnor of 1*1*2 ,or d1_1*d2_2*d3_2
!T_{ 0,-0.5,-0.5} a Tesnor of 2*1*2 ,or d1_2*d2_1*d3_2
!T_{ 0, 0.5, 0.5} a Tesnor of 2*1*2 ,or d1_2*d2_2*d3_3
!T_{ 1,-0.5, 0.5} a Tesnor of 1*1*2 ,or d1_3*d2_1*d3_3
!T_{ 1, 0.5, 1.5} a Tesnor of 1*1*1 ,or d1_3*d2_2*d3_3
	logical,save::testprint=.false.
	interface inde_counter
		module procedure inde_counter_real
		module procedure inde_counter_int
	end interface
	interface SymTensorQN
		module procedure SymTensorQN1
		module procedure SymTensorQN2
		module procedure SymTensorQN1_name
		module procedure SymTensorQN2_name
	end interface	
	interface Symresetdim
		module procedure Symresetdim1
		module procedure Symresetdim2
	end interface	
	interface symTenproduct
		module procedure symTenproduct_noName
		module procedure symTenproduct_Name!In put dimension as character
		module procedure symTenproduct_name_int
		module procedure symTenproduct_int_name
		module procedure symTenproduct_old
		module procedure symTenproduct_name_rename1
		module procedure symTenproduct_name_rename2
	end interface	
	interface SymDegeneracy
		module procedure SymDegeneracy1
		module procedure SymDegeneracy2
		module procedure SymDegeneracy3
		module procedure SymDegeneracy4
	end interface
	interface setbolckData
		module procedure setbolckData1
		module procedure setbolckData2
		module procedure setbolckData3
		module procedure setbolckData4
	end interface	
	interface symeye
		module procedure symeye_one
		module procedure symeye_one_dim
	end interface	
	interface Symgenerate
		module procedure Symgenerate_U1_Dim
		module procedure Symgenerate_U1
		module procedure Symgenerate_U1_QN1
		module procedure Symgenerate_U1_QN2
	end interface	
	interface symSVDcutoff
		module procedure symSVDcutoff1
		module procedure symSVDcutoff2
	end interface	
	interface SymfuseTensor
		module procedure fuseTensor_noName
		module procedure fuseTensor_Name
		module procedure fuseTensorDim_noName
		module procedure fuseTensorDim_Name
		module procedure fuseTensor_Name_array
		module procedure fuseTensor_noName_array
	end interface	
	
	interface setSymTensorName
		module procedure setSymTensorName2!set all the name in Tensor call setTensorName(T,'A') or setTensorName(T,'A12')
		module procedure setSymTensorName3!modify the name of ith index
		module procedure setSymTensorName4!modify the index whose name is old_w to new_w
	end interface
	
	interface SymTenNameorder
		module procedure SymTenNameorder1
		module procedure SymTenNameorder2
	end interface
!*******************   permutation   ****************************	
	interface operator(.p.)
		module procedure Sympermute_rank2 !permute the tensor whose rank is 2,if rank=1 ,do nothing
		module procedure Sympermute_rank3 !permute the tensor whose rank is 3
		module procedure Sympermutation	 !permute the tensor of any rank,give the new order of the dimension
												 !		If operate on a big Tensor,use permute_rank2 or permute_rank3,
												 !	 they are faster.if rank>3,use contract to reshape
		module procedure Sympermutation_name!input a character(:) as the new order
	end interface
	interface operator(.pf.)
		module procedure Sympermutefo!permute the inde index to the first
										  !T_{1,2,3,..,i,..,n},permutefo(T,i)=_{i,1,2,3,..,i-1,i+1,..,n}
										  !T_{1,2,3,..,j,..,i,.,k,...,n},permutefo_vec(T,(/i,j,k/))=_{i,j,k,1,2,3,...,n}
										  ! note the output will decompose all the dimension
		module procedure Sympermutefo_vec
		module procedure Sympermutefo_name
		module procedure Sympermutefo_vec_name
	end interface
	interface operator(.pb.)!	T_{1,2,3,..,i,..,n},permuteback(T,i)=_{1,2,3,..,i-1,i+1,..,n,i}
									!!T_{1,2,3,..,j,..,i,.,k,...,n},permuteback_vec(T,(/i,j,k/))=_{1,2,3,...,n,i,j,k}
									! note the output will decompose all the dimension
		module procedure Sympermuteback
		module procedure Sympermuteback_vec
		module procedure Sympermuteback_name
		module procedure Sympermuteback_vec_name
	end interface
	interface operator(.pi.)!T_{1,2,3,..,i,..,n},permuteInde(T,i)=_{2,3,..,i,1,i+1,..,n}
									! note the output will decompose all the dimension
		module procedure SympermuteInde
	end interface
	
	interface operator(.pbi.)!T_{1,2,3,..,i,..,n},permuteInde(T,i)=_{1,2,3,..,i-1,n,i,i+1,..,n-1}
									! note the output will decompose all the dimension
		module procedure SympermutebackInde
	end interface		
	
!**************** contract   ****************
!		combine two index of the Tensor,which is index1 index1+1,..,index2-1,index2
!		if index2 larger than rank,the function will contract the index of index1 -> rank
	interface operator(.c.)
		module procedure Symcontract!combine two index of the Tensor,which is con_index and con_index+1
		module procedure Symcontractsv!combine two index of the Tensor,input vector specify the to be combined index
	end interface
!************************************************************	
!	 compose_decompse do the contract or decompose on the Tensor 
!	(T.cd.Operation)Operation is a matrix,the first element of every row specify
!		1:contracts
!		2:decompose
!		3:decomposeAll
!	other element of every row is input parameter of the function
!		this function run for big T,do as less as possible operation on the
!	 Tensor_Data of T.
	interface operator(.cd.)
		module procedure Symcompose_decompse
	end interface
!dimOperation(Tensor,Operation):
!	 dimOperations do the contract or decompose on the Tensor with
!	  no change on the Tensot_data in the Tensor.
!	Operation is a matrix or vector,the first element of every row specify
!		1:contracts
!		2:decompose
!		3:decomposeAll
!	other element of every row is input parameter of the function
!		this subroutine run for big Tensor,do as less as possible operation on the
!	 Tensor_Data of Tensor.
	interface symdimOperation
		module procedure symdimOperations
		module procedure symdimOperation1
	end interface
!************************************************************		
! decompose the de_index index of the Tensor into n(1),n(2)
!		for example the de_index index is (1,2,3,4,..inde,inde+1,...rank)
!		(1,2,3,4,..inde,inde+1,...rank)-->(1,2,3,4,..inde),(inde+1,...rank)		
!		if inde larger than rank ,the function will return no change	
	interface operator(.dc.)
		module procedure Symdecompose1
		module procedure SymdecomposeAll
		module procedure Symdecomposev!decompose index of the Tensor,input vector specify the de_index and index
	end interface	
	interface operator(.i.)
		module procedure SymElement
		module procedure SymElement2
	end interface	
	interface operator(.iname.)!the index name
		module procedure SymoutIndexName
		module procedure SymoutAllIndexName
	end interface
	interface operator(.Tname.)!the Tensor name
		module procedure SymoutTensorName
		module procedure SymoutAllTensorName
	end interface
	
	interface operator(.dim.)
		module procedure SymgetTenDim_Namei!output the ith dimension of the Tensor, input the name of the dimension
													!If can not find , output 0
		module procedure SymgetTenDim_i!output the i dimension of the Tensor,output an integer
		module procedure SymgetTenDim!output the Dimension of Tensor,return type(Dimension)
											!It can use to commit the assignment to a vector
											!vector=type(Dimension),if the vector is allocatable
	end interface
	interface operator(.subDim.)
		module procedure SymgetTenSubDim!output the ith dimension of the Tensor,return type(Dimension)
		module procedure SymgetTenSubDim_name
		module procedure SymgetTenDim!output the Dimension of Tensor,return type(Dimension)
											!It can use to commit the assignment to a vector
											!vector=type(Dimension)
	end interface
	
	interface operator(+)
		module procedure Symadd
	end interface
	interface operator(-)
		module procedure Symminus
	end interface
	
	interface operator(*)
		module procedure SymProductTensor
		module procedure Symmultiply_number!A Tensor times a number,T*num
		module procedure Symmultiply_number_!num*T
		module procedure Symmultiply_real!A Tensor times a real*8 number
		module procedure Symmultiply_real_!a real*8 number times a Tensor
	end interface
	interface operator(/)
		module procedure Symdivide_number
		module procedure Symdivide_numberReal
	end interface
	
	interface operator(.x.)
		module procedure Symdot! dot product conjugating the first vector,The Tensor will be regard as a vector
	end interface
	
	interface operator(.numx.)
		module procedure numdot!the same as dot product but do not conjugating the first vector,The Tensor will be regard as a vector
	end interface
	
	interface operator(.H.)
		module procedure SymHtranspose! conjugateTranspose! one or two dimension case
											! If input rank=1,the result is the same as conjugate
	end interface
	interface operator(.HR.)
		module procedure SymHtransposeReverse! conjugateTranspose and reverse the rule of the SymTensor! one or two dimension case
											! If input rank=1,the result is the same as conjugate
	end interface
	interface operator(.con.)
		module procedure Symconjugate! conjugate
	end interface
	interface operator(.conR.)
		module procedure SymconjugateReverse! conjugate and reverse the rule of the SymTensor
	end interface
	interface operator(.R.)
		module procedure reverseTensorRule! reverse the rule of the SymTensor the same as reverseRule
	end interface
	interface operator(.equ.)
		module procedure equal_array_real
		module procedure Symequal_of_Tensor
		!module procedure equal_array_int
	end interface
	
	interface assignment(=)
		module procedure assignmentSymTen!T1=T2 ,both T1 and T2 are Tensor
		module procedure SymcopyTensordim1
		module procedure SymcopyTensordim2
		module procedure SymcopyTensordim3
		module procedure U1Tensor2Ten
		module procedure SymassignmentTenArray
	end interface
contains
	logical function U1Rule(QNumber,rule)
		real*4,intent(in)::QNumber(:)
		integer,intent(in)::rule(:)
		real*8::temp
		temp=sum(QNumber*rule)
		if(abs(temp).lt.zero_error_)then
			U1Rule=.true.
		else
			U1Rule=.false.
		end if
		return
	end function
	logical function equal_array_real(array1,array2)result(res)
		real*4,intent(in)::array1(:),array2(:)
		integer::i,len1,len2
		res=.true.
		len1=size(array1)
		len2=size(array2)
		if(len1.ne.len2)then
			res=.false.
			return
		end if
		do i=1,len1
			if(abs(array1(i)-array2(i)).gt.zero_error_)then
				res=.false.
				return
			end if
		end do
		return
	end function
	logical function equal_array_int(array1,array2)result(res)
		integer,intent(in)::array1(:),array2(:)
		integer::i,len1,len2
		res=.true.
		len1=size(array1)
		len2=size(array2)
		if(len1.ne.len2)then
			res=.false.
			return
		end if
		do i=1,len1
			if(array1(i).ne.array2(i))then
				res=.false.
				return
			end if
		end do
		return
	end function
!!if[s1,s2,...,sn]=[max_s1,max_s2,..,max_sn] output false and return
!else output ture and do the code below
![s1,s2,...,sn]-->[s1+1,s2,..,sn],delta=1
!if s1+1>max_s1,s1=min_s1 and s2=s2+1
!if s2+1>max_s2,s2=min_s2 and s3=s3+1
!...
!inde are [s1,s2,...,sn]
!minindex are [min_s1,min_s2,..,min_sn]
!maxindex are [max_s1,max_s2,..,max_sn]
	logical function inde_counter_real(inde,minindex,maxindex,delta) result(inde_counter)
		real*4,intent(inout)::inde(:)
		real*4,intent(in)::minindex(:),maxindex(:),delta
		integer::indexlen,i
		indexlen=size(inde)
		if(inde.equ.maxindex) then
			inde_counter=.false.
			return
		end if
		i=1
		inde_counter=.true.
		do i=1,indexlen
			inde(i)=inde(i)+delta
			if(inde(i).gt.maxindex(i))then
				inde(i)=minindex(i)
			else
				exit
			end if
		end do
		return
	end function
	logical function inde_counter_int(inde,minindex,maxindex,delta)result(inde_counter)
		integer,intent(inout)::inde(:)
		integer,intent(in)::minindex(:),maxindex(:),delta
		integer::indexlen,i
		indexlen=size(inde)
		if(equal_array_int(inde,maxindex)) then
			inde_counter=.false.
			return
		end if
		i=1
		inde_counter=.true.
		do i=1,indexlen
			inde(i)=inde(i)+delta
			if(inde(i).gt.maxindex(i))then
				inde(i)=minindex(i)
			else
				exit
			end if
		end do
		return
	end function
!In U(1), if quantum number is S, then the dimension of this freedom is 2S+1
	subroutine U1QNtoMaxIndex(outvec,QN)
		integer,intent(inout)::outvec(:)
		real*4,intent(in)::QN(:)
		outvec=2.*QN+1.
		return
	end subroutine	
	!***************** assignment *********************
	subroutine assignmentSymTen(T,T2)
		type(SymTensor),intent(inout) ::T
		type(SymTensor),intent(in) :: T2
		integer::length,i
		if(T2%flag) then
			T%rank=T2%rank
			T%TenDim=T2%TenDim
			length=T2%TotalData
			if(allocated(T%block)) then
				if(size(T%block).lt.length) then!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!If old data longer than new,do not allocate
					deallocate(T%block)
					allocate(T%block(length))
				end if
			else
				allocate(T%block(length))
			end if
			do i=1,length
				T%block(i)=T2%block(i)
			end do
			T%totalData=T2%totalData
			T%totalblock=T2%totalblock
			T%flag=T2%flag
			return
		end if
		T%flag=T2%flag
		return
	end subroutine
	subroutine SymassignmentTenArray(T,T2)
		type(SymTensor),intent(inout) ::T(:)
		type(SymTensor),intent(in) :: T2(:)
		integer::length,i
		length=size(T2)
		if(size(T).lt.length)then
			write(*,*)"ERROR in assignment of two SymTensor array "
			write(*,*)"T1(:)=T2(:),size(T1)<size(T2)"
			write(*,*)size(T),length
			stop
		end if
		do i=1,length
			T(i)=T2(i)
		end do
		return
	end subroutine
!Do not allocate the Tensor data
!copy the Tensor_Data to a vector	
	subroutine SymcopyTensordim1(Vec,T)
		type(Tensor),intent(out) ::Vec(:)
		type(SymTensor),intent(in) :: T
		integer::length,i
		length=T%TotalData
		if(size(Vec).lt.length) then
			write(*,*)"The array can not store the Data"
			write(*,*)"Length of the array",size(Vec)
			write(*,*)"Length of the Tensor",length
			stop
		end if
		call SZCOPY(length,T%block,1,Vec,1)
		return
	end subroutine
!if rank=2,then copy the Tensor_Data to a mat
!Do not allocate the Tensor data
! mat is a matrix	
	subroutine SymcopyTensordim2(Mat,T)
		type(Tensor),intent(out) ::Mat(:,:)
		type(SymTensor),intent(in) :: T
		integer::m,n
		if(SymgetRank(T).ne.2)then
			write(*,*)"T is not a rank-2 Tensor"
			write(*,*)"Can not copy to a matrix"
			stop
		end if
		m=T.dim.1
		n=T.dim.2
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2)
			write(*,*)"size of the Tensor",m,n
			stop
		end if
		call SZCOPY(m*n,T%block,1,Mat,1)
		return
	end subroutine
!if rank=3,then copy the Tensor_Data to a mat
!Do not allocate the Tensor data
!mat is a 3 dimension array	
	subroutine SymcopyTensordim3(Mat,T)
		type(Tensor),intent(out) ::Mat(:,:,:)
		type(SymTensor),intent(in) :: T
		integer::m,n,l
		if(SymgetRank(T).ne.3)then
			write(*,*)"T is not a rank-3 Tensor"
			write(*,*)"Can not copy to the array of dimension 3"
			stop
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3)
			write(*,*)"size of the Tensor",m,n,l
			stop
		end if
		call SZCOPY(size(Mat),T%block,1,Mat,1)
		return
	end subroutine	
	subroutine SZCOPY(N,A,incA,B,incB)!B=A  !zcopy for block data,Not finished but it could run
		integer,intent(in)::N,INCA,incB
		type(Tensor),intent(inout)::B(*)
		type(Tensor),intent(in)::A(*)
		integer::i
		if( incA.eq.1 .and. incB.eq.1 ) then
			do i=1,N
				B(i)=A(i)
			end do
			return
		end if
		write(*,*)"Do not finish this part"
		stop
	end 	subroutine
	subroutine SymstoreTenData(ST,TData,length)
		type(SymTensor),intent(inout) :: ST
		type(Tensor),intent(in)::TData(*)
		integer,intent(in)::length
		if(allocated(ST%block)) then
			if(size(ST%block).lt.length) then!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!If old data longer than new,do not allocate
				deallocate(ST%block)
				allocate(ST%block(length))
			end if
		else
			allocate(ST%block(length))
		end if
		call SZCOPY(length,TData,1,ST%block,1)
		return
	end subroutine
!*************** cleanTensor  *****************
	subroutine cleanSymTensor(T)
		Type(SymTensor),intent(inout)::T
		integer::i
		T%rank=0
		call cleanSymDimension(T%TenDim)
		if(allocated(T%block)) then
			do i=1,size(T%block)
				call cleanTensor(T%block(i))
			end do
			deallocate(T%block)
		end if
		T%totalData=0
		T%totalblock=0
		T%flag=.false.
		return
	end subroutine
	subroutine allocateSymTensor(T,dimen)
		Type(SymTensor),intent(inout)::T
		type(SymDimension),intent(in)::dimen
		if(Symgetflag(T))then
			write(*,*)"Can not allocate SymTensor"
			write(*,*)"Should cleanSymTensor fiset"
			stop
		end if
		T%rank=symDimsize(dimen)
		T%TenDim=dimen
		T%totalData=Symouttotaldata(dimen)
		T%totalBlock=0
		allocate(T%block(T%totalData))
		T%flag=.true.
		return
	end subroutine
!*********************  getRank	 **********************
	integer function SymgetRank(T)
		type(SymTensor),intent(in) :: T
		SymgetRank=T%rank
	end function
!*********************  getTotalData	 **********************
	integer function SymgetTotalData(T)
		type(SymTensor),intent(in) :: T
		SymgetTotalData=T%totalData
	end function	
	integer function SymgetTotalblock(T)
		type(SymTensor),intent(in) :: T
		SymgetTotalblock=T%totalblock
	end function	
!*********************  getFlag	 **********************
	logical function SymgetFlag(T)
		type(SymTensor),intent(in) :: T
		SymgetFlag=T%flag
	end function		
!*********************  get Tensor dimenison *********************
	integer function SymgetTenDim_i(T,inde)
		type(SymTensor),intent(in) :: T
		integer,intent(in) :: inde
		SymgetTenDim_i=SymDim_i(T%TenDim,inde)
		return
	end function
	integer function SymgetTenDim_Namei(T,w)
		type(SymTensor),intent(in) :: T
		character(len=*),intent(in)::w
		integer :: inde
		inde=SymNameorder(T%TenDim,w)
		if(inde.eq.0)then
			SymgetTenDim_Namei=0
			return
		end if
		SymgetTenDim_Namei=SymDim_i(T%TenDim,inde)
		return
	end function
	type(SymDimension) function SymgetTenDim(T)
		type(SymTensor),intent(in) :: T
		SymgetTenDim=T%TenDim
		return
	end function
	type(SymDimension) function SymgetTenSubDim(T,inde)
		type(SymTensor),intent(in)::T
		integer,intent(in)::inde
		SymgetTenSubDim=T%TenDim.sub.inde
		return
	end function
	type(SymDimension) function SymgetTenSubDim_name(T,w)
		type(SymTensor),intent(in)::T
		CHARACTER(len=*),intent(in)::w
		integer::inde
		inde=symTenNameorder(T,w)
		SymgetTenSubDim_name=T%TenDim.sub.inde
		return
	end function
	
	real*4 function  SymTensorQN1(T,ith,jth)
		type(SymTensor),intent(in)::T
		integer,intent(in)::ith,jth
		SymTensorQN1=SymQNum(T%TenDim,ith,jth)
		return
	end function
	real*4 function  SymTensorQN2(T,ith)
		type(SymTensor),intent(in)::T
		integer,intent(in)::ith
		SymTensorQN2=SymQNum(T%TenDim,ith)
		return
	end function
real*4 function  SymTensorQN1_name(T,w,jth)
		type(SymTensor),intent(in)::T
		CHARACTER(len=*),intent(in)::w
		integer,intent(in)::jth
		integer::inde
		inde=symTenNameorder(T,w)
		SymTensorQN1_name=SymQNum(T%TenDim,inde,jth)
		return
	end function
	real*4 function  SymTensorQN2_name(T,w)
		type(SymTensor),intent(in)::T
		CHARACTER(len=*),intent(in)::w
		integer::inde
		inde=symTenNameorder(T,w)
		SymTensorQN2_name=SymQNum(T%TenDim,inde)
		return
	end function
!the dimensionis[2,2,3] the quantum number are [(-0.5,0.5),(-0.5,0.5),(-1,0,1)]
!the degeneracy are [(1,1),(1,1),(1,2,1)]
!output the degeneracy the dimension [1,1,3],That is [1,1,1]
	subroutine SymDegeneracy1(outvec,T,inde)
		type(SymTensor),intent(in)::T
		integer,intent(in)::inde(:)
		integer,intent(out)::outvec(:)
		call SymDimDegeneracy(outvec,T%TenDim,inde)
		return
	end subroutine
!the dimensionis[2,2,3] the quantum number are [(-0.5,0.5),(-0.5,0.5),(-1,0,1)]
!the degeneracy are [(1,1),(1,1),(1,2,1)]
!output the degeneracy the dimension [-0.5,-0.5,1],That is [1,1,1]
	subroutine SymDegeneracy2(outvec,T,QN)
		type(SymTensor),intent(in)::T
		integer,intent(out)::outvec(:)
		real*4,intent(in)::QN(:)
		call SymDimDegeneracy(outvec,T%TenDim,QN)
		return
	end subroutine
	subroutine SymDegeneracy3(outvec,T,ith,QN)
		type(SymTensor),intent(in)::T
		integer,intent(out)::outvec(:)
		integer,intent(in)::ith(:)
		real*4,intent(in)::QN(:)
		call SymDimDegeneracy(outvec,T%TenDim,ith,QN)
		return
	end subroutine
!the dimensionis[2,2,3] the quantum number are [(-0.5,0.5),(-0.5,0.5),(-1,0,1)]
!the degeneracy are [(1,1),(1,1),(1,2,1)]
!output the degeneracy the dimension 3,That is [1,2,1]
	subroutine SymDegeneracy4(outvec,T,ith)
		type(SymTensor),intent(in)::T
		integer,intent(out)::outvec(:)
		integer,intent(in)::ith
		call SymDimDegeneracy(outvec,T%TenDim,ith)
		return
	end subroutine
	
	integer function  SymTensorindex(T,ith,QN)
		type(SymTensor),intent(in)::T
		integer,intent(in)::ith
		real*4,intent(in)::QN
		SymTensorindex=Symindexi(T%TenDim,ith,QN)
		return
	end function
	
! reverse the U(1) symmetry rule in the Tensor
! examplex rule [1,1,-1] --> [-1,-1,1]
! the data will not change
	subroutine reverseRule(ST)	
		type(SymTensor),intent(inout) :: ST
		call reverseDimrule(ST%TenDim)
		return
	end subroutine
	type(SymTensor) function reverseTensorRule(ST)	
		type(SymTensor),intent(in) :: ST
		reverseTensorRule=ST
		call reverseDimrule(reverseTensorRule%TenDim)
		return
	end function
!*********************** TensorName   **********************
	subroutine setSymTensorName2(T,TensorName)
		type(SymTensor),intent(inout)::T
		CHARACTER(len=*),intent(in)::TensorName
		call SymDimName(T%TenDim,TensorName)
		return
	end subroutine
	subroutine cleanSymTensorDim(T)
		type(SymTensor),intent(inout)::T
		call cleanSymDimensionName(T%TenDim)
		return
	end subroutine
	subroutine setSymTensorName3(T,ith,w)
		type(SymTensor),intent(inout)::T
		integer,intent(in)::ith
		CHARACTER(len=*),intent(in)::w
		call Symresetname(T%TenDim,ith,w)
		return
	end subroutine
	subroutine setSymTensorName4(T,oldw,neww)
		type(SymTensor),intent(inout)::T
		CHARACTER(len=*),intent(in)::oldw,neww
		call Symresetname(T%TenDim,oldw,neww)
		return
	end subroutine
	CHARACTER(len=len_of_Name+len_of_Name) function SymoutIndexName(T,ith)
		type(SymTensor),intent(in)::T
		integer,intent(in)::ith
		SymoutIndexName=SymoutName(T%TenDim,ith)
		return
	end function
	function SymoutAllIndexName(T)
		CHARACTER(len=len_of_Name+len_of_Name),allocatable::SymoutAllIndexName(:)
		type(SymTensor),intent(in)::T
		integer::i
		allocate(SymoutAllIndexName(T%rank))
		do i=1,T%rank
			SymoutAllIndexName(i)=SymoutName(T%TenDim,i)
		end do
		return
	end function
	CHARACTER(len=len_of_Name) function SymoutTensorName(T,ith)
		type(SymTensor),intent(in)::T
		integer,intent(in)::ith
		SymoutTensorName=SymoutNameTen(T%TenDim,ith)
		return
	end function
	function SymoutAllTensorName(T)
		CHARACTER(len=len_of_Name),allocatable::SymoutAllTensorName(:)
		type(SymTensor),intent(in)::T
		integer::i
		allocate(SymoutAllTensorName(T%rank))
		do i=1,T%rank
			SymoutAllTensorName(i)=SymoutNameTen(T%TenDim,i)
		end do
		return
	end function
	integer function SymTenNameorder1(T,w)
		type(SymTensor),intent(in) :: T
		character(len=*),intent(in)::w
		integer :: inde
		SymTenNameorder1=SymNameorder(T%TenDim,w)
		return
	end function
	function SymTenNameorder2(T,w)
		integer,allocatable::SymTenNameorder2(:)
		type(SymTensor),intent(in) :: T
		character(len=*),intent(in)::w(2)
		integer :: inde
		allocate(SymTenNameorder2(size(w)))
		SymTenNameorder2=SymNameorder(T%TenDim,w)
		return
	end function	
!input dimension [1,1,2,1,3,1,1,4,1]
!output dimenison [2,3,4]	
	subroutine SymRNTensorDim(T)
		type(SymTensor),intent(inout)::T
		type(SymDimension)::dimen
		dimen=SymRNDim(T%TenDim)
		T%TenDim=dimen
		T%rank=SymDimsize(dimen)
		return
	end subroutine
	integer function addressToIndes(T,Adim)
		type(SymTensor),intent(in) :: T
		integer,intent(in) :: Adim(:)
		integer,allocatable::Tdim(:)
		integer :: i,Dimlen
		Tdim=T%TenDim
		Dimlen=size(TDim)
		if(Dimlen.eq.1) then
			addressToIndes=Adim(1)
			return
		end if
		addressToIndes=0
		do i=Dimlen,2,-1
			addressToIndes=addressToIndes+(Adim(i)-1)*product(TDim(1:(i-1)))
		end do
		addressToIndes=addressToIndes+Adim(1)
		return 
	end function
	integer function addressToIndes2(TDim,Adim)
		integer,intent(in) :: TDim(:)!max dimension
		integer,intent(in) :: Adim(:)
		integer :: i,Dimlen
		Dimlen=size(TDim)
		if(Dimlen.eq.1) then
			addressToIndes2=Adim(1)
			return
		end if
		addressToIndes2=0
		do i=Dimlen,2,-1
			addressToIndes2=addressToIndes2+(Adim(i)-1)*product(TDim(1:(i-1)))
		end do
		addressToIndes2=addressToIndes2+Adim(1)
		return 
	end function	
	subroutine IndesToaddress(TDim,Adim,inde)
		integer,intent(in) :: TDim(:)!max dimension
		integer,intent(in) :: inde
		integer,allocatable,intent(inout)::Adim(:)
		integer :: i,k,lenDim,modi,j
		logical::goon
		lenDim=size(TDim)
		if(size(Adim).ne.lenDim) then
			if(allocated(Adim)) then
				deallocate(Adim)
			end if
			allocate(Adim(lenDim))
		end if
		i=1
		k=inde
		goon=.true.
		do while (goon)
			if(i.eq.lenDim) then! if[1]
				if(k.gt.Tdim(i)) then
					write(*,*)"ERROR in IndesToaddress"
					write(*,*)"stop"
					stop
				end if
				modi=mod(k,Tdim(i))
				k=k/Tdim(i)
				if(modi.eq.0) then
					Adim(i)=TDim(i)
				else
					Adim(i)=modi
				end if
				goon=.false.
			else
				modi=mod(k,Tdim(i))
				k=k/Tdim(i)
				if(modi.eq.0) then
					Adim(i)=TDim(i)
				else
					Adim(i)=modi
				end if
				if(k.eq.0)then
					goon=.false.
					do j=i+1,lenDim
						Adim(j)=1
					end do
				else
					if(modi.ne.0) then
						k=k+1
					end if
				end if
			end if! if[1]
			i=i+1
		end do
		return
	end subroutine			
!*********************  generate *********************
!		generate a Tensor with random number
!degeneracy(i,:) are the degeneracy of QN(i). QN are quantum number(max quantum number corrding to non-zero term)
!rule: symmetry rule of U(1)
	type(SymTensor) function Symgenerate_U1(QN,degeneracy,rule,region) result (T)
		real*4,intent(in) :: QN(:)
		integer,intent(in)::degeneracy(:,:),rule(:)
		real*8,optional,intent(in)::region(2)
		complex*16::minnum,maxnum,delmum
		integer :: rank,totalData,i,maxindex(size(QN)),inde(size(QN)),lenDim,address
		integer ::deg(size(QN))
		logical::goon
		real*4::QN_i(size(QN)),minQN(size(QN))
		if(present(region))then
			minnum=region(1)
			maxnum=region(2)
		else
			minnum=0d0
			maxnum=1d0
		end if
		delmum=maxnum-minnum
		lenDim=size(QN)
		T%rank=lenDim
!Dimension
		call U1QNtoMaxIndex(maxindex,QN)
		totalData=product(maxindex)
		T%totalData=totalData
		T%TenDim=maxindex
		call setDimQN(T%TenDim,QN,degeneracy,rule)
!block data		
		allocate(T%block(totalData))
		goon=.true.
		minQN=-1d0*QN
		QN_i=minQN
		T%totalblock=0
		do while(goon)
			if(U1Rule(QN_i,rule))then
				call Symindex(inde,T%TenDim,QN_i)
				address=addressToIndes2(maxindex,inde)
				call SymDimDegeneracy(deg,T%TenDim,inde)
				T%block(address)=( delmum*generate(deg) )-minnum
				
				T%totalblock=T%totalblock+1
			end if
			goon=inde_counter(QN_i,minQN,QN,1.)
		end do
		T%flag=.true.
		return
	end function
	type(SymTensor) function Symgenerate_U1_Dim(dimen,region) result (T)
		type(SymDimension),intent(in)::dimen
		real*8,optional,intent(in)::region(2)
		real*4,allocatable :: QN(:)
		integer,allocatable::degeneracy(:,:),rule(:)
		integer::lenDim,lenDeg,i
		real*4::maxQN
		lenDim=SymDimSize(dimen)
		allocate(QN(lenDim))
		allocate(rule(lenDim))
		call SymQN(QN,dimen)
		call SymRule(rule,dimen)
		maxQN=maxval(QN)
		lenDeg=2*maxQN+1
		allocate(degeneracy(lenDim,lenDeg))
		do i=1,lenDim
			call SymDimDegeneracy(degeneracy(i,:),dimen,i)
		end do
		T=Symgenerate_U1(QN,degeneracy,rule,region)
		return
	end function
!QN and Deg
!0.5:       1  1
!1  :      1  2  1
!1.5:     1  3  3  1
!2  :    1  4  6  4  1		
!2.5:   1 5 10 10  5  1 
!3  :  1 6 15 20 15 6  1
!3.5: 1 7 21 35 35 21 7 1
!These QNs or N*Deg give a good result on MPS
!input QN,Num
!output a random symTensor with degeneracy deg(i)=Num(i)*C_{N-1}^{i-1},i=0,N-1
!N is the lenof QN
	type(SymTensor) function Symgenerate_U1_QN1(QN,rule,region,Num) result (T)
		real*8,intent(in)::region(2)
		real*4,intent(in)::QN(:)
		integer,optional,intent(in)::Num(:)
		integer,intent(in)::rule(:)
		type(SymDimension)::Dimen
		Dimen=SymU1Dim(QN,rule,Num)
		T=Symgenerate_U1_Dim(Dimen)
		return
	end function
		type(SymTensor) function Symgenerate_U1_QN2(QN,rule,Num) result (T)
		real*4,intent(in)::QN(:)
		integer,intent(in)::rule(:)
		integer,optional,intent(in)::Num(:)
		type(SymDimension)::Dimen
		Dimen=SymU1Dim(QN,rule,Num)
		T=Symgenerate_U1_Dim(Dimen)
		return
	end function
!Input a Tensor,translate it to a U(1) symmetry Tensor,according to the QN(quantum number) , and rule
!if rulecheck=true, code will check if it is a U(1) symmetry Tensor
!If not, Do not check, just pick up the element of the non-zero term according to the QN(quantum number) and rule
!The degeneracy of the output are 1
	type(SymTensor) function Tensor2U1Ten(Ten,QN,rule,rulecheck_) result (T)
		type(Tensor),intent(in)::Ten
		real*4,intent(in) :: QN(:)
		integer,intent(in)::rule(:)
		logical,optional,intent(in)::rulecheck_
		integer :: rank,totalData,i,maxindex(size(QN)),inde(size(QN)),lenDim,address,blockdim(size(QN))
		logical::goon,rulecheck
		real*4::QN_i(size(QN)),minQN(size(QN))
		if(.not. present(rulecheck_)) then
			rulecheck=.false.
		else
			rulecheck=rulecheck_
		end if
		lenDim=size(QN)
		T%rank=lenDim
!Dimension
		call U1QNtoMaxIndex(maxindex,QN)
		totalData=product(maxindex)
		T%totalData=totalData
		T%TenDim=maxindex
		call setDimQN(T%TenDim,QN,rule)
!block data		
		allocate(T%block(totalData))
		goon=.true.
		minQN=-1d0*QN
		QN_i=minQN
		T%totalblock=0
		blockdim=1
		do while(goon)
			if(U1Rule(QN_i,rule))then!if[1]
				call Symindex(inde,T%TenDim,QN_i)
				address=addressToIndes2(maxindex,inde)
				T%block(address)=Ten.i.inde
				call resetdim(T%block(address),blockdim)
				T%totalblock=T%totalblock+1
			else!if[1]
				if(rulecheck) then!if[2]
					call Symindex(inde,T%TenDim,QN_i)
					if(cdabs(Ten.i.inde).gt.zero_error_)then!if[3]
						write(*,*)"It is not a symmety Tensor"
						stop
					end if!if[3]
				end if!if[2]
			end if!if[1]
			goon=inde_counter(QN_i,minQN,QN,1.)
		end do
		T%flag=.true.
		return
	end function
!Transfer the SymTensot to Tensor
	subroutine U1Tensor2Ten(T,ST)
		type(Tensor),intent(inout)::T
		type(SymTensor),intent(in)::ST
		type(Tensor)::block
		integer::address,blockinde
		logical::goon,goon2
		real*4,allocatable::QN_i(:),minQN(:),QN(:)
		integer,allocatable::imin(:),imax(:),inde(:)
		integer,allocatable::newinde(:),intdim(:),maxindex(:)
		allocate(intdim(SymgetRank(ST)))
		allocate(imin(SymgetRank(ST)))
		allocate(imax(SymgetRank(ST)))
		allocate(inde(SymgetRank(ST)))
		allocate(newinde(SymgetRank(ST)))
		allocate(QN_i(SymgetRank(ST)))
		allocate(minQN(SymgetRank(ST)))
		allocate(QN(SymgetRank(ST)))
		allocate(maxindex(SymgetRank(ST)))
		call QN2nonSymmaxdim(ST%TenDim,intdim)
		call cleanTensor(T)
		call allocateTensor(T,intdim)
		T=0d0
		goon=.true.
		call SymQN(QN,ST%TEnDim)
		call U1QNtoMaxIndex(maxindex,QN)
		minQN=-1d0*QN
		QN_i=minQN
		do while(goon)
			call Symindex(inde,ST%TenDim,QN_i)
			address=addressToIndes2(maxindex,inde)
			block=ST%block(address)
			call QN2nonSyminde(ST%TenDim,QN_i,imin,imax)
			goon2=.true.
			blockinde=1
			newinde=imin
			do while(goon2)
				if(getFlag(block))then
					call modify(T,newinde,block%Tensor_Data(blockinde))
					blockinde=blockinde+1
				end if
				goon2=inde_counter(newinde,imin,imax,1)
			end do
			goon=inde_counter(QN_i,minQN,QN,1.)
		end do
		return
	end subroutine
	subroutine SymTMprint(ST,printtype)
		type(SymTensor),intent(in) :: ST
		integer,optional,intent(in)::printtype
		type(Tensor)::T
		T=ST
		if(present(printtype)) then
			call TMprint(T,printtype)
		else
			call TMprint(T)
		end if
		return
	end subroutine
	subroutine SymTprint(ST,printtype)
		type(SymTensor),intent(in) :: ST
		integer,optional,intent(in)::printtype
		type(Tensor)::T
		T=ST
		write(*,*)"**********************"
		write(*,*)"***** START PRINT ****"
		write(*,*)"**********************"
		write(*,*)"turn the SymTensor into Tensor"
		if(present(printtype)) then
			call Tprint(T,printtype)
		else
			call Tprint(T)
		end if
		call SymDprint0(ST%TenDim)
		write(*,*)"**********************"
		write(*,*)"*****  END PRINT  ****"
		write(*,*)"**********************"
		write(*,*)" "
		write(*,*)" "
		write(*,*)" "
		return
	end subroutine
	subroutine SymBTMprint(ST,printtype)
		type(SymTensor),intent(in) :: ST
		integer,optional,intent(in)::printtype
		integer::i
		real*4,allocatable::QN_i(:)
		integer,allocatable::inde(:),maxdim(:)
		allocate(QN_i(SymgetRank(ST)))
		allocate(inde(SymgetRank(ST)))
		allocate(maxdim(SymgetRank(ST)))
		maxdim=.dim.ST
		write(*,*)"**********************"
		write(*,*)"***** START PRINT ****"
		write(*,*)"**********************"
		do i=1,ST%TotalData
			if(getFlag(ST%block(i)))then
				call IndesToaddress(maxdim,inde,i)
				call SymQN(QN_i,ST%TenDim,inde)
				write(*,*)"The block of",QN_i
				if(present(printtype)) then
					call TMprint(ST%block(i),printtype)
				else
					call TMprint(ST%block(i))
				end if
				write(*,*)"  ===   "
			end if
		end do
		call SymDprint0(ST%TenDim)
		write(*,*)"**********************"
		write(*,*)"*****  END PRINT  ****"
		write(*,*)"**********************"
		write(*,*)" "
		write(*,*)" "
		write(*,*)" "
		return
	end subroutine
	subroutine SymBTprint(ST,printtype)
		type(SymTensor),intent(in) :: ST
		integer,optional,intent(in)::printtype
		integer::i
		real*4,allocatable::QN_i(:)
		integer,allocatable::inde(:),maxdim(:)
		allocate(QN_i(SymgetRank(ST)))
		allocate(inde(SymgetRank(ST)))
		allocate(maxdim(SymgetRank(ST)))
		maxdim=.dim.ST
		write(*,*)"**********************"
		write(*,*)"***** START PRINT ****"
		write(*,*)"**********************"
		do i=1,ST%TotalData
			if(getFlag(ST%block(i)))then
				call IndesToaddress(maxdim,inde,i)
				call SymQN(QN_i,ST%TenDim,inde)
				write(*,*)"The block of",QN_i
				if(present(printtype)) then
					call Tprint(ST%block(i),printtype)
				else
					call Tprint(ST%block(i))
				end if
				write(*,*)"  ===   "
			end if
		end do
		call SymDprint0(ST%TenDim)
		write(*,*)"**********************"
		write(*,*)"*****  END PRINT  ****"
		write(*,*)"**********************"
		write(*,*)" "
		write(*,*)" "
		write(*,*)" "
		return
	end subroutine
	
	subroutine SymTDprint(ST)
		type(SymTensor) :: ST
		if(ST%flag) then
			write(*,*)"**********************"
			write(*,*)"***** START PRINT ****"
			write(*,*)"**********************"
			write(*,*) "The rank of the SymTensor is"
			write(*,*) ST%rank
			write(*,*) "The number of  data of the SymTensor is"
			write(*,*) ST%totalData
			write(*,*) "The number of  blocks of the SymTensor is"
			write(*,*) ST%totalblock
			call SymDprint0(ST%TenDim)
			write(*,*)"**********************"
			write(*,*)"*****  END PRINT  ****"
			write(*,*)"**********************"
			write(*,*)" "
			write(*,*)" "
			write(*,*)" "
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
	subroutine SymBTDprint(ST)
		type(SymTensor),intent(in) :: ST
		integer::i
		real*4,allocatable::QN_i(:)
		integer,allocatable::inde(:),maxdim(:)
		allocate(QN_i(SymgetRank(ST)))
		allocate(inde(SymgetRank(ST)))
		allocate(maxdim(SymgetRank(ST)))
		maxdim=.dim.ST
		write(*,*)"**********************"
		write(*,*)"***** START PRINT ****"
		write(*,*)"**********************"
		do i=1,ST%TotalData
			if(getFlag(ST%block(i)))then
				call IndesToaddress(maxdim,inde,i)
				call SymQN(QN_i,ST%TenDim,inde)
				write(*,*)"The block of",QN_i
				call TDprint(ST%block(i),.true.)
				write(*,*)"  ===   "
			end if
		end do
		write(*,*) "The rank of the SymTensor is"
		write(*,*) ST%rank
		write(*,*) "The number of  data of the SymTensor is"
		write(*,*) ST%totalData
		write(*,*) "The number of  blocks of the SymTensor is"
		write(*,*) ST%totalblock
		call SymDprint0(ST%TenDim)
		write(*,*)"**********************"
		write(*,*)"*****  END PRINT  ****"
		write(*,*)"**********************"
		write(*,*)" "
		write(*,*)" "
		write(*,*)" "
		return
	end subroutine
!cccccccccccccccc add          cccccccccccccccccc
	type(SymTensor) function Symadd(T1,T2)result(add)
		type(SymTensor),intent(in) :: T1,T2
		type(SymDimension)::dim1,dim2
		integer::testint,i
		if((.not.T1%flag).or.(.not.T2%flag))then
			write(*,*)"There is no data in the SymTensor,(+)"
			stop
		end if
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*)"The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		add%Totaldata=T1%totaldata
		add%Totalblock=T1%Totalblock
		add%flag=T1%flag
		add%TenDim=T1%TenDim
		add%rank=T1%rank
		allocate(add%block(add%Totaldata))
		testint=0
		do i=1,add%Totaldata
			if(getFlag(T1%block(i)).and.getFlag(T2%block(i)))then
				add%block(i)=T1%block(i)+T2%block(i)
				testint=testint+1
			end if
		end do
		!if(testint.ne.add%Totalblock) then
		!	write(*,*)"ERROR in (+)"
		!	stop
		!end if
		return
	end function
!cccccccccccccccc    minus       cccccccccccccccccc		
	type(SymTensor) function Symminus(T1,T2)result(minus)
		type(SymTensor),intent(in) :: T1,T2
		type(SymDimension)::dim1,dim2
		integer::i,testint
		if((.not.T1%flag).or.(.not.T2%flag))then
			write(*,*)"There is no data in the Tensor,(-)"
			stop
		end if
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*)"The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		minus%Totaldata=T1%totaldata
		minus%totalblock=T1%totalblock
		minus%flag=T1%flag
		minus%TenDim=T1%TenDim
		minus%rank=T1%rank
		allocate(minus%block(minus%Totaldata))
		testint=0
		do i=1,minus%Totaldata
			if(getFlag(T1%block(i)).and.getFlag(T2%block(i)))then
				minus%block(i)=T1%block(i)-T2%block(i)
				testint=testint+1
			else
				if(getFlag(T1%block(i)))then
					minus%block(i)=T1%block(i)-dcmplx(0d0,0d0)
				else if(getFlag(T2%block(i)))then
					minus%block(i)=dcmplx(0d0,0d0)-T2%block(i)
				end if
				testint=testint+1
			end if
		end do
		!if(testint.ne.minus%Totalblock) then
		!	write(*,*)"ERROR in (-)"
		!	write(*,*)testint,minus%Totalblock
		!	stop
		!end if
		return
	end function
!cccccccccccccccc divide_number          cccccccccccccccccc	
	type(SymTensor) function Symdivide_number(T1,num)result(divide_number)
		type(SymTensor),intent(in) :: T1
		complex*16,intent(in) :: num
		integer::i
		if(.not.T1%flag)then
			write(*,*)"There is no data in the Tensor,(/)"
			stop
		end if
		divide_number%rank=T1%rank
		divide_number%TenDim=T1%TenDim
		divide_number%totalData=T1%totalData
		divide_number%totalblock=T1%totalblock
		allocate(divide_number%block(divide_number%Totaldata))
		do i=1,divide_number%Totaldata
			if(getFlag(T1%Block(i)))then
				divide_number%block(i)=T1%Block(i)/num
			end if
		end do
		divide_number%flag=.true.
		return
	end function
!ccccccccccccccccc divide_number          cccccccccccccccccc	
	type(SymTensor) function Symdivide_numberReal(T1,num)result(divide_numberReal)
		type(SymTensor),intent(in) :: T1
		real*8,intent(in) :: num
		integer::i
		if(.not.T1%flag)then
			write(*,*)"There is no data in the Tensor,(/)"
			stop
		end if
		divide_numberReal%rank=T1%rank
		divide_numberReal%TenDim=T1%TenDim
		divide_numberReal%totalData=T1%totalData
		divide_numberReal%totalblock=T1%totalblock
		allocate(divide_numberReal%block(divide_numberReal%Totaldata))
		do i=1,divide_numberReal%Totaldata
			if(getFlag(T1%Block(i)))then
				divide_numberReal%block(i)=T1%Block(i)/num
			end if
		end do
		divide_numberReal%flag=.true.
		return
	end function
!cccccccccccccccc multiply_number          cccccccccccccccccc	
	type(SymTensor) function Symmultiply_number(T1,num)result(multiply_number)
		type(SymTensor),intent(in) :: T1
		complex*16,intent(in) ::   num
		integer::i
		if(.not.T1%flag)then
			write(*,*)"There is no data in the Tensor,(*)"
			stop
		end if
		multiply_number%rank=T1%rank
		multiply_number%TenDim=T1%TenDim
		multiply_number%totalData=T1%totalData
		multiply_number%totalblock=T1%totalblock
		allocate(multiply_number%Block(multiply_number%Totaldata))
		do i=1,multiply_number%Totaldata
			if(getFlag(T1%Block(i)))then
				multiply_number%block(i)=T1%Block(i)*num
			end if
		end do
		multiply_number%flag=.true.
		return
	end function
	type(SymTensor) function Symmultiply_number_(num,T1)result(multiply_number_)
		type(SymTensor),intent(in) :: T1
		complex*16,intent(in) ::   num
		integer::i
		if(.not.T1%flag)then
			write(*,*)"There is no data in the Tensor,(*)"
			stop
		end if
		multiply_number_%rank=T1%rank
		multiply_number_%TenDim=T1%TenDim
		multiply_number_%totalData=T1%totalData
		multiply_number_%totalblock=T1%totalblock
		allocate(multiply_number_%Block(multiply_number_%Totaldata))
		do i=1,multiply_number_%Totaldata
			if(getFlag(T1%Block(i)))then
				multiply_number_%block(i)=num*T1%Block(i)
			end if
		end do
		multiply_number_%flag=.true.
		return
	end function
!cccccccccccccccc multiply_real          cccccccccccccccccc	
	type(SymTensor) function Symmultiply_real(T1,num)result(multiply_real)
		type(SymTensor),intent(in) :: T1
		real*8,intent(in) ::   num
		integer::i
		if(.not.T1%flag)then
			write(*,*)"There is no data in the Tensor,(*)"
			stop
		end if
		multiply_real%rank=T1%rank
		multiply_real%TenDim=T1%TenDim
		multiply_real%totalData=T1%totalData
		multiply_real%totalblock=T1%totalblock
		allocate(multiply_real%Block(multiply_real%Totaldata))
		do i=1,multiply_real%Totaldata
			if(getFlag(T1%Block(i)))then
				multiply_real%block(i)=T1%Block(i)*num
			end if
		end do
		multiply_real%flag=.true.
		return
	end function
	type(SymTensor) function Symmultiply_real_(num,T1)result(multiply_real_)
		type(SymTensor),intent(in) :: T1
		real*8,intent(in) ::   num
		integer::i
		if(.not.T1%flag)then
			write(*,*)"There is no data in the Tensor,(*)"
			stop
		end if
		multiply_real_%rank=T1%rank
		multiply_real_%TenDim=T1%TenDim
		multiply_real_%totalData=T1%totalData
		multiply_real_%totalblock=T1%totalblock
		allocate(multiply_real_%Block(multiply_real_%Totaldata))
		do i=1,multiply_real_%Totaldata
			if(getFlag(T1%Block(i)))then
				multiply_real_%block(i)=num*T1%Block(i)
			end if
		end do
		multiply_real_%flag=.true.
		return
	end function	
	subroutine productSymTest(indim1,indim2)
		type(SymDimension),intent(in)::indim1,indim2
		type(SymDimension)::dim1,dim2
		integer::i,dimsize
		integer,allocatable::rule1(:),rule2(:)
		dim1=SymDimDecomposeAll(indim1)
		dim2=SymDimDecomposeAll(indim2)
		dimsize=SymDimSize(dim1)
		if(dimsize.ne.SymDimSize(dim2))then
			write(*,*)"Cannot do the product"
			write(*,*)"dimension"
			write(*,*)dimsize,SymDimSize(dim2)
			stop
		end if
		allocate(rule1(dimsize))
		allocate(rule2(dimsize))
		call SymRule(rule1,dim1)
		call SymRule(rule2,dim2)
		do i=1,dimsize
			if((rule1(i)*rule2(i)).ne.-1)then
				write(*,*)"Cannot do the product"
				write(*,*)"The symmetry rule"
				write(*,*)rule1
				write(*,*)rule2
				stop
			end if
		end do
		return
	end subroutine 
	subroutine SymDiagonalTest(dim1,dim2)
		type(SymDimension),intent(in)::dim1,dim2
		integer::N1,N2
		N1=SymDimSize(dim1)
		N2=SymDimSize(dim2)
		if(N1.ne.1 .or. N2.ne.1 )then
			write(*,*) "it is not a matrix"
			stop
		end if
		call SymRule(N1,dim1,1)
		call SymRule(N2,dim2,1)
		if((N1*N2).ne.-1 )then
			write(*,*) "it is not a Diagonal matrix"
			stop
		end if
		return
	end subroutine 
	
	type(SymTensor) function SymProductTensor(ST1,ST2) result(Res)
		type(SymTensor),intent(in) :: ST1,ST2
		integer::rank1,rank2,flag,T1m,T1n,T2m,T2n,T1l,T2l
		type(SymDimension)::D1,D2,newD,D11,D22
		integer::i,totaldata,totalblock,row(3),col(3)
		real*4::QN1,QN2,minQN
		complex*16::NData
		if(.not.ST1%flag) then
			write(*,*)"ERROR in (*)"
			write(*,*)"there is no data is the first SymTensor"
			write(*,*)"stop"
			stop
		end if
		if(.not.ST2%flag) then
			write(*,*)"ERROR in (*)"
			write(*,*)"there is no data is the second SymTensor"
			write(*,*)"stop"
			stop
		end if	
		rank1=SymgetRank(ST1)
		rank2=SymgetRank(ST2)
		D1=.subDim.ST1
		D2=.subDim.ST2
		row=1
		col=1
		if((rank1.eq.1).and.(rank2.eq.1)) then
			flag=1
		else if((rank1.eq.1).and.(rank2.ge.2)) then
			flag=2
		else if((rank1.ge.2).and.(rank2.eq.1)) then
			flag=3	
		else if((rank1.ge.2).and.(rank2.ge.2)) then
			flag=4
		else
			write(*,*)"ERROR in ProductTensor",rank1,rank2
			stop
		end if
		select case (flag)
			case (1)
				T1m=1
				T1n=ST1.dim.1
				T2m=ST2.dim.1
				T2n=1
				if(T1m.ne.T2n) then
					write(*,*)"Do no finish this part"
					stop
				end if
				call productSymTest(D1,D2)
				allocate(Res%block(1))
				NData=Szdotu(T1n,ST1%BLOCK,1,ST2%BLOCK,1)
				Res%BLOCK(1)=(/NData/)
				Res%totalData=1
				Res%totalblock=1
				Res%TenDim=SymU1Dim(0.,(/Res%BLOCK(1).dim.1/),0)
				Res%rank=1
				Res%flag=.true.
				!Ndata=ZDOTU(T1m,T1%Tensor_Data,1,T2%Tensor_Data,1)
				!ProductTensor=set((/1/),(/Ndata/))
				return
			case(2)
				if(rank2.ge.2) then
					newD=D2.sub.2
					do i=3,rank2
						newD=newD+(D2.sub.i)
					end do
					D2=SymDimConstract(D2,2,rank2)
				end if
				T1m=1
				T1n=D1.i.1
				T2m=D2.i.1
				T2n=D2.i.2
				if((D1.i.1) .ne. T2m) then
					write(*,*)"Do no finish this part"
					stop
				end if
				call productSymTest(D1,D2.sub.1)
				allocate(Res%block(T2n))
				call SZGEMM(row,col,T1m,T2n,T1n,ST1%BLOCK,T1m,T1n,&
					ST2%BLOCK,T2m,T2n,Res%BLOCK,T1m,T2n,totalblock)
				Res%totalData=T2n
				Res%totalblock=totalblock
				Res%TenDim=newD
				Res%rank=SymDimSize(newD)
				Res%flag=.true.
				return
			case (3)
				if(rank1.ge.2) then
					newD=D1.sub.1
					do i=2,rank1-1
						newD=newD+(D1.sub.i)
					end do
					D1=SymDimConstract(D1,1,rank1-2)
				end if
				T1m=D1.i.1
				T1n=D1.i.2
				T2m=D2.i.1
				T2n=1
				if((D2.i.1) .ne. T1n) then
					write(*,*)"Do no finish this part"
					stop
				end if
				call productSymTest(D1.sub.2,D2)
				allocate(REs%block(T1m))
				call SZGEMM(row,col,T1m,T2n,T1n,ST1%BLOCK,T1m,T1n,&
					ST2%BLOCK,T2m,T2n,Res%BLOCK,T1m,T2n,totalblock)
				Res%totalData=T1m
				Res%totalblock=totalblock
				Res%TenDim=newD
				Res%rank=SymDimSize(newD)
				Res%flag=.true.
				return
			case(4)
				if(rank1.ge.2) then
					newD=D1.sub.1
					do i=2,rank1-1
						newD=newD+(D1.sub.i)
					end do
					D1=SymDimConstract(D1,1,rank1-2)
				end if
				if(rank2.ge.2) then
					do i=2,rank2
						newD=newD+(D2.sub.i)
					end do
					D2=SymDimConstract(D2,2,rank2)
				end if
				call productSymTest(D1.sub.2,D2.sub.1)
				if((D1.i.2).ne.(D2.i.1)) then!if1
					D11=D1.sub.2
					D22=D2.sub.1
					if(.not.Symif_original_dim(D11) .or. .not.Symif_original_dim(D22)) then
						write(*,*)"ERROR in SymProductTensor_dif_dim"
						write(*,*)"The last index ST1 and the fist index of ST2 should be original dim"
						stop
					end if
					QN1=SymQNum(D11,1)
					QN2=SymQNum(D22,1)
					if(QN1.le.QN2)then!if2
						minQN=SymQNum(D11,1,1)
						Row(2)=Symindexi(D22,1,minQN)
					else
						minQN=SymQNum(D22,1,1)
						col(1)=Symindexi(D11,1,minQN)
					end if!if2
				end if!if1
				T1m=D1.i.1
				T1n=D1.i.2
				T2m=D2.i.1
				T2n=D2.i.2
				totaldata=T1m*T2n
				allocate(Res%block(totalData))
				!(Row,Col,M,N,K,A,LDRA,LDCA,B,LDRB,LDCB,C,LDRC,LDCC,totalblock)
				call SZGEMM(row,col,T1m,T2n,T1n,ST1%BLOCK,T1m,T1n,&
					ST2%BLOCK,T2m,T2n,Res%BLOCK,T1m,T2n,totalblock)
				Res%totalData=totaldata
				Res%totalblock=totalblock
				Res%TenDim=newD
				Res%rank=SymDimSize(newD)
				Res%flag=.true.
		case default 
			write(*,*) "ERROR in ProductTensor,no such data"
			stop
		end 	select
		return
	end function
	
		
!  Row(3):start Row of A,B,C, Col(3):start Col of A,B,C
!    Row,Col determine which part of element will be use
!	M the munber of row of A to use for product
!  N the munber of Col of B to use for product
!  K the munber of col of A of Row of B to use for product
!  A, B, C the matrixes
!	LDRA,LDCA the size of A
!	LDRB,LDCB the size of B
!	LDRC,LDCC the size of C
!	totalblock: in output the total block of C
!  design for the U(1) symmetry product of
!
!  s=2.5     s=0.5      s=1.5   s=0.5        2.5     0.5
!   -->--[A]-->--   *   -->--[B]-->--   =   -->--[C]-->-- 
!
!  A is a matrix of LDRA,LDCA
!  B is a matrix of LDRB,LDCB
!  B is a matrix of LDRC,LDCC
!  A(Row(1):Row(1)+M,col(1):col(1)+K) * B( Row(2):Row(2)+K,col(2):col(2)+N)
!     =C(Row(3):Row(3)+M,col(3):col(3)+N)
!
	subroutine SZGEMM(Row,Col,M,N,K,A,LDRA,LDCA,B,LDRB,LDCB,C,LDRC,LDCC,totalblock)!ZGEMM for block data,Not finished but it could run
	  	INTEGER,intent(in)::M,N,k,LDRA,LDCA,LDRB,LDCB,LDRC,LDCC
	  	integer,intent(in)::Row(3),Col(3)
	  	type(Tensor),intent(in)::A(ldra,*),B(ldrb,*)
	  	type(Tensor),intent(inout)::C(ldrc,*)
	  	integer,intent(inout)::totalblock
	  	integer::i,j,l,rowA,colA,rowB,colB,rowC,colC
	  	integer::RA,CA,RB,CB,RC,CC
	  	type(Tensor)::zero
	  	totalblock=0
	  	rowA=Row(1)-1
	  	rowB=Row(2)-1
	  	rowC=Row(3)-1
	  	colA=Col(1)-1
	  	colB=Col(2)-1
	  	colC=Col(3)-1
	  	do i=1,M!do2
  			RA=rowA+i
  			RC=rowC+i
		  	do j=1,N!do1
		  		CB=colB+j
		  		CC=ColC+j
	  			do l=1,k!do3
	  				CA=colA+l
	  				RB=rowB+l
	  				if(getFlag(A(RA,CA)).and.getFlag(B(RB,CB)))then
	  					if(getflag(C(RC,CC)))then
	  						C(RC,CC)=C(RC,CC)+A(RA,CA)*B(RB,CB)
	  					else
	  						C(RC,CC)=A(RA,CA)*B(RB,CB)
	  					end if
	  				end if
	  			end do!do3
	  			if(getflag(C(RC,CC))) totalblock=totalblock+1
	  		end do!do2
	  	end do!do1
	end subroutine
	COMPLEX*16 FUNCTION Szdotu(N,ZX,INCX,ZY,INCY)!zdotu for block data
		integer,intent(in)::INCX,INCY,N
		type(Tensor),intent(in)::ZX(*),ZY(*)
		integer::i,ix,iy
		Szdotu=(0.,0.)
		if(n.le.0) return
		if(incx.eq.1 .and. incy.eq.1) then
			do i=1,n
				if(getflag(zx(i)).and.getFlag(zy(i))) then
					Szdotu = Szdotu + ( zx(i) .numx. zy(i) )
				end if
			end do
		else
			ix=1
			iy=1
			IF (incx.LT.0) ix = (-n+1)*incx + 1
			IF (incy.LT.0) iy = (-n+1)*incy + 1
			DO i = 1,n
				if(getflag(zx(ix)).and.getFlag(zy(iy))) then
					Szdotu = Szdotu + (zx(ix) .numx. zy(iy))
				end if
				ix = ix + incx
				iy = iy + incy
			end do
		end if
		return
	end FUNCTION
		COMPLEX*16 FUNCTION SzdotC(N,ZX,INCX,ZY,INCY)!zdotc for block data,! dot product conjugating the first vector,The Tensor will be regard as a vector
		integer,intent(in)::INCX,INCY,N
		type(Tensor),intent(in)::ZX(*),ZY(*)
		integer::i,ix,iy
		SzdotC=(0.,0.)
		if(n.le.0) return
		if(incx.eq.1 .and. incy.eq.1) then
			do i=1,n
				if(getflag(zx(i)).and.getFlag(zy(i))) then
					SzdotC = SzdotC + ( zx(i) .x. zy(i) )
				end if
			end do
		else
			ix=1
			iy=1
			IF (incx.LT.0) ix = (-n+1)*incx + 1
			IF (incy.LT.0) iy = (-n+1)*incy + 1
			DO i = 1,n
				if(getflag(zx(ix)).and.getFlag(zy(iy))) then
					SzdotC = SzdotC + (zx(ix) .x. zy(iy))
				end if
				ix = ix + incx
				iy = iy + incy
			end do
		end if
		return
	end FUNCTION
!ccccccccccccccccc permute_rank3   cccccccccccccccccc		
	type(SymTensor) function Sympermute_rank3(T1,index_not_permute)result(permute_rank3)
		type(SymTensor),intent(in) :: T1
		integer,intent(in) ::   index_not_permute
		integer ::i,newdimen(3),totaldata
		integer :: dimen(3)
		type(SymDimension) ::newTDim
		type(Tensor),allocatable :: Tdata(:,:,:),newdata(:,:,:)
		if(T1%rank.ne.3) then
			write(*,*)"ERROR in permute_rank3"
			write(*,*)"stop"
			stop
			return
		end if
		dimen=T1%TenDim
		allocate(Tdata(dimen(1),dimen(2),dimen(3)))
		Tdata=T1
		if(index_not_permute.eq.1) then
			newdimen(1)=dimen(1)
			newdimen(2)=dimen(3)
			newdimen(3)=dimen(2)
			newTDim=SymDimpermute(T1%TenDim,(/1,3,2/))
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(1)
				newdata(i,:,:)=transpose(Tdata(i,:,:))
			end do
			call SymstoreTenData(permute_rank3,newdata,T1%totalData)
			permute_rank3%rank=3
			permute_rank3%totalData=T1%totalData
			permute_rank3%totalblock=T1%totalblock
			permute_rank3%TenDim=newTDim
			permute_rank3%flag=.true.
		end if
		if(index_not_permute.eq.2) then
			newdimen(1)=dimen(3)
			newdimen(2)=dimen(2)
			newdimen(3)=dimen(1)
			newTDim=SymDimpermute(T1%TenDim,(/3,2,1/))
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(2)
				newdata(:,i,:)=transpose(Tdata(:,i,:))
			end do
			call SymstoreTenData(permute_rank3,newdata,T1%totalData)
			permute_rank3%rank=3
			permute_rank3%totalData=T1%totalData
			permute_rank3%TenDim=newTDim
			permute_rank3%flag=.true.
		end if
		if(index_not_permute.eq.3) then
			newdimen(1)=dimen(2)
			newdimen(2)=dimen(1)
			newdimen(3)=dimen(3)
			newTDim=SymDimpermute(T1%TenDim,(/2,1,3/))
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(3)
				newdata(:,:,i)=transpose(Tdata(:,:,i))
			end do
			call SymstoreTenData(permute_rank3,newdata,T1%totalData)
			permute_rank3%rank=3
			permute_rank3%totalData=T1%totalData
			permute_rank3%TenDim=newTDim
			permute_rank3%flag=.true.
		end if
		do i=1,permute_rank3%totalData
			if(getflag(permute_rank3%block(i))) then
				permute_rank3%block(i)=permute_rank3%block(i).p.index_not_permute
			end if
		end do
		return
	 end function

!cccccccccccccccc permute_rank2   cccccccccccccccccc			 
	type(SymTensor) function Sympermute_rank2 (T)result(permute_rank2)
		type(SymTensor),intent(in) :: T
		type(Tensor),allocatable :: Tdata(:,:),newdata(:,:)
		integer ::dimen(2),i
		type(SymDimension) ::newTDim
		if(T%rank.eq.1) then
			permute_rank2=T
			return
		end if
		if(T%rank.gt.2) then
			write(*,*)"ERROR in Dpermute_rank2"
			write(*,*)"stop"
			stop
			return
		end if
		dimen=T%TenDim
		newTDim=SymDimpermute(T%TenDim,(/2,1/))
		allocate(Tdata(dimen(1),dimen(2)))
		Tdata= T
		allocate(newdata(dimen(2),dimen(1)))
		newdata=transpose(Tdata)
		call SymstoreTenData(permute_rank2,newdata,T%totalData)
		permute_rank2%rank=2
		permute_rank2%totalData=T%totalData
		permute_rank2%TenDim=newTDim
		permute_rank2%totalblock=T%totalblock
		permute_rank2%flag=.true.
		do i=1,permute_rank2%totalData
			if(getflag(permute_rank2%block(i))) then
				permute_rank2%block(i)=.p.permute_rank2%block(i)
			end if
		end do
		return
	end function
	subroutine Sympermutation_data3(T1data,index_not_permute,dimens,totaldata)
		Type(Tensor),intent(inout) :: T1data(:)
		integer,intent(in) ::   index_not_permute
		type(symDimension),intent(inout) ::dimens
		integer::dimen(3)
		integer ::i,newdimen(3),totaldata
		Type(Tensor),allocatable :: newdata(:,:,:),Tdata(:,:,:)
		dimen=dimens
		allocate(Tdata(dimen(1),dimen(2),dimen(3)))
		call SZCOPY(totaldata,T1data,1,Tdata,1)
		if(index_not_permute.eq.1) then
			newdimen(1)=dimen(1)
			newdimen(2)=dimen(3)
			newdimen(3)=dimen(2)
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(1)
				newdata(i,:,:)=transpose(Tdata(i,:,:))
			end do
			call SZCOPY(totaldata,newdata,1,T1data,1)
			dimens=symDimpermute(dimens,(/1,3,2/))
		end if
		if(index_not_permute.eq.2) then
			newdimen(1)=dimen(3)
			newdimen(2)=dimen(2)
			newdimen(3)=dimen(1)
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(2)
				newdata(:,i,:)=transpose(Tdata(:,i,:))
			end do
			call SZCOPY(totaldata,newdata,1,T1data,1)
			dimens=SymDimpermute(dimens,(/3,2,1/))
		end if
		if(index_not_permute.eq.3) then
			newdimen(1)=dimen(2)
			newdimen(2)=dimen(1)
			newdimen(3)=dimen(3)
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(3)
				newdata(:,:,i)=transpose(Tdata(:,:,i))
			end do
			call SZCOPY(totaldata,newdata,1,T1data,1)
			dimens=SymDimpermute(dimens,(/2,1,3/))
		end if
		do i=1,totalData
			if(getflag(T1data(i))) then
				T1data(i)=T1data(i).p.index_not_permute
			end if
		end do
		return
	 end subroutine
	subroutine Sympermutation_data2 (T1data,dimens,totaldata)
		Type(Tensor),intent(inout) :: T1data(:)
		type(SymDimension),intent(inout) ::dimens
		integer,intent(inout)::totaldata
		Type(Tensor),allocatable::Tdata(:,:),newdata(:,:)
		integer::dimen(2),i
		if(SymDimsize(dimens).eq.1) then
			return
		end if
		dimen=dimens
		allocate(Tdata(dimen(1),dimen(2)))
		call SZCOPY(totaldata,T1data,1,Tdata,1)
		allocate(newdata(dimen(2),dimen(1)))
		newdata=transpose(Tdata)
		call SZCOPY(totaldata,newdata,1,T1data,1)
		dimens=symDimpermute(dimens,(/2,1/))
		do i=1,totalData
			if(getflag(T1data(i))) then
				T1data(i)=.p.T1data(i)
			end if
		end do
		return
	end subroutine
!cccccccccccccccc permutation   cccccccccccccccccc
	type(SymTensor) function Sympermutation(T,newOrder)result(permutation)
		type(SymTensor),intent(in) :: T
		integer,intent(in)::newOrder(:)
		integer,allocatable ::inde(:)
		integer::lenOrder,i,totaldata,j
		type(SymDimension)::dimen		
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(permutation)"
			stop
		end if
		lenorder=size(newOrder)-1
		allocate(inde(lenorder))
		totaldata=T%totalData
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		allocate(permutation%block(totaldata))
		call SZCOPY(totaldata,T%block,1,permutation%block,1)
		dimen=T%TenDim
		do i=lenorder,1,-1
			call Sympermutefo_data(permutation%block,inde(i),dimen,totaldata)
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		permutation%totalData=totaldata
		permutation%totalBlock=T%totalBlock
		call Symresetdim(permutation,dimen)
		permutation%flag=.true.
		do i=1,totalData
			if(getflag(permutation%block(i))) then
				permutation%block(i)=permutation%block(i).p.newOrder
			end if
		end do
		return
	end function
	type(SymTensor) function Sympermutation_name(T,newOrderchar)result(permutation)
		type(SymTensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::newOrderchar(:)
		integer,allocatable::newOrder(:)
		integer,allocatable ::inde(:)
		integer::lenOrder,i,totaldata,j
		type(SymDimension)::dimen		
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(permutation)"
			stop
		end if
		allocate(newOrder(size(newOrderchar)))
		newOrder=SymNameorder(T%TenDim,newOrderchar)
		lenorder=size(newOrder)-1
		allocate(inde(lenorder))
		totaldata=T%totalData
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		allocate(permutation%block(totaldata))
		call SZCOPY(totaldata,T%block,1,permutation%block,1)
		dimen=T%TenDim
		do i=lenorder,1,-1
			call Sympermutefo_data(permutation%block,inde(i),dimen,totaldata)
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		permutation%totalData=totaldata
		call Symresetdim(permutation,dimen)
		permutation%totalBlock=T%totalBlock
		permutation%flag=.true.
		do i=1,totalData
			if(getflag(permutation%block(i))) then
				permutation%block(i)=permutation%block(i).p.newOrder
			end if
		end do
		return
	end function
	subroutine sympermutefo_data(Tdata,inde,dimen,totaldata)
		Type(Tensor),intent(inout)::Tdata(:)
		type(SymDimension),intent(inout) ::dimen
		integer,intent(inout)::totaldata
		integer,intent(in)::inde
		integer::rank,num,i
		integer::oper(2,3)
		rank=SymDimsize(dimen)
		if(inde.gt.rank) then
			write(*,*)"ERROR IN permutefo"
			stop
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permutefo_data"
			write(*,*)"index=",inde
			stop
		end if
		if(inde.eq.1) then
			return
		end if
		if(inde.eq.rank) then
			dimen=SymDimConstract(dimen,1,rank-2)
			call dimoperation(TData,(/1,1,rank-1/))
			call Sympermutation_data2(Tdata,dimen,totaldata)
			dimen=SymDimDecomposeAll(dimen)
			call dimoperation(TData,(/3/))
			return
		end if
		num=inde-2
		dimen=SymDimConstract(dimen,1,num)
		num=rank-3
		dimen=SymDimConstract(dimen,3,num)
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,3,rank/)
		call dimoperation(TData,oper)
		call Sympermutation_data3(Tdata,3,dimen,totaldata)
		dimen=SymDimDecomposeAll(dimen)
		call dimoperation(TData,(/3/))
		return
	end subroutine
	subroutine Sympermuteback_data(Tdata,inde,dimen,totaldata)
		Type(Tensor),intent(inout)::Tdata(:)
		type(SymDimension),intent(inout) ::dimen
		integer,intent(inout)::totaldata
		integer,intent(in)::inde
		integer::rank,num
		integer::oper(2,3)
		rank=SymDimsize(dimen)
		if(inde.eq.rank) then
			return
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permuteback_data"
			write(*,*)"index=",inde
			stop
		end if
		if(inde.eq.1) then
			dimen=SymDimConstract(dimen,2,rank-2)
			call dimoperation(TData,(/1,2,rank/))
			call Sympermutation_data2(Tdata,dimen,totaldata)
			dimen=SymDimDecomposeAll(dimen)
			call dimoperation(TData,(/3/))
			return
		end if
		num=inde-2
		dimen=SymDimConstract(dimen,1,num)
		num=rank-3
		dimen=SymDimConstract(dimen,3,num)
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,3,rank/)
		call dimoperation(TData,oper)
		call Sympermutation_data3(Tdata,1,dimen,totaldata)
		dimen=SymDimDecomposeAll(dimen)
		call dimoperation(TData,(/3/))
		return
	end subroutine
!*****************permutefo******************************
!		T_{1,2,3,..,i,..,n},permutefo(T,i)=_{i,1,2,3,..,i-1,i+1,..,n}
!
	type(SymTensor) function Sympermutefo(T,inde)result(permutefo)
		type(SymTensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		integer::oper(2,3)
		rank=SymgetRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permutefo"
			write(*,*)inde,rank
			write(*,*)"stop"
			stop
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permutefo"
			write(*,*)"index=",inde
			stop
		end if
		if(inde.eq.1) then
			permutefo=.dc.T
			return
		end if
		if(inde.eq.rank) then
			permutefo=Symcontracts(T,1,rank-1)
			permutefo=.p.permutefo
			permutefo%TenDim=SymDimDecomposeAll(permutefo%TenDim)
			permutefo%rank=SymDimSize(permutefo%TenDim)
			return
		end if
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,3,rank/)
		permutefo=T.cd.oper
		permutefo=permutefo.p.3
		permutefo%TenDim=SymDimDecomposeAll(permutefo%TenDim)
		permutefo%rank=SymDimSize(permutefo%TenDim)
		return
	end function
	type(SymTensor) function Sympermutefo_name(T,indechar)result(permutefo)
		type(SymTensor),intent(in)::T
		character(len=*),intent(in)::indechar
		integer::inde
		integer::rank
		integer::oper(2,3)
		rank=SymgetRank(T)
		inde=SymNameorder(T%TenDim,indechar)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permutefo"
			write(*,*)inde,rank
			write(*,*)"stop"
			stop
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permutefo"
			write(*,*)"index=",inde
			stop
		end if
		if(inde.eq.1) then
			permutefo=.dc.T
			return
		end if
		if(inde.eq.rank) then
			permutefo=Symcontracts(T,1,rank-1)
			permutefo=.p.permutefo
			permutefo%TenDim=SymDimDecomposeAll(permutefo%TenDim)
			permutefo%rank=SymDimSize(permutefo%TenDim)
			return
		end if
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,3,rank/)
		permutefo=T.cd.oper
		permutefo=permutefo.p.3
		permutefo%TenDim=SymDimDecomposeAll(permutefo%TenDim)
		permutefo%rank=SymDimSize(permutefo%TenDim)
		return
	end function
!*****************permutefo******************************
!		T_{1,2,3,..,j,..,i,.,k,...,n},permutefo(T,(/i,j,k/))=_{i,j,k,1,2,3,...,n}
!
	type(SymTensor) function Sympermutefo_vec(T,vec_)result(permutefo_vec)
		type(SymTensor),intent(in)::T
		integer,intent(in)::vec_(:)
		integer::rank,lenVec,totalData,i,j
		integer::oper(2,3)
		integer,allocatable::vec(:)
		type(SymDimension)::dimen		
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(permutefo_vec)"
			stop
		end if
		rank=SymgetRank(T)
		lenVec=size(vec_)
		totalData=T%totalData
		dimen=T%TenDim
		allocate(vec(lenVec))
		vec=vec_
		allocate(permutefo_vec%block(totalData))
		call SZCOPY(totaldata,T%block,1,permutefo_vec%block,1)
		do i=lenVec,1,-1
			call Sympermutefo_data(permutefo_vec%block,vec(i),dimen,totaldata)
			do j=1,i-1
				if(vec(j).lt.vec(i))then
					vec(j)=vec(j)+1
				end if
			end do
		end do
		permutefo_vec%totalData=totaldata
		call Symresetdim(permutefo_vec,dimen)
		permutefo_vec%flag=.true.
		permutefo_vec%totalblock=T%totalblock
		return
	end function
	type(SymTensor) function Sympermutefo_vec_name(T,indechar)result(permutefo_vec)
		type(SymTensor),intent(in)::T
		character(len=*),intent(in)::indechar(:)
		integer::rank,lenVec,totalData,i,j
		integer::oper(2,3)
		integer,allocatable::vec(:)
		type(SymDimension)::dimen			
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(permutefo_vec)"
			stop
		end if
		rank=SymgetRank(T)
		lenVec=size(indechar)
		totalData=T%totalData
		dimen=T%TenDim
		allocate(vec(lenVec))
		vec=SymNameorder(T%TenDim,indechar)
		allocate(permutefo_vec%block(totalData))
		call SZCOPY(totaldata,T%block,1,permutefo_vec%block,1)
		do i=lenVec,1,-1
			call Sympermutefo_data(permutefo_vec%block,vec(i),dimen,totaldata)
			do j=1,i-1
				if(vec(j).lt.vec(i))then
					vec(j)=vec(j)+1
				end if
			end do
		end do
		permutefo_vec%totalData=totaldata
		call Symresetdim(permutefo_vec,dimen)
		permutefo_vec%flag=.true.
		permutefo_vec%totalblock=T%totalblock
		return
	end function
!*****************permuteback******************************
!		T_{1,2,3,..,i,..,n},permuteback(T,i)=_{1,2,3,..,i-1,i+1,..,n,i}
!
	type(SymTensor) function Sympermuteback(T,inde)result(permuteback)
		type(SymTensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		integer::oper(2,3)
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(permuteback)"
			stop
		end if
		rank=SymgetRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permuteback"
			write(*,*)"stop"
			stop
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permuteback,Can not find the name"
			write(*,*)"index",inde
			stop
		end if
		if(inde.eq.rank) then
			permuteback=.dc.T
			return
		end if
		if(inde.eq.1) then
			permuteback=Symcontracts(T,2,rank)
			permuteback=.p.permuteback
			call Symdimoperation(permuteback,(/3/))
			return
		end if
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,3,rank/)
		permuteback=T.cd.oper
		permuteback=permuteback.p.1
		call Symdimoperation(permuteback,(/3/))
		return
	end function
	type(SymTensor) function Sympermuteback_name(T,indechar)result(permuteback)
		type(SymTensor),intent(in)::T
		character(len=*),intent(in)::indechar
		integer::inde
		integer::rank
		integer::oper(2,3)
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(permuteback)"
			stop
		end if
		rank=SymgetRank(T)
		inde=SymNameorder(T%TenDim,indechar)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permuteback"
			write(*,*)"stop"
			stop
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permuteback,Can not find the name"
			write(*,*)"index",inde
			stop
		end if
		if(inde.eq.rank) then
			permuteback=.dc.T
			return
		end if
		if(inde.eq.1) then
			permuteback=Symcontracts(T,2,rank)
			permuteback=.p.permuteback
			call Symdimoperation(permuteback,(/3/))
			return
		end if
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,3,rank/)
		permuteback=T.cd.oper
		permuteback=permuteback.p.1
		call Symdimoperation(permuteback,(/3/))
		return
	end function
!*****************permuteback******************************
!		T_{1,2,3,..,i,..,n},permuteback(T,i)=_{1,2,3,..,i-1,i+1,..,n,i}
!
	type(SymTensor) function Sympermuteback_vec(T,vec_)result(permuteback_vec)
		type(SymTensor),intent(in)::T
		integer,intent(in)::vec_(:)
		integer::rank,lenVec,totalData,i,j
		integer::oper(2,3)
		integer,allocatable::vec(:)
		type(SymDimension)::dimen		
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(permuteback_vec)"
			stop
		end if
		rank=SymgetRank(T)
		lenVec=size(vec_)
		totalData=T%totalData
		dimen=T%TenDim
		allocate(vec(lenVec))
		vec=vec_
		allocate(permuteback_vec%block(totalData))
		call SZCOPY(totaldata,T%block,1,permuteback_vec%block,1)
		do i=1,lenVec
			call Sympermuteback_data(permuteback_vec%block,vec(i),dimen,totaldata)
			do j=lenVec,i+1,-1
				if(vec(j).gt.vec(i))then
					vec(j)=vec(j)-1
				end if
			end do
		end do
		permuteback_vec%totalData=totaldata
		call symresetdim(permuteback_vec,dimen)
		permuteback_vec%flag=.true.
		permuteback_vec%totalblock=T%totalblock
		return
	end function
	type(SymTensor) function Sympermuteback_vec_name(T,indechar)result(permuteback_vec)
		type(SymTensor),intent(in)::T
		character(len=*),intent(in)::indechar(:)
		integer::rank,lenVec,totalData,i,j
		integer::oper(2,3)
		integer,allocatable::vec(:)
		type(SymDimension)::dimen		
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(permuteback_vec)"
			stop
		end if
		rank=SymgetRank(T)
		lenVec=size(indechar)
		totalData=T%totalData
		dimen=T%TenDim
		allocate(vec(lenVec))
		vec=SymNameorder(T%TenDim,indechar)
		allocate(permuteback_vec%block(totalData))
		call SZCOPY(totaldata,T%block,1,permuteback_vec%block,1)
		do i=1,lenVec
			call Sympermuteback_data(permuteback_vec%block,vec(i),dimen,totaldata)
			do j=lenVec,i+1,-1
				if(vec(j).gt.vec(i))then
					vec(j)=vec(j)-1
				end if
			end do
		end do
		permuteback_vec%totalData=totaldata
		call symresetdim(permuteback_vec,dimen)
		permuteback_vec%flag=.true.
		permuteback_vec%totalblock=T%totalblock
		return
	end function
!****************** permuteInde**************************
!		T_{1,2,3,..,i,..,n},permuteInde(T,i)=_{2,3,..,i,1,i+1,..,n}
!
	type(symTensor) function sympermuteInde(T,inde)result(permuteInde)
		type(symTensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		integer::oper(2,3)
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(permuteInde)"
			stop
		end if
		rank=symgetRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permuteInde",inde,rank
			write(*,*)"stop"
			stop
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permuteInde"
			write(*,*)"index",inde
			stop
		end if
		if(inde.eq.1) then
			permuteInde=.dc.T
			return
		end if
		if(inde.eq.rank) then
			permuteInde=symcontracts(T,2,rank)
			permuteInde=.p.permuteInde
			permuteInde%TenDim=symDimDecomposeAll(permuteInde%TenDim)
			permuteInde%rank=symDimSize(permuteInde%TenDim)
			return
		end if
		oper(1,:)=(/1,2,inde/)
		oper(2,:)=(/1,3,rank/)
		permuteInde=T.cd.oper
		permuteInde=permuteInde.p.3
		permuteInde%TenDim=symDimDecomposeAll(permuteInde%TenDim)
		permuteInde%rank=symDimSize(permuteInde%TenDim)
		return
	end function
!****************** permutebackInde**************************
!		T_{1,2,3,..,i,..,n},permutebackInde(T,i)=_{1,2,3,..,i-1,n,i,i+1,..,n-1}
!
	type(symTensor) function sympermutebackInde(T,inde)result(permutebackInde)
		type(symTensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		integer::oper(2,3)
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(permutebackInde)"
			stop
		end if
		rank=symgetRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permutebackInde",inde,rank
			write(*,*)"stop"
			stop
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permutebackInde"
			write(*,*)"index",inde
			stop
		end if
		if(inde.eq.rank) then
			permutebackInde=.dc.T
			return
		end if
		if(inde.eq.1) then
			permutebackInde=symcontracts(T,1,rank-1)
			permutebackInde=.p.permutebackInde
			call symdimoperation(permutebackInde,(/3/))
			return
		end if
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,2,rank-inde+1/)
		permutebackInde=T.cd.oper
		permutebackInde=permutebackInde.p.1
		call symdimoperation(permutebackInde,(/3/))
		return
	end function
!******************  Tenproduct  *********************
!	T1:[i1,i2,i3,i4,i5,i6,i7,i8]
!	T2:[j1,j2,j3,j4,j5,j6,j7,j8,j9,j10]
!	i1=(/5,1,2/)
!	i2=(10,3,5,6)
!	then the result will be T1'*T2'
!	where ,
!	T1'=[i3,i4,i6,i7,i8,(i5*i1*i2)]
!	T2'=[(j10*j3*j5*j6),j1,j2,j4,j7,j8,j9]
!	output T1'*T2
! 	input Tensor should be in its original dimenison,there is no contract on it
	type(symTensor) function symTenproduct_noName(T1_,i1,T2_,i2) result(T)
		type(symTensor),intent(in) :: T1_,T2_
		integer,intent(in) :: i1(:),i2(:)
		type(symTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.(symif_original_dim(T1_%TenDim).and.symif_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in symTenproduct"
			write(*,*)"stop"
			stop
		end if
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		call symdimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		call symdimoperation(T2,(/1,1,leni2/))
		
		T=T1 * T2
		return
	end function
	type(SymTensor) function SymTenproduct_name(T1_,name1,T2_,name2) result(T)
		type(SymTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:),name2(:)
		integer :: i1(size(name1)),i2(size(name2))
		type(SymTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.(Symif_original_dim(T1_%TenDim).and.Symif_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		i1=SymNameorder(T1_%TenDim,name1)
		i2=SymNameorder(T2_%TenDim,name2)
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		call Symdimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		call Symdimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(SymTensor) function SymTenproduct_name_rename1(T1_,name1,T2_,name2,rename,whichname) result(T)
		type(SymTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:),name2(:),rename(2)
		integer,intent(in)::whichname
		integer :: i1(size(name1)),i2(size(name2))
		type(SymTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		character*50::oldName,newName
		if(.not.(Symif_original_dim(T1_%TenDim).and.Symif_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		oldName=rename(1)
		newName=rename(2)
		i1=SymNameorder(T1_%TenDim,name1)
		i2=SymNameorder(T2_%TenDim,name2)
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		if(whichname.eq.1)then
			call Symresetname(T1%TenDim,oldName,newName)
		end if
		call Symdimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		if(whichname.eq.2)then
			call Symresetname(T2%TenDim,oldName,newName)
		end if
		call Symdimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(SymTensor) function SymTenproduct_name_rename2(T1_,name1,rename1,T2_,name2,rename2) result(T)
		type(SymTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:),name2(:),rename1(2),rename2(2)
		integer :: i1(size(name1)),i2(size(name2))
		type(SymTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		character*50::oldName,newName
		if(.not.(Symif_original_dim(T1_%TenDim).and.Symif_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		i1=SymNameorder(T1_%TenDim,name1)
		i2=SymNameorder(T2_%TenDim,name2)
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		oldName=rename1(1)
		newName=rename1(2)
		call Symresetname(T1%TenDim,oldName,newName)
		call Symdimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		oldName=rename2(1)
		newName=rename2(2)
		call Symresetname(T2%TenDim,oldName,newName)
		call Symdimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(SymTensor) function SymTenproduct_int_name(T1_,i1,T2_,name2) result(T)
		type(SymTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name2(:)
		integer,intent(in)::i1(:)
		integer :: i2(size(name2))
		type(SymTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.(Symif_original_dim(T1_%TenDim).and.Symif_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		i2=SymNameorder(T2_%TenDim,name2)
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		call Symdimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		call Symdimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(SymTensor) function SymTenproduct_name_int(T1_,name1,T2_,i2) result(T)
		type(SymTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:)
		integer,intent(in)::i2(:)
		integer :: i1(size(name1))
		type(SymTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.(Symif_original_dim(T1_%TenDim).and.Symif_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		i1=SymNameorder(T1_%TenDim,name1)
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		call Symdimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		call Symdimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(SymTensor) function SymTenproduct_old(T1_,T2_,i1_,i2_) result(T)
		type(SymTensor),intent(in) :: T1_,T2_
		integer,intent(in) :: i1_(:),i2_(:)
		integer :: i,j,k,D2,decompoD1,decompoD2!Di use for decompese
		type(SymTensor) :: T1,T2
		integer :: decompoD1keep,decompoD2keep,rank1,&
			Tdim1(T1_%rank),Tdim2(T2_%rank),rank2,oper(2,3)
		integer::leni1,leni2,i1(size(i1_)),i2(size(i2_))
		if(.not.(Symif_original_dim(T1_%TenDim).and.Symif_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1_(1)
		i1(1)=i1_(1)
		do j=2,leni1
			if(i1(1).gt.i1_(j)) then
				i1(j)=i1_(j)
			else
				i1(j)=i1_(j)-1
			end if
		end do
		
		T2=T2_.pf.i2_(leni2)
		i2(leni2)=i2_(leni2)
		do j=leni2-1,1,-1
			if(i2(leni2).gt.i2_(j)) then
				i2(j)=i2_(j)+1
			else
				i2(j)=i2_(j)
			end if
		end do
		!T1
		do i=2,leni1
		
			T1=T1.pb.i1(i)
			
			do j=i+1,leni1
				if(i1(i).gt.i1(j)) then
					i1(j)=i1(j)
				else
					i1(j)=i1(j)-1
				end if
			end do
			
		end do
		call Symdimoperation(T1,(/1,rank1-leni1+1,rank1/))
		!T2
		do i=leni2-1,1,-1
			T2=T2.pf.i2(i)
			
			do j=i-1,1,-1
				if(i2(i).gt.i2(j)) then
					i2(j)=i2(j)+1
				else
					i2(j)=i2(j)
				end if
			end do
		
		end do
		call Symdimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
!**************** contract   ****************
!		combine two index of the Tensor,which is con_index and con_index+1
	type(SymTensor) function Symcontract(T1,con_index)
		integer,intent(in) :: con_index
		type(SymTensor),intent(in) :: T1
		type(SymDimension) ::newdim
		integer::i
		newdim=SymDimConstract(T1%TenDim,con_index,1)
		call SymstoreTenData(Symcontract,T1%block,T1%totalData)
		Symcontract%rank=SymDimSize(newDim)
		Symcontract%totalData=T1%totalData
		Symcontract%totalblock=T1%totalblock
		Symcontract%TenDim=newdim
		Symcontract%flag=.true.
		do i=1,Symcontract%totalData
			if(getflag(Symcontract%block(i)))then
				call dimoperation(Symcontract%block(i),(/1,con_index,con_index+1/))
			end if
		end do
		return
	end function
!**************** contracts   ****************
!		combine two index of the Tensor,which is index1 index1+1,..,index2-1,index2
!		if index2 larger than rank,the function will contract the index of index1 -> rank
	type(SymTensor) function Symcontracts(T1,index1,index2)
		integer,intent(in) :: index1,index2
		type(SymTensor),intent(in) :: T1
		type(SymDimension) ::newdim
		integer ::num,i
		num=index2-index1
		newdim=SymDimConstract(T1%TenDim,index1,num)
		call SymstoreTenData(Symcontracts,T1%block,T1%totalData)
		Symcontracts%rank=SymDimSize(newDim)
		Symcontracts%totalData=T1%totalData
		Symcontracts%totalblock=T1%totalblock
		Symcontracts%TenDim=newdim
		Symcontracts%flag=.true.
		do i=1,Symcontracts%totalData
			if(getflag(Symcontracts%block(i)))then
				call dimoperation(Symcontracts%block(i),(/1,index1,index2/))
			end if
		end do
		return
	end function			
	type(SymTensor) function Symcontractsv(T1,vector)
		integer,intent(in) ::vector(2)
		type(SymTensor),intent(in) :: T1
		type(SymDimension) ::newdim
		integer ::num,i
		num=vector(2)-vector(1)
		newdim=SymDimConstract(T1%TenDim,vector(1),num)
		call SymstoreTenData(Symcontractsv,T1%block,T1%totalData)
		Symcontractsv%rank=SymDimSize(newDim)
		Symcontractsv%totalData=T1%totalData
		Symcontractsv%totalblock=T1%totalblock
		Symcontractsv%TenDim=newdim
		Symcontractsv%flag=.true.
		do i=1,Symcontractsv%totalData
			if(getflag(Symcontractsv%block(i)))then
				call dimoperation(Symcontractsv%block(i),(/1,vector(1),vector(2)/))
			end if
		end do
		return
	end function		
!*****************  decompose  *****************
! decompose the de_index index of the Tensor into n(1),n(2)
!		for example the de_index index is (1,2,3,4,..inde,inde+1,...rank)
!		(1,2,3,4,..inde,inde+1,...rank)-->(1,2,3,4,..inde),(inde+1,...rank)		
!		if inde larger than rank ,the function will return no change	
	type(SymTensor) function Symdecompose(T1,de_index,inde)
		type(SymTensor),intent(in) :: T1
		integer,intent(in) :: de_index,inde
		type(SymDimension) :: newDim
		integer::i
		newDim=SymDimDecompose(T1%TenDim,de_index,inde)
		call SymstoreTenData(Symdecompose,T1%block,T1%totalData)
		Symdecompose%rank=SymDimSize(newDim)
		Symdecompose%totalData=T1%totalData
		Symdecompose%totalblock=T1%totalblock
		Symdecompose%TenDim=newdim
		Symdecompose%flag=.true.
		do i=1,Symdecompose%totalData
			if(getflag(Symdecompose%block(i)))then
				call dimoperation(Symdecompose%block(i),(/2,de_index,inde/))
			end if
		end do
		return
	end function
	type(SymTensor) function Symdecomposev(T1,vector)
		type(SymTensor),intent(in) :: T1
		integer,intent(in) ::vector(2)
		integer:: de_index,inde,i
		type(SymDimension) :: newDim
		de_index=vector(1)
		inde=vector(2)
		newDim=SymDimDecompose(T1%TenDim,de_index,inde)
		call SymstoreTenData(Symdecomposev,T1%block,T1%totalData)
		Symdecomposev%rank=SymDimSize(newDim)
		Symdecomposev%totalData=T1%totalData
		Symdecomposev%totalblock=T1%totalblock
		Symdecomposev%TenDim=newdim
		Symdecomposev%flag=.true.
		do i=1,Symdecomposev%totalData
			if(getflag(Symdecomposev%block(i)))then
				call dimoperation(Symdecomposev%block(i),(/2,de_index,inde/))
			end if
		end do
		return
	end function
			
	type(SymTensor) function Symdecompose1(T1,de_index)
		type(SymTensor),intent(in) :: T1
		integer,intent(in) :: de_index
		integer::i
		type(SymDimension) :: newDim
		newDim=SymDimDecompose(T1%TenDim,de_index,1)
		call SymstoreTenData(Symdecompose1,T1%block,T1%totalData)
		Symdecompose1%rank=SymDimSize(newDim)
		Symdecompose1%totalData=T1%totalData
		Symdecompose1%totalblock=T1%totalblock
		Symdecompose1%TenDim=newdim
		Symdecompose1%flag=.true.
		do i=1,Symdecompose1%totalData
			if(getflag(Symdecompose1%block(i)))then
				call dimoperation(Symdecompose1%block(i),(/2,de_index,1/))
			end if
		end do
		return
	end function
			
	type(SymTensor) function SymdecomposeAll(T1)
		type(SymTensor),intent(in) :: T1
		integer:: i
		type(SymDimension) :: newDim
		newDim=SymDimDecomposeAll(T1%TenDim)
		call SymstoreTenData(SymdecomposeAll,T1%block,T1%totalData)
		SymdecomposeAll%rank=SymDimSize(newDim)
		SymdecomposeAll%totalData=T1%totalData
		SymdecomposeAll%totalblock=T1%totalblock
		SymdecomposeAll%TenDim=newdim
		SymdecomposeAll%flag=.true.
		do i=1,SymdecomposeAll%totalData
			if(getflag(SymdecomposeAll%block(i)))then
				call dimoperation(SymdecomposeAll%block(i),(/3/))
			end if
		end do
		return
	end function
!**********************  compose_decompse  *********************
!		compose_decompse do the contract or decompose on the Tensor
!		Operar is a matrix,the first element of every row specify
!			1:contracts
!			2:decompose
!			3:decomposeAll
!		other element of every row is input parameter of the function
	type(SymTensor) function Symcompose_decompse(T1,Operar)result(compose_decompse)
		integer,intent(in)::Operar(:,:)
		type(SymTensor),intent(in) :: T1
		integer::i,j,size1
		integer :: index1,index2
		integer ::num
		integer :: de_index,inde
		type(Symdimension)::dimen
		if(.not.T1%flag)then
			write(*,*)"There is no data in the Tensor,(Symcompose_decompse)"
			stop
		end if
		size1=size(Operar,1)
		dimen=T1%TenDim
		do i=1,size1
			if(Operar(i,1).eq.1) then
				index1=Operar(i,2)
				index2=Operar(i,3)
				num=index2-index1
				dimen=SymDimConstract(dimen,index1,num)
			else if(Operar(i,1).eq.2) then
				de_index=Operar(i,2)
				inde=Operar(i,3)
				dimen=SymDimDecompose(dimen,de_index,inde)
			else if(Operar(i,1).eq.3) then
				dimen=SymDimDecomposeAll(dimen)
			else
				write(*,*) "error in compose_decompse"
				stop
			end if
		end do
		
		!call SymstoreTenData(compose_decompse,T1%block,T1%totalData)
		allocate(compose_decompse%block(T1%totalData))
		do i=1,T1%totalData
			if(getFlag(T1%block(i)))then
					compose_decompse%block(i)=T1%block(i).cd.Operar
			end if
		end do
		
		compose_decompse%TenDim=dimen
		compose_decompse%rank=SymDimSize(dimen)
		compose_decompse%flag=.true.
		compose_decompse%totalData=T1%totalData
		compose_decompse%totalblock=T1%totalblock
		return
	end function
!********************  dimOperations  ****************************
!		 dimOperations do the contract or decompose on the Tensor with
!	  no change on the Tensot_data in the Tensor.
!		Operar is a matrix,the first element of every row specify
!		 1:contracts
!		 2:decompose
!		 3:decomposeAll
!		other element of every row is input parameter of the function
	subroutine symdimOperations(T1,Operar)
		integer,intent(in)::Operar(:,:)
		type(symTensor),intent(inout) :: T1
		integer::i,j,size1
		integer :: index1,index2
		integer ::num
		integer :: de_index,inde
		if(.not.T1%flag)then
			write(*,*)"There is no data in the Tensor,(compose_decompse)"
			stop
		end if
		size1=size(Operar,1)
		do i=1,size1
			if(Operar(i,1).eq.1) then
				index1=Operar(i,2)
				index2=Operar(i,3)
				num=index2-index1
				T1%TenDim=symDimConstract(T1%TenDim,index1,num)
				T1%rank=symDimSize(T1%TenDim)
			else if(Operar(i,1).eq.2) then
				de_index=Operar(i,2)
				inde=Operar(i,3)
				T1%TenDim=symDimDecompose(T1%TenDim,de_index,inde)
				T1%rank=symDimSize(T1%TenDim)
			else if(Operar(i,1).eq.3) then
				T1%TenDim=symDimDecomposeAll(T1%TenDim)
				T1%rank=symDimSize(T1%TenDim)
			else
				write(*,*) "error in compose_decompse"
				stop
			end if
		end do
		call dimOperation(T1%block,Operar)
		return
	end subroutine
	subroutine symdimOperation1(T1,Operat)
		integer,intent(in)::Operat(:)
		type(symTensor),intent(inout) :: T1
		integer,allocatable::Operat2(:,:)
		allocate(Operat2(1,size(Operat)))
		Operat2(1,:)=Operat
		call symdimOperations(T1,Operat2)
		return
	end subroutine
!***************** equal_of_Tensor *****************
	logical function Symequal_of_Tensor(T1,T2)result(equal_of_Tensor)
		type(SymTensor),intent(in) :: T1,T2
		integer :: i
		equal_of_Tensor=.true.
		if(T1%rank.ne.T2%rank) then
			equal_of_Tensor=.false.
			return
		end if
		if(.not.Symequal_of_dim(T1%TenDim,T2%TenDim)) then
			equal_of_Tensor=.false.
			return
		end if
		do i=1,T1%TotalData
			if(getFlag(T1%block(i)).and.(.not.getFlag(T2%block(i))))then
				equal_of_Tensor=.false.
				return
			end if
			if(getFlag(T2%block(i)).and.(.not.getFlag(T1%block(i))))then
				equal_of_Tensor=.false.
				return
			end if
			if(getFlag(T2%block(i)).and.(getFlag(T1%block(i))))then
				if(.not.(T2%block(i).equ.T1%block(i)))then
					equal_of_Tensor=.false.
					return
				end if
			end if
		end do
		return
	end function
!
!  s1
! -->--
!      \     s3
!       |--->--- here is the case of s1+s1-s3=0
! -->--/
!   s2
!
!  use for fusing two index of a Tensor into one index
!    the quantum number of the two index are s1 and s2
!    fuse into the quantum number of s3=|s1|+|s2|
!  The Tensor to be fuse product the two index with this
!  fuse Tensor will come out what we want.
!       __   s1                           __  
!      |  |-->--  -->--\                 |  |
! -->--|  |      *     |-->--s3  =  -->--|  |-->-- s3
!      |__|-->--  -->--/                 |__|
!             s2
	type(SymTensor) function fuseTensorDim_noName(dimensi,ith,jth,outrule,newQN)result(fuse)
		type(SYmdimension),intent(in)::dimensi
		integer,intent(in)::ith,jth,outrule
		real*4,optional,intent(in)::newQN
		real*4::QN_i(3),QN(3),minQN(3)
		integer::rule(3),maxindex(3),i,address,inde(3),deg(3),LDR,length
		integer,allocatable::degeneracy(:)
		real*4,allocatable::QNorder(:,:)
		integer,allocatable::identity_index(:,:)
		integer::temp1,temp2
		call SymRule(rule(1),dimensi,ith)
		call SymRule(rule(2),dimensi,jth)
		rule(1)=-1*rule(1)! the rule of the fuse is the opposite direction of the dimen
		rule(2)=-1*rule(2)! the rule of the fuse is the opposite direction of the dimen
		rule(3)=outrule
		QN(1)=SymQNum(dimensi,ith)
		QN(2)=SymQNum(dimensi,jth)
		if(present(newQN))then
			QN(3)=newQN
		else
			QN(3)=QN(1)+QN(2)
		end if
		length=2*QN(3)+1
		allocate(degeneracy(length))
		length=100*length!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! may be wrong
		allocate(QNorder(length,3))
		allocate(identity_index(length,3))
!Dimension		
		call U1QNtoMaxIndex(maxindex,QN)
		fuse%totalData=product(maxindex)
		fuse%TenDim=maxindex
		fuse%rank=3
		call SymDimDegeneracy(degeneracy,dimensi,ith)
		call set_QN1(fuse%TenDim,1,QN(1),degeneracy(1:(fuse%TenDim.i.1)),rule(1))
		call SymDimDegeneracy(degeneracy,dimensi,jth)
		call set_QN1(fuse%TenDim,2,QN(2),degeneracy(1:(fuse%TenDim.i.2)),rule(2))
!the last Dimension
		call fuse_degeneracy(degeneracy,QNorder,identity_index,length,fuse%TenDim,QN,rule)
		!(outdeg,outQN,QNorder,identity_index,dimen,ith,jth,QN,rule)
		
	!	write(*,*)QN(3),degeneracy((fuse%TenDim.i.3)),size(degeneracy)!The length of degeneracy is too short!!!!!!
		call set_QN1(fuse%TenDim,3,QN(3),degeneracy(1:(fuse%TenDim.i.3)),rule(3))
		
		allocate(fuse%block(fuse%totalData))
	!	call SymDprint0(fuse%TenDim)
	!	read(*,*)
		fuse%totalblock=length
		do i=1,length
			call Symindex(inde,fuse%TenDim,QNorder(i,:))
			address=addressToIndes2(maxindex,inde)
			call SymDimDegeneracy(deg,fuse%TenDim,inde)
			fuse%block(address)=zeroTen(deg)
			LDR=deg(1)*deg(2)
			temp1=identity_index(i,1)+identity_index(i,3)-1
			temp2=identity_index(i,2)+identity_index(i,3)-1
!			write(*,*)gettotalData(fuse%block(address)),temp1*temp2,identity_index(i,3)
			call fuse_combination_identity_data(fuse%block(address)%Tensor_Data&
					,identity_index(i,1),identity_index(i,2),identity_index(i,3),LDR)
		end do
		fuse%flag=.true.
		return
	end function
	type(SymTensor) function fuseTensor_noName(ST,ith,jth,outrule,newQN)result(fuse)
		type(SYmTensor),intent(in)::ST
		integer,intent(in)::ith,jth,outrule
		real*4,optional,intent(in)::newQN
		real*4::QN_i(3),QN(3),minQN(3)
		integer::rule(3),maxindex(3),i,address,inde(3),deg(3),LDR,length
		real*4,allocatable::QNorder(:,:)
		integer,allocatable::identity_index(:,:),degeneracy(:)
		call SymRule(rule(1),ST%TenDim,ith)
		call SymRule(rule(2),ST%TenDim,jth)
		rule(1)=-1*rule(1)! the rule of the fuse is the opposite direction of the dimen
		rule(2)=-1*rule(2)! the rule of the fuse is the opposite direction of the dimen
		rule(3)=outrule
		QN(1)=SymQNum(ST%TenDim,ith)
		QN(2)=SymQNum(ST%TenDim,jth)
		if(present(newQN))then
			QN(3)=newQN
		else
			QN(3)=QN(1)+QN(2)
		end if
		length=2*QN(3)+1
		allocate(degeneracy(length))
		length=100*length
		allocate(QNorder(length,3))
		allocate(identity_index(length,3))
!Dimension		
		call U1QNtoMaxIndex(maxindex,QN)
		fuse%totalData=product(maxindex)
		fuse%TenDim=maxindex
		fuse%rank=3
		call SymDimDegeneracy(degeneracy,ST%TenDim,ith)
		call set_QN1(fuse%TenDim,1,QN(1),degeneracy(1:(fuse%TenDim.i.1)),rule(1))
		call SymDimDegeneracy(degeneracy,ST%TenDim,jth)
		call set_QN1(fuse%TenDim,2,QN(2),degeneracy(1:(fuse%TenDim.i.2)),rule(2))
!the last Dimension
		call fuse_degeneracy(degeneracy,QNorder,identity_index,length,fuse%TenDim,QN,rule)
		!(outdeg,outQN,QNorder,identity_index,dimen,ith,jth,QN,rule)
		call set_QN1(fuse%TenDim,3,QN(3),degeneracy(1:(fuse%TenDim.i.3)),rule(3))
		
		allocate(fuse%block(fuse%totalData))
		fuse%totalblock=length
		do i=1,length
			call Symindex(inde,fuse%TenDim,QNorder(i,:))
			address=addressToIndes2(maxindex,inde)
			call SymDimDegeneracy(deg,fuse%TenDim,inde)
			fuse%block(address)=zeroTen(deg)
			LDR=deg(1)*deg(2)
			call fuse_combination_identity_data(fuse%block(address)%Tensor_Data&
					,identity_index(i,1),identity_index(i,2),identity_index(i,3),LDR)
		end do
		fuse%flag=.true.
		return
	end function
	type(SymTensor) function fuseTensor_Name(ST,name1,name2,newname,outrule,newQN)result(fuse)
		type(SYmTensor),intent(in)::ST
		integer,intent(in)::outrule
		CHARACTER(len=*),intent(in)::name1,name2,newname
		real*4,optional,intent(in)::newQN
		integer::ith,jth
		ith=SymTenNameOrder(ST,name1)
		jth=SymTenNameOrder(ST,name2)
		fuse=fuseTensor_noName(ST,ith,jth,outrule,newQN)
		call SetSymTensorName(fuse,1,name1)
		call SetSymTensorName(fuse,2,name2)
		call SetSymTensorName(fuse,3,newname)
		return
	end function
	type(SymTensor) function fuseTensorDim_Name(dimen,name1,name2,newname,outrule,newQN)result(fuse)
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::outrule
		CHARACTER(len=*),intent(in)::name1,name2,newname
		real*4,optional,intent(in)::newQN
		integer::ith,jth
		ith=SymNameorder(dimen,name1)
		jth=SymNameorder(dimen,name2)
		fuse=fuseTensorDim_noName(dimen,ith,jth,outrule,newQN)
		call SetSymTensorName(fuse,1,name1)
		call SetSymTensorName(fuse,2,name2)
		call SetSymTensorName(fuse,3,newname)
		return
	end function
	type(SymTensor) function fuseTensor_noName_array(ST,ith,outrule)result(fuse)
		type(SymTensor),intent(in)::ST
		integer,intent(in)::outrule,ith(:)
		type(Symdimension)::Dimen
		type(SymTensor)::tempfuse
		integer::i,lenindex
		logical::first_flag
		lenindex=size(ith)
		dimen=(ST.subDim.ith(1))
		first_flag=.true.
		do i=2,lenindex
		!	write(*,*)"here1",i
			dimen=dimen+(ST.subDim.ith(i))
			if(first_flag)then
				fuse=fuseTensordim_noName(dimen,1,2,outrule)
				first_flag=.false.
			else
		!		call SymDprint0(dimen)
				
				tempfuse=fuseTensordim_noName(dimen,1,2,outrule)
		!		stop
				fuse=SymTenproduct(fuse,(/SymgetRank(fuse)/),tempfuse,(/1/))
			end if
		!	write(*,*)"here2",i
			dimen=fuse.subdim.SymgetRank(fuse)
		!	write(*,*)"here3",i
		end do
		return
	end function
	type(SymTensor) function fuseTensor_Name_array(ST,dimname,newname,outrule)result(fuse)
		type(SymTensor),intent(in)::ST
		integer,intent(in)::outrule
		CHARACTER(len=*),intent(in)::dimname(:),newname
		type(Symdimension)::Dimen
		type(SymTensor)::tempfuse
		integer::i,lenindex
		logical::first_flag
		lenindex=size(dimname)
		dimen=(ST.subDim.dimname(1))
		first_flag=.true.
		do i=2,lenindex
			dimen=dimen+(ST.subDim.dimname(i))
			if(first_flag)then
				fuse=fuseTensorDim_Name(dimen,dimname(i-1),dimname(i),newname,outrule)
				first_flag=.false.
			else
				tempfuse=fuseTensorDim_Name(dimen,newname,dimname(i),newname,outrule)
				fuse=SymTenproduct(fuse,(/SymgetRank(fuse)/),tempfuse,(/1/))
			end if
			dimen=fuse.subdim.SymgetRank(fuse)
		end do
		return
	end function
			
			
			
!calculate the degeneracy for the newdiemn(fuse)
!example:QN are 0.5 and 0.5 ,new QN will be 1. 0.5 + 0.5 -> 1
! and the degeneracy are 0.5:[2,3], and 0.5:[4,5]
!-0.5 -0.5 = -1 ,
!   degeneracy for -1. is degeneracy(-0.5)*degeneracy(-0.5)+degeneracy(-1)
!   and degeneracy(-0.5)=2,degeneracy(-0.5)=4,degeneracy(-1)=0
!   so degeneracy(-1)=8
!+0.5 -0.5 = 0. 
!   degeneracy for 0. is degeneracy(+0.5)*degeneracy(-0.5)+degeneracy(0.)
!   and degeneracy(+0.5)=3,degeneracy(-0.5)=4,degeneracy(-1)=0
!   so degeneracy(0.)=12
!-0.5 +0.5 = 0. 
!   degeneracy for 0. is degeneracy(-0.5)*degeneracy(+0.5)+degeneracy(0.)
!   and degeneracy(-0.5)=2,degeneracy(+0.5)=5,degeneracy(0.)=12
!   so degeneracy(0.)=2*4+12=20
!+0.5 +0.5 = +1 ,
!   degeneracy for +1. is degeneracy(+0.5)*degeneracy(+0.5)+degeneracy(+1)
!   and degeneracy(+0.5)=3,degeneracy(+0.5)=15,degeneracy(+1)=0
!   so degeneracy(-1)=15
!At last the degeneracy are
!QN:  [-1., 0. , 1.]
!deg: [8  , 20 , 15]
!QNorder:
!/ -0.5,-0.5,-1 \
!| +0.5,-0.5, 0 |
!| -0.5,+0.5, 0 |
!\ +0.5,+0.5,+1 /
!identity_index,The two index of QN=0.5and QN=0.5 will be regard as one index
!/ 1 ,1 ,8  \ 
!| 1 ,1 ,12 |
!| 1 ,13,8  |
!\ 1 ,1 ,15 /
!This data used for building the identity matrix in fuse
!for example if one line in identity_index are
! [1,1,3] mean the identity is the fist 1d. is in (1,1) of the matrix are there are 3 1.
! /1 0 0 \
! |0 1 0 |
! \0 0 1 /
! if one line in identity_index are
! [2,1,2] mean the identity is the fist 1d. is in (2,1) of the matrix are there are 2 1.
! /0 0 0 \
! |0 1 0 |
! \0 0 1 /
	subroutine fuse_degeneracy(outdeg,QNorder,identity_index,length,dimen,QN,rule)
		integer,intent(inout)::outdeg(:),length!the length of outdeg in output
		real*4,allocatable,intent(inout)::QNorder(:,:)
		integer,allocatable,intent(inout)::identity_index(:,:)
		type(SYmDimension),intent(in)::dimen
		real*4,allocatable::outQN(:)
		integer,intent(in)::rule(3)
		real*4,intent(in)::QN(3)
		real*4::QN_i(3),minQN(3)
		logical::goon
		integer::i,deg_i,QNdeg(2),QNcouter
		allocate(outQN(size(outdeg)))
		minQN=-1*QN
		goon=.true.
		QN_i=minQN
		outQN=QN(3)+1!initial outQN that will be met in the code,if cannot find,it will rewrite
		QNcouter=1!the index of outQN that is store data
		outdeg=0
		length=0
		i=1
		do while(goon)
			if(U1Rule(QN_i,rule))then!if[1]
				deg_i=fuse_degeneracy_QNindex(outQN,QN_i(3),QNcouter)
				if(deg_i.gt.size(outdeg))write(*,*)"ERROR1"
				call SymDimDegeneracy(QNdeg,dimen,(/1,2/),QN_i(1:2))
				if(i.gt.size(QNorder,1))write(*,*)"ERROR2"
				QNorder(i,:)=QN_i
				identity_index(i,:)=(/1,outdeg(deg_i)+1,product(QNdeg)/)
				length=length+1
				outdeg(deg_i)=outdeg(deg_i)+identity_index(i,3)
				i=i+1
			end if!if[1]
			goon=inde_counter(QN_i,minQN,QN,1.)	
			if(i.gt. 1d6)stop
		end do
		return
	end subroutine
!outdeg=[2,2,3]
!outQN=[-1.,0.,1.]
!input -1, output 1,the first element in outQN	
!use in fuse_degeneracy
	integer function fuse_degeneracy_QNindex(outQN,QN_i,ith)
		real*4,intent(inout)::outQN(:),QN_i
		integer,intent(inout)::ith
		integer::i
		do i=1,ith-1
			if(abs(QN_i-outQN(i)).lt.zero_error_)then
				fuse_degeneracy_QNindex=i
				return
			end if
		end do
		fuse_degeneracy_QNindex=ith
		outQN(ith)=QN_i
		ith=ith+1
		return
		stop
	end function
!fusedata :
! 0 0 0 0
! 0 0 0 0
!if row,col,length_i,LDR are 1,1,2,2
!then output
! 1 0 0 0
! 0 1 0 0
!if row,col,length_i,LDR are 1,3,2,2
! 0 0 1 0
! 0 0 0 1
!input fusedata should be all 0., use to modify Tensor%Tensor_Data
	subroutine fuse_combination_identity_data(fusedata,row,col,length_i,LDR)
		integer,intent(in)::row,col,length_i,LDR
		complex*16,intent(inout)::fusedata(LDR,*)
		integer::i
		do i=0,length_i-1
			fusedata(row+i,col+i)=1d0
		end do
		return
	end subroutine
	
	subroutine symSVDcutoff1(T,U,s,V,minNewQN,degeneracy) !in U(1) symmetry , T must be a partitioned matrix
		type(SymTensor),intent(in)::T
		type(SymTensor),intent(out)::U,s,V
		real*4,intent(in)::minNewQN
		integer,intent(in)::degeneracy(:)
		integer::minTi(2),minUi(2),minsi(2),minVi(2),Tm,Tn,Um,Un,sm,sn,Vm,Vn,totalblock,newRule,newblock
		type(SymDimension) :: T1Dim,T2Dim,sdim,newDim
		real*4::QN1,QN2,minQN
		if(T%rank.ne.2) then
			if(.not.T%flag) then
				write(*,*)"ERROR in svd"
				write(*,*)"There is no data in the Tensor"
				write(*,*)"stop"
				stop
			end if
			write(*,*) "Input Tensor should be 2 dimension in svd"
			write(*,*)"stop"
			stop
		end if
		QN1=SymTensorQN(T,1)
		QN2=SymTensorQN(T,2)
		 
			Un=2*minNewQN+1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!U(1) symmetry determine the diemnison
			newDim=(/Un/)
			call SymRule(newRule,T%TenDim,2)
			call setDimQN(newDim,1,minNewQN,degeneracy,newRule)
			T1Dim=T.subdim.1
			T1Dim=T1Dim+newDim
			
			sdim=newDim
			call reverseDimrule(newDim,1)
			sdim=newDim+sdim
			T2Dim=newDim
			T2Dim=T2Dim+(T.subdim.2)
			
		
		call allocateSymTensor(U,T1Dim)
		call allocateSymTensor(s,sDim)
		call allocateSymTensor(V,T2Dim)
		
		minQN=-1.*minNewQN
		minTi(1)=SymTensorindex(T,1,minQN)
		minTi(2)=SymTensorindex(T,2,minQN)
		minUi(1)=SymTensorindex(U,1,minQN)
		minUi(2)=SymTensorindex(U,2,minQN)
		minsi=1
		minVi(1)=SymTensorindex(V,1,minQN)
		minVi(2)=SymTensorindex(V,2,minQN)
		Tm=T.dim.1
		Tn=T.dim.2
		Um=T1Dim.i.1
		Un=T1Dim.i.2
		sm=sdim.i.1
		sn=sdim.i.2
		Vm=T2Dim.i.1
		Vn=T2Dim.i.2
		totalblock=sm
		call symSVDcutoff_svd(totalblock,newblock,minTi,minUi,minSi,minVi,T%block,&
				Tm,Tn,U%block,Um,Un,s%block,sm,sn,V%block,Vm,Vn,degeneracy) 
		U%TotalBlock=newblock
		s%TotalBlock=newblock
		V%TotalBlock=newblock
		return
	end subroutine	
!minNewQN and Deg
!0.5:       1  1
!1  :      1  2  1
!1.5:     1  3  3  1
!2  :    1  4  6  4  1		
!2.5:   1 5 10 10  5  1 
!3  :  1 6 15 20 15 6  1
!3.5: 1 7 21 35 35 21 7 1
!These QNs or N*Deg give a good result on MPS	
	subroutine symSVDcutoff2(T,U,s,V,maxNewQN_,Num_) !in U(1) symmetry , T must be a partitioned matrix
		type(SymTensor),intent(in)::T
		type(SymTensor),intent(out)::U,s,V
		real*4,intent(in)::maxNewQN_
		integer,optional,intent(in)::Num_
		integer::minTi(2),minUi(2),minsi(2),minVi(2),Tm,Tn,Um,Un,sm,sn,Vm,Vn,totalblock,newRule,newblock
		integer::tempint,i,minindex,Num
		logical::flag
		type(SymDimension) :: T1Dim,T2Dim,sdim,newDim
		integer,allocatable::degeneracy(:)
		real*4::QN1,QN2,minQN,maxNewQN
		if(T%rank.ne.2) then
			if(.not.T%flag) then
				write(*,*)"ERROR in svd"
				write(*,*)"There is no data in the Tensor"
				write(*,*)"stop"
				stop
			end if
			write(*,*) "Input Tensor should be 2 dimension in svd"
			write(*,*)"stop"
			stop
		end if
		QN1=SymTensorQN(T,1)
		QN2=SymTensorQN(T,2)
		if(present(Num_))then
			Num=Num_
		else
			Num=1
		end if
		if(QN1.lt.QN2)then!find the minQN, cut off if maxNewQN less than it
			minQN=QN1
			minindex=1
		else
			minQN=QN2
			minindex=2
		end if
		call SymRule(newRule,T%TenDim,2)
		if(maxNewQN_.lt.minQN)then!if T is [-2.5,0.5],maxNewQN=1.,cut to [-2.5,0.5]
											!if T is [-3.,2.],maxNewQN=1.,cut to [-3.,1.]
											!if T is [-3.,2.],maxNewQN=1.5,cut to [-3,1]
			tempint=2*minQN+1      !if T is [-2.5,0.5],maxNewQN=1.5,cut to [-2.5,1.5]
			flag=.false.
			do i=tempint,1,-1
				if(SymTensorQN(T,minindex,i).eq.maxNewQN_)then
					flag=.true.
					exit
				end if
			end do
			if(flag)then
				maxNewQN=maxNewQN_
			else
				maxNewQN=maxNewQN_-0.5
			end if
			newDim=SymU1Dim2((/maxNewQN/),(/newRule/),(/Num/))
		else
			maxNewQN=minQN
			newDim=T.subDim.minindex
			if(minindex.eq.1)then
				call reverseDimrule(newDim)
			end if
		end if
		tempint=2*maxNewQN+1
		allocate(degeneracy(tempint))
		call SymDimDegeneracy(degeneracy,newDim,1)
		T1Dim=T.subdim.1
		T1Dim=T1Dim+newDim
		
		sdim=newDim
		call reverseDimrule(newDim,1)
		sdim=newDim+sdim
		T2Dim=newDim
		T2Dim=T2Dim+(T.subdim.2)
		
			
			
		
		call allocateSymTensor(U,T1Dim)
		call allocateSymTensor(s,sDim)
		call allocateSymTensor(V,T2Dim)
		
		minQN=-1.*maxNewQN
		minTi(1)=SymTensorindex(T,1,minQN)
		minTi(2)=SymTensorindex(T,2,minQN)
		minUi(1)=SymTensorindex(U,1,minQN)
		minUi(2)=SymTensorindex(U,2,minQN)
		minsi=1
		minVi(1)=SymTensorindex(V,1,minQN)
		minVi(2)=SymTensorindex(V,2,minQN)
		Tm=T.dim.1
		Tn=T.dim.2
		Um=T1Dim.i.1
		Un=T1Dim.i.2
		sm=sdim.i.1
		sn=sdim.i.2
		Vm=T2Dim.i.1
		Vn=T2Dim.i.2
		totalblock=sm
		call symSVDcutoff_svd(totalblock,newblock,minTi,minUi,minSi,minVi,T%block,&
				Tm,Tn,U%block,Um,Un,s%block,sm,sn,V%block,Vm,Vn,degeneracy) 
		U%TotalBlock=newblock
		s%TotalBlock=newblock
		V%TotalBlock=newblock
		return
	end subroutine		
	subroutine symSVDcutoff_svd(lenblock_,newblock,minTi,minUi,minSi,minVi,T,LDRT,LDCT,U,LDRU,LDCU,s,LDRS,LDCS,V,LDRV,LDCV,degeneracy) 
		integer,intent(in)::minTi(2),minUi(2),minsi(2),minVi(2),lenblock_
		integer,intent(in)::LDRT,LDCT,LDRU,LDCU,LDRS,LDCS,LDRV,LDCV
		type(Tensor),intent(in)::T(LDRT,LDCT)
		type(Tensor),intent(out)::U(LDRU,LDCU),s(LDRS,LDCS),V(LDRV,LDCV)
		integer,intent(out)::newblock
		integer,intent(in)::degeneracy(:)
		integer::i,j,rowT,colT,rowU,colU,rowS,colS,rowV,colV,lenblock
		rowT=minTi(1)-1
		colT=minTi(2)-1
		rowU=minUi(1)-1
		colU=minUi(2)-1
		rowS=minSi(1)-1
		colS=minSi(2)-1
		rowV=minVi(1)-1
		colV=minVi(2)-1
		lenblock=lenblock_
		i=1
		j=1
		newblock=0
		do while(i.le.lenblock)
			if(getflag(T(rowT+i,colT+i))) then
				call SVDcutoff(T(rowT+i,colT+i),U(rowU+i,colU+i),s(rowS+i,colS+i),V(rowV+i,colV+i),degeneracy(j))
				s(rowS+i,colS+i)=eye(s(rowS+i,colS+i))
				j=j+1
				newblock=newblock+1
			else
				lenblock=lenblock+1
			end if
			i=i+1
		end do
		return
	end subroutine
! modify all the dimension in the Tensor
! Will not test symmetry
	subroutine Symresetdim1(Ten,QN,rule)
		type(SymTensor),intent(inout)::Ten
		integer,intent(in)::rule(:)
		real*4,intent(in)::QN(:)
		integer::test
		Ten%TenDim=SymU1Dim(QN,rule)
		Ten%rank=size(QN)
		test=product(2*QN+1)
		if(test.ne.Ten%Totaldata) then
			write(*,*)"ERROR in resetdim1"
			write(*,*)test,Ten%Totaldata
			call SymTDprint(Ten)
			write(*,*)"stop"
			stop
		end if
		return
	end subroutine	
	subroutine Symresetdim2(Ten,dimen)
		type(SymTensor),intent(inout)::Ten
		type(Symdimension),intent(in)::dimen
		if(Symouttotaldata(dimen).ne.Ten%Totaldata) then
			write(*,*)"ERROR in resetdim2"
			write(*,*)Symouttotaldata(dimen),Ten%Totaldata
			call SymTDprint(Ten)
			write(*,*)"stop"
			stop
		end if
		Ten%TenDim=dimen
		Ten%rank=SymDimSize(Dimen)
		return
	end subroutine	
!***************   dot   ******************
!		return <phi1|phi2>
!	dot product conjugating the first vector,The Tensor will be regard as a vector
	complex*16 function Symdot(phi1,phi2)
		Type(SymTensor),intent(in)::phi1,phi2
		integer::N1,N2
		N1=SymgetTotalData(phi1)
		N2=SymgetTotalData(phi2)
		if(N1.ne.N2) then
			write(*,*)"ERROR in Symdot,do not finish this part"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			stop
		end if
		if(.not.phi1%flag)then
			write(*,*)"There is no data in the Tensor,(Symdot)"
			stop
		end if
		Symdot=SZDOTC(N1,phi1%block,1,phi2%block,1)
		RETURN
	end function
!the same as dot product but do not conjugating the first vector,The Tensor will be regard as a vector
	complex*16 function numdot(phi1,phi2)
		Type(SymTensor),intent(in)::phi1,phi2
		integer::N1,N2
		N1=SymgetTotalData(phi1)
		N2=SymgetTotalData(phi2)
		if(N1.ne.N2) then
			write(*,*)"ERROR in Symdot,do not finish this part"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			stop
		end if
		if(.not.phi1%flag)then
			write(*,*)"There is no data in the Tensor,(Symdot)"
			stop
		end if
		numdot=SZDOTU(N1,phi1%block,1,phi2%block,1)
		RETURN
	end function
	

!*****************  norm   *****************
!			return  <phi|phi>
	real*8 function symnorm2(Tr)
		type(SymTensor),intent(in) :: Tr
		if(.not.SymgetFlag(Tr))then
			write(*,*)"There is no data in the Tensor,(norm2)"
			stop
		end if
		symnorm2=SDZNRM2(Tr%TotalData,Tr%block,1)
		return
	end function	
!		return  sqrt(<phi|phi>	)
	real*8 function symnorm(Tr)
		type(SymTensor),intent(in) :: Tr
		if(.not.SymgetFlag(Tr))then
			write(*,*)"There is no data in the Tensor,(norm)"
			stop
		end if
		symnorm=SDZNRM2(Tr%TotalData,Tr%block,1,.true.)
		return
	end function				
	real*8 FUNCTION SDZNRM2(N,X,INCX,sqrtflag)!znrm2 for block data,Not finished but it could run
		integer,intent(in)::incx,n
		type(Tensor),intent(in)::X(*)
		logical,optional,intent(in)::sqrtflag
		integer::i
		SDZNRM2=0d0
		do i=1,N
			if(getFlag(X(i)))then
				SDZNRM2=SDZNRM2+norm2(X(i))
			end if
		end do
		if(present(sqrtflag))then
			SDZNRM2=dsqrt(SDZNRM2)
		end if
		return
	end function
!*****************  Element   *****************
	type(Tensor) function SymElement(T,Tdim)
		type(SymTensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		integer::Inde
		if(.not.T%flag)then
			write(*,*)"There is no data in the SymTensor,(.i.)"
			stop
		end if
		inde=addressToIndes(T,Tdim)
		SymElement=T%block(inde)
		return
	end function		  
	type(Tensor) function SymElement2(T,inde)
		type(SymTensor),intent(in) ::T
		integer,intent(in)::inde
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(.i.)"
			stop
		end if
		SymElement2=T%block(inde)
		return
	end function
!*****************  Htranspose  *****************
	type(SymTensor) function SymHtranspose(T)result(Htranspose)
		type(SymTensor),intent(in) :: T
		integer::rank,m,n,i
		rank=T%rank
		if(rank.eq.1) then
			m=T.dim. 1
			allocate(Htranspose%block(m))
			do i=1,m
				if(getflag(T%block(i)))then
					Htranspose%block(i)=.H.T%block(i)
				end if
			end do
			Htranspose%rank=T%rank
			Htranspose%totalData=T%totalData
			Htranspose%totalblock=T%totalblock
			Htranspose%TenDim=T%TenDim
			Htranspose%flag=.true.
		else if(rank.eq.2) then
			Htranspose=.p.T
			do i=1,Htranspose%TotalData
				if(getflag(Htranspose%block(i)))then
					Htranspose%block(i)=.con.Htranspose%block(i)
				end if
			end do
		else
			write(*,*) "SymTensor should be 1 or 2 dimension"
			write(*,*) "ERROR in SymHtranspose,stop"
			stop
		end if
		return
	end function
!Htranspose and reverseRule
	type(SymTensor) function SymHtransposeReverse(T)result(Htranspose)
		type(SymTensor),intent(in) :: T
		integer::rank,m,n,i
		rank=T%rank
		if(rank.eq.1) then
			m=T.dim. 1
			allocate(Htranspose%block(m))
			do i=1,m
				if(getflag(T%block(i)))then
					Htranspose%block(i)=.H.T%block(i)
				end if
			end do
			Htranspose%rank=T%rank
			Htranspose%totalData=T%totalData
			Htranspose%totalblock=T%totalblock
			Htranspose%TenDim=T%TenDim
			call reverseDimrule(Htranspose%TenDim)
			Htranspose%flag=.true.
		else if(rank.eq.2) then
			Htranspose=.p.T
			call reverseDimrule(Htranspose%TenDim)
			do i=1,Htranspose%TotalData
				if(getflag(Htranspose%block(i)))then
					Htranspose%block(i)=.con.Htranspose%block(i)
				end if
			end do
		else
			write(*,*) "SymTensor should be 1 or 2 dimension"
			write(*,*) "ERROR in SymHtranspose,stop"
			stop
		end if
		return
	end function
!*****************  conjugate  *****************
	type(SymTensor) function Symconjugate(T)result(conjugate)
		type(SymTensor),intent(in) :: T
		integer :: m,i
		m=T%totalData
		allocate(conjugate%block(m))
		do i=1,m
			if(getFlag(T%block(i)))then
				conjugate%block(i)=.con.T%block(i)
			end if
		end do
		conjugate%rank=T%rank
		conjugate%totalData=T%totalData
		conjugate%totalblock=T%totalblock
		conjugate%TenDim=T%TenDim
		conjugate%flag=.true.
		return
	end function
!conjugate and reverseRule
	type(SymTensor) function Symconjugatereverse(T)result(conjugate)
		type(SymTensor),intent(in) :: T
		integer :: m,i
		m=T%totalData
		allocate(conjugate%block(m))
		do i=1,m
			if(getFlag(T%block(i)))then
				conjugate%block(i)=.con.T%block(i)
			end if
		end do
		conjugate%rank=T%rank
		conjugate%totalData=T%totalData
		conjugate%totalblock=T%totalblock
		conjugate%TenDim=T%TenDim
		call reverseDimrule(conjugate%TenDim)
		conjugate%flag=.true.
		return
	end function
!*****************  one_com  *****************
	type(SymTensor) function symeye_one(QN,degeneracy)
		real*4,intent(in) :: QN(:)
		integer,intent(in)::degeneracy(:,:)
		type(SymDimension)::dimen
		integer::maxdim(2)
		if(size(QN).gt.2)then
			write(*,*)"ERROR in symeye_one"
			stop
		end if
		maxdim=2*QN+1
		dimen=SymU1Dim(QN,degeneracy,(/1,-1/))
		call allocateSymTensor(symeye_one,dimen)
		symeye_one%totalblock=minval(maxdim)
		call setbolckData1(symeye_one%block,maxdim(1),maxdim(2),degeneracy)
		return
	end function
	type(SymTensor) function symeye_one_dim(dimen)
		type(SymDimension),intent(in)::dimen
		integer::maxdim(2)
		real*4::QN(2)
		integer,allocatable::degeneracy(:,:)
		call SymDiagonalTest(dimen.sub.1,dimen.sub.2)
		QN(1)=SymQNum(dimen,1)
		QN(2)=SymQNum(dimen,2)
		maxdim=2*QN+1
		allocate(degeneracy(2,maxval(maxdim)))
		call SymDimDegeneracy4(degeneracy(1,:),dimen,1)
		call SymDimDegeneracy4(degeneracy(2,:),dimen,2)
		call allocateSymTensor(symeye_one_dim,dimen)
		symeye_one_dim%totalblock=minval(maxdim)
		call setbolckData1(symeye_one_dim%block,maxdim(1),maxdim(2),degeneracy)
		return
	end function
	subroutine setbolckData1(block,LDR,LDC,degeneracy)!the SymTensor is a matrix, set the element on its diagonal,identity
		integer,intent(in)::LDR,LDC,degeneracy(:,:)
		type(Tensor),intent(inout)::block(LDR,LDC)
		integer::i
		if(size(degeneracy,2).ne.LDR) then
			write(*,*)"ERROR in setbolckData1"
			stop
		end if
		do i=1,LDR
			block(i,i)=eye(degeneracy(1,i),degeneracy(2,i))
		end do
		return
	end subroutine
	subroutine setbolckData2(block,LDR,LDC,degeneracy,blockdata)!the SymTensor is a matrix, set the element on its diagonal
		integer,intent(in)::LDR,LDC,degeneracy(:,:)
		type(Tensor),intent(inout)::block(LDR,LDC)
		type(Tensor),intent(in)::blockdata(:)
		integer::i
		if(size(degeneracy,2).ne.LDR) then
			write(*,*)"ERROR in setbolckData1"
			stop
		end if
		do i=1,LDR
			if((blockdata(i).dim.1).ne.degeneracy(1,i) .or. &
					(blockdata(i).dim.2).ne.degeneracy(2,i)) then
				write(*,*)"ERROR in setbolckData2"
				stop
			end if
			block(i,i)=blockdata(i)
		end do
		return
	end subroutine
	subroutine setbolckData3(block,LDR,LDC,ith,degeneracy,blockdata)!the SymTensor is a matrix, set the element on its diagonal
		integer,intent(in)::LDR,LDC,degeneracy(:),ith
		type(Tensor),intent(inout)::block(LDR,LDC)
		type(Tensor),intent(in)::blockdata
		if((blockdata.dim.1).ne.degeneracy(1) .or. &
				(blockdata.dim.2).ne.degeneracy(2)) then
			write(*,*)"ERROR in setbolckData3"
			stop
		end if
		block(ith,ith)=blockdata
		return
	end subroutine
	subroutine setbolckData4(block,LDR,LDC,degeneracy,blockdata)!the SymTensor is a matrix, set the element on its diagonal
		integer,intent(in)::LDR,LDC,degeneracy(:,:)
		type(Tensor),intent(inout)::block(LDR,LDC)
		type(Tensor),intent(in)::blockdata
		integer::i
		if(size(degeneracy,2).ne.LDR) then
			write(*,*)"ERROR in setbolckData1"
			stop
		end if
		do i=1,LDR
			if((blockdata.dim.1).ne.degeneracy(1,i) .or. &
					(blockdata.dim.2).ne.degeneracy(2,i)) then
				write(*,*)"ERROR in setbolckData4"
				stop
			end if
			block(i,i)=blockdata
		end do
		return
	end subroutine
!*****************  maxElement  *****************
	real*8 function SymmaxElement(T,block)
		type(SymTensor),intent(in) :: T
		logical,optional,intent(in)::block
		real*8::testmax,maxi
		integer::i
		logical::flag
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(maxElement)"
			stop
		end if
		flag=.false.
		do i=1,T%TotalData
			if(getFlag(T%block(i)))then!if1
				maxi=maxElement(T%block(i))
				if(flag)then!if2
					if(maxi.gt.testmax)then
						testmax=maxi
					end if
				else
					testmax=maxi
					flag=.true.
				end if!if2
			end if!if1
		end do
		if(present(block).and.block)then
			SymmaxElement=testmax
		else
			SymmaxElement=max(testmax,0d0)
		end if
		return
	end function
!*****************  maxRealElement  *****************
	real*8 function SymmaxRealElement(T,block)
		type(SymTensor),intent(in) :: T
		logical,optional,intent(in)::block
		real*8::testmax,maxi
		integer::i
		logical::flag
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(maxRealElement)"
			stop
		end if
		flag=.false.
		do i=1,T%TotalData
			if(getFlag(T%block(i)))then!if1
				maxi=maxRealElement(T%block(i))
				if(flag)then!if2
					if(maxi.gt.testmax)then
						testmax=maxi
					end if
				else
					testmax=maxi
					flag=.true.
				end if!if2
			end if!if1
		end do
		if(present(block).and.block)then
			SymmaxRealElement=testmax
		else
			SymmaxRealElement=max(testmax,0d0)
		end if
		return
	end function
!*****************  minElement  *****************
	real*8 function SymminElement(T,block)
		type(SymTensor),intent(in) :: T
		logical,optional,intent(in)::block
		real*8::testmin,mini
		integer::i
		logical::flag
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(minElement)"
			stop
		end if
		flag=.false.
		do i=1,T%TotalData
			if(getFlag(T%block(i)))then!if1
				mini=minElement(T%block(i))
				if(flag)then!if2
					if(mini.lt.testmin)then
						testmin=mini
					end if
				else
					testmin=mini
					flag=.true.
				end if!if2
			end if!if1
		end do
		if(present(block).and.block)then
			SymminElement=testmin
		else
			SymminElement=max(testmin,0d0)
		end if
		return
	end function
!*****************  minRealElement  *****************
	real*8 function SymminRealElement(T,block)
		type(SymTensor),intent(in) :: T
		logical,optional,intent(in)::block
		real*8::testmin,mini
		integer::i
		logical::flag
		if(.not.T%flag)then
			write(*,*)"There is no data in the Tensor,(minRealElement)"
			stop
		end if
		flag=.false.
		do i=1,T%TotalData
			if(getFlag(T%block(i)))then!if1
				mini=minRealElement(T%block(i))
				if(flag)then!if2
					if(mini.lt.testmin)then
						testmin=mini
					end if
				else
					testmin=mini
					flag=.true.
				end if!if2
			end if!if1
		end do
		if(present(block).and.block)then
			SymminRealElement=testmin
		else
			SymminRealElement=max(testmin,0d0)
		end if
		return
	end function	

	type(SymTensor) function Symexpm(H)result(expm)!input H shoulf be a matrix
		type(SymTensor),intent(in)::H
		integer::hdim1,hdim2
		if(SymgetRank(H).ne.2) then
			write(*,*)"ERROR in expm"
			write(*,*)"input Tensor should be a matrix"
			stop
		end if
		hdim1=H.dim.1
		hdim2=H.dim.2
		if(hdim1.ne.hdim2) then
			write(*,*)"ERROR in expm"
			stop
		end if
		call allocateSymTensor(expm,H%TenDim)
		expm%TotalBlock=H%TotalBlock
		call expm_Data(expm%block,hdim1,H%block)
		return
	end function
	subroutine expm_Data(expmHData,LDH,HData)
		integer,intent(in)::LDH
		type(Tensor),intent(in)::HData(LDH,LDH)
		type(Tensor),intent(inout)::expmHData(LDH,LDH)
		integer::i
		do i=1,LDH
			if(getFlag(HData(i,i)))then
				expmHData(i,i)=expm(HData(i,i))
			end if
		end do
		return
	end subroutine
		
	subroutine SVDtest()
		type(SymDimensionName)::dimenname,dimenname2
		type(SymDimension)::dimen
		type(SymTensor)::ST1,ST2,ST3,ST4,ST5,SU,Ss,SV,fuse1,fuse2,fuse3,fuse4,fuse5,SH
		type(Tensor)::T1,Sx,Sy,Sz,T2,T3,T4,T5,T6,H,U,S,V
		real*4::Qnum(3),time1,time2,S1,S2
		integer::rule(3),i,inde(3),NS1,NS2,Ncut
		logical::flag
		integer,allocatable::Adim(:),degeneracy(:,:)
		S1=0.5
		S2=S1+0.5
		NS1=2*S1+1
		NS2=2*S2+1
		allocate(degeneracy(3,NS2))
		degeneracy(1,:)=2
		degeneracy(2,:)=1
		degeneracy(3,:)=2
		ST1=Symgenerate_U1((/S2,0.5,S1/),degeneracy,(/1,1,-1/))
		ST2=Symgenerate_U1((/S1,0.5,S2/),degeneracy,(/1,1,-1/))
		T1=ST1
		T2=ST2
		Ncut=T1.dim.3
		
		ST3=ST1*ST2![1.,0.5,0.5',-1.']
		T3=T1*T2
		
		call pauli_matrix(Sx,Sy,Sz,1d0)
		H=(Sx.xx.Sx)+(Sy.xx.Sy)+(Sz.xx.Sz)
		H=H.c.1
		H=H.c.2
		H=expm(dcmplx(-1d-2,-1d-2)*H)
		H=.dc.H
		SH=Tensor2U1Ten(H,(/0.5,0.5,0.5,0.5/),(/-1,-1,1,1/),.true.)![-0.5,-0.5,0.5,0.5]
		
		ST3=SymTenproduct(ST3,(/2,3/),SH,(/1,2/))![1,-1,0.5,0.5]
		T3=Tenproduct(T3,(/2,3/),H,(/1,2/))
		T3=T3.pb.2![1,0.5,0.5,-1]
		
		fuse1=SymfuseTensor(ST3,1,3,1)![-1,-0.5,1.5]
		ST3=SymTenproduct(fuse1,(/1,2/),ST3,(/1,3/))![1.5,-1,0.5]
		fuse2=SymfuseTensor(ST3,2,3,-1)![1,-0.5,-1.5]
		ST3=SymTenproduct(ST3,(/2,3/),fuse2,(/1,2/))![1.5,-1.5]
		T3=ST3
		call symSVDcutoff(ST3,SU,Ss,SV,s1,degeneracy(1,:))		
		ST4=SU*Ss*SV![1.5,-1.5]
		
		T4=ST4
		call TMprint(T4-T3,5)
		call SymTMprint(ST4-ST3,5)
		stop
		call reverseRule(fuse2)![-1,0.5,1.5]
		ST4=SymTenproduct(ST4,(/2/),fuse2,(/3/))![1.5,-1,0.5]
		call reverseRule(fuse1)![1,0.5,-1.5]
		ST4=fuse1*ST4![1,0.5,-1,0.5]
		ST4=ST4.pb.3![1,0.5,0.5,-1]
		T5=ST4
		T5=T5.c.1
		T5=T5.c.2
		
		T3=T3.c.1
		T3=T3.c.2
		call SVDcutoff(T3,U,s,V,Ncut)		
		T4=U*eye(s)*V
		call TMprint(T4-T5,3)
	end subroutine
	subroutine testl()
		type(SymDimensionName)::dimenname,dimenname2
		type(SymDimension)::dimen
		type(SymTensor)::ST1,ST2,ST3,ST4,ST5,SU,Ss,SV,fuse1,fuse2,fuse3,fuse4,fuse5,SH
		type(Tensor)::T1,Sx,Sy,Sz,T2,T3,T4,T5,T6,H,U,S,V
		real*4::Qnum(3),time1,time2,S1,S2
		integer::rule(3),i,inde(3),NS1,NS2,Ncut
		logical::flag
		integer,allocatable::Adim(:),degeneracy(:,:)
		complex*16::rn
		real*8::num
		call pauli_matrix(Sx,Sy,Sz,0.5d0)
		H=(Sx.xx.Sx)+(Sy.xx.Sy)+(Sz.xx.Sz) ![1,1',2,2']
		SH=Tensor2U1Ten(H,(/0.5,0.5,0.5,0.5/),(/1,1,-1,-1/))
		H=H.c.1
		H=H.c.2
		H=expm(dcmplx(1,-1)*H)
		fuse1=SymfuseTensor(SH,1,2,1)![-0.5,-0.5,1]
		SH=SymTenproduct(fuse1,(/1,2/),SH,(/1,2/))![1,-0.5,0.5]
		fuse2=SymfuseTensor(SH,2,3,-1)![0.5,0.5,-1]
		SH=SymTenproduct(SH,(/2,3/),fuse2,(/1,2/))![1,-1]
		SH=Symexpm(dcmplx(1,-1)*SH)![1,-1]
		call reverseRule(fuse1)![0.5,0.5,-1]
		call reverseRule(fuse2)![-0.5,-0.5,1]
		SH=SymTenproduct(SH,(/2/),fuse2,(/3/))![1,-0.5,-0.5]
		SH=SymTenproduct(fuse1,(/3/),SH,(/1/))![0.5,0.5,-0.5,-0.5]
		call SymTMprint(SH,3)
		call TMprint(.dc.H,3)
		stop
		
		S1=0.5
		S2=S1+0.5
		NS1=2*S1+1
		NS2=2*S2+1
		allocate(degeneracy(3,NS2))
		degeneracy(1,:)=1
		degeneracy(2,:)=1
		degeneracy(3,:)=1
		ST1=Symgenerate_U1((/S2,0.5,S1/),degeneracy,(/1,1,-1/),(/-1d0,1d0/))
		ST2=ST1+ST1
		call SymTMprint(ST2)
		num=SymmaxElement(ST2)
		write(*,*)num
		num=SymminElement(ST2)
		write(*,*)num
		num=SymmaxRealElement(ST2)
		write(*,*)num
		num=SymminRealElement(ST2)
		write(*,*)num
		
	end subroutine
end module


