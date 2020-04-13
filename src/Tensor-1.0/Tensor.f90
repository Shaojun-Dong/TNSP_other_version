!************************************************************
!************* START OF Tensor **********************
!************************************************************
!*****************     worning    *******************
!			1.When calling lapack and blas subroutine,The matrix is input as 
!		one dimension vector.It gose no wrong in fortran90,but I do not 
!		test other version of fortran.
!			In fortran,the data store in memory as:A_{1,1}->A_{2,1}->A_{3,1}
!		->...A{M,1}->A_{1,2}->A_{2,2}->..->A_{M,2}->..->A_{M-1,N}->A_{M,N}.
!		But in matlab or C,it is as:A{0,0}->A_{0,1}->A_{0,2}->..A{0,N}
!		->A{1,0}->..->A{M,N-1}->A_{M,N}
!			2.A is a allocatable array and B is array too,when use A=B before
!		allocating storage space for A,A will get no data.This happen in 
!		fortran90 but do not in fortran77.The interface assignments is to  
!		solve this program.
!			3.A+B=C,where A is complex*16,B is real*8,and C is complex*16.It goes
!		no wrong in fortran90.This code is in add_DTen1 and add_DTen2
!			4.do not write a function that return a Tensorlink or Tensorlist,as the
!		data in the link or list will not be deallocate after calling the function.
!		one could write a surbroutine instead
!			5.The code is fortran90 version.
!*****************     note      *******************
!			1.The product of Tensor only allow rank<=2,if larger than 2,
!		use .p. and contract(or .c.) to reshape it.
!			2.permutation(.p.)  for rank<=3 is faster then the one with input a vector
!		,if larger then 3,use contract(or .c.) to reshape it to rank=3 or rank=2,and
!		then do the permutation,and at last decompose the Tensor back to the original
!		rank.One could reshape to any dimension through (.p.vectot),but it is slow
!		than (.p.number).
!			3.to compile the code one should link the files to lapack and blas.
!			4.svd ,eng and expm could be more faster by define the input matrix of the 
!		lapack subroutine as a vector of the output Tensor,becase the data is already
!		store as a vector,when copy to a matrix,it waste time.This save the time of 
!		assignment of the Tensor
!			5.These files call arpack lapack and blas . only eigen use arpack .
!		gfortran name.f90 -o name -larpack -llapack -lrefblas 
!			6.Send email to tyrants@qq.com to report any bugs.
!***********************************************************
!		bugs do not fix yet:
!			permutation:most sitiation,it gose right
!			expm:if the input matrix is not hermitian,gose wrong			
!
!		new function:Tensorlist,I just do a little test on it,2014.3.13
!						push_back , .i.vec , = ,link_i, are ok
module Tensor_complex
	use Tensor_real
	implicit none
!****************************************************	
!*********		define data type	***************
	type  Tensor
		integer,private :: rank!length of Dimension.TenDim%Dimsize
		type(Dimension),private	:: TenDim!dimenison
		complex*16,allocatable,private :: Tensor_data(:)!The data store in one dimenson array
		integer,private :: totalData!length of Tensor_data
		logical*4,private :: flag=.false. !if flag=true ,it means there are data in Tensor
	end type Tensor
! 	Tensorlink and Tensornode used as vector
!	Tensorlink store the head and the last element of Tensornode
!	Tensornode is used to store Tensor
!	indice(:) is used as index such as H_{1,2},H_{1,2,3} and so on
	type Tensornode
		integer,allocatable:: indice(:)
		type(Tensor) ::Ten
		type(Tensornode),private,pointer :: next
	end type Tensornode
		
	type Tensorlink
		integer,private :: length=0
		type(Tensornode),pointer,private  :: head
		type(Tensornode),pointer,private  :: Tend
	end type Tensorlink	
!	-----------
	type Tensorlist
		integer,private :: length=0
		type(listnode),private,pointer :: head
		type(listnode),pointer,private  :: Lend
	end type Tensorlist
	
	type listnode
		integer,allocatable:: indice(:)
		type(Tensorlink) ::link
		type(listnode),private,pointer :: next
	end type listnode
	
!*********		End define		*********************


!****************************************************
!*********		define private data	***************	
	real*8,private :: cputime(2)
	complex*16,parameter,private :: II=(0,1),EE=(1,0)! II**2=-1
	real*8,parameter,private :: large_number_warning=1d80!use for checking
	real*8,parameter,private :: zero_num=1d-20!
	logical,private,save::check_flag=.false.!use for checking
	complex*16,private::ZDOTU!dot product
	complex*16,private::ZDOTC!dot product,conjugating the first vector
	real*8,private::DZNRM2!norm :sqrt(A**H * A)
	integer,save,private::max_ncv=500!use in funciton of eigen2 and eigen_max_min2 ncv<=max_ncv
	External::ZDOTU,DZNRM2,ZDOTC!blas's function
	private::storeTenData
	private::addressToIndes
	private::addressToIndes2
	private::IndesToaddress
	private::compose_decompse_subroutine1,compose_decompse_subroutine2
	private::set
!*********		End define	***************		
		
	interface TMprint
		module procedure TMprint1
		module procedure TMprint2
	end interface
	interface Tprint
		module procedure Tprint1
		module procedure Tprint2
	end interface
	
	interface sqrt
		module procedure sqrtTen
	end interface
	
	interface eig!eigenvalue
		module procedure eng
		module procedure eng2
	end interface
	
	interface eigen
!eigenvalue near some value,say s0
!ncv is the Number of Lanczos vectors to be generated.2 <= NCV-num_of_eig and NCV <= N		
!tol:Stopping criteria
		module procedure eigen1!(T,s0,num_of_eig,ncv,tol,outtype)
		module procedure eigen2!(T,s0,num_of_eig,outtype)
		module procedure eigen3!(T,s0,num_of_eig,),outtype=1
!	if outtype=1 eigenvalues
!	if outtype=2 eigenvalues and eigenvectors
!	if outtype=3 eigenvalues,but using |T-s0|^2 to find the eigenvalue near s0
!	if outtype=4 eigenvalues and eigenvectors,but using |T-s0|^2 to find the eigenvalue near s0		
		
		module procedure eigen_max_min1!(T,WHICH,num_of_eig,ncv,tol,outvector)
		module procedure eigen_max_min2!(T,WHICH,num_of_eig,outvector)
		module procedure eigen_max_min3!(T,WHICH,num_of_eig),outvector=.false.
!				WHICH   Character*2.  (INPUT)
!          'LM' -> want the NEV eigenvalues of largest magnitude.
!          'SM' -> want the NEV eigenvalues of smallest magnitude.
!          'LR' -> want the NEV eigenvalues of largest real part.
!          'SR' -> want the NEV eigenvalues of smallest real part.
!          'LI' -> want the NEV eigenvalues of largest imaginary part.
!          'SI' -> want the NEV eigenvalues of smallest imaginary part.
	end interface
	
	interface set!create a Tensor,giving dimension and Tensor_Data
					  !buildTen can do all the thing of set so I define it as private
		module procedure set1 
		module procedure set2 
	end interface
	
!	copy the array to the T%Tensor_Data,array can be 1,2 or 3 dimension
	interface storeTenData
		module procedure storeTenData_dim1
		module procedure storeTenData_dim2
		module procedure storeTenData_dim3
	end interface
	
	interface	buildTen!build a Tensor
! 	note when input a vector in to the tensor of rank=2
! for example buildTen((/3,2/),(/a,b,c,d,e,f/)),then the Tensor return is
!	/ a d \
!	| b e |
!	\ c f /
		module procedure buildTen1
		module procedure buildTen2
		module procedure buildTen3
		module procedure buildTen4
		module procedure buildTenreal4
		module procedure buildTen5
		module procedure buildTenreal5
		module procedure buildTen6
		module procedure buildTenreal6
		module procedure buildTen7
		module procedure buildTen8
		module procedure buildTen9
	end interface
	
	interface modify
		module procedure modifylink
		module procedure modifylist
		module procedure modifyTen_val
		module procedure modifyTen_dat
		module procedure modifyTen_dat2
		module procedure modifyTen_dat3
	end interface
!************************************************************	
!operateTensor(A,B,oA,oB,productflag):
!	oA,oB is a matrix or vector,the first element of every row specify
!		1:contracts
!		2:decompose
!		3:decomposeAll
!	other element of every row is input parameter of the function
!		at last A*B .when flag=1,result store in A with B no change when productflag=2,
!	 the result will store in B
!		this subroutine run for big inoutA and B,do as less as possible operation on the
!	 Tensor_Data of inoutA and B.
!************************************************************	
	interface operateTensor
		module procedure operateTensor1
		module procedure operateTensor2
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
	interface dimOperation
		module procedure dimOperations
		module procedure dimOperation1
	end interface
	
	
	interface resetdim
		module procedure resetdim1
		module procedure resetdim2
	end interface
!************************************************************		
	
	interface eye
!  return a matrix of diag[s1,s2,s3,...s_{max},0,0...]	,a m*n matrix
		module procedure eye_com!input s(:) is complex*16
		module procedure eye_real!input s(:) is real*8
		module procedure one_com!identity matrix
		module procedure one_Ten!return a identity Tensor 
		module procedure eye_Ten1!input Tensor,output diag[T.i.1 , T.i.2 , ... , , T.i.n],rank of T should be 1
		module procedure eye_Ten2!input Tensor,output a matrix of m*n: diag[T.i.1 , T.i.2 , ... , , T.i.n],rank of T should be 1
	end interface
!***********   interface for Tensorlink  and Tesorlist *************************	
	interface push_back
		module procedure push_backTen
		module procedure push_backnode
		module procedure push_backI
		module procedure push_backnum
		module procedure push_backnumI
		
		module procedure push_backlink
		module procedure push_backlinknode
		module procedure push_backlinkI
		
		
	end interface
	
	interface push_forward
		module procedure push_fo
		module procedure push_fonode
		module procedure push_foI
		module procedure push_fonum
		module procedure push_fonumI
		
		module procedure push_folink
		module procedure push_folinknode
		module procedure push_folinkI
	end interface
	
	interface node_i
		module procedure node_i_int
		module procedure node_i_vec
		
		module procedure listnode_i_int
		module procedure listnode_i_vec
	end interface
	
	interface link_i!Do not use Link=list.i.1,because the nodes in Link are pointer,doing so will create rubbish data in memory
		module procedure Li!output the i Tensorlink in the Tensorlist
		module procedure LiI!output the index Tensorlink in the Tensorlist,index is a vector
	end interface
!*******************   permutation   ****************************	
	interface operator(.p.)
		module procedure permute_rank2 !permute the tensor whose rank is 2
		module procedure permute_rank3 !permute the tensor whose rank is 3
		module procedure permutation	 !permute the tensor of any rank,give the new order of the dimension
												 !		If operate on a big Tensor,use permute_rank2 or permute_rank3,
												 !	 they are faster.if rank>3,use contract to reshape
	end interface
	interface operator(.pf.)
		module procedure permutefo!permute the inde index to the first
										  !T_{1,2,3,..,i,..,n},permutefo(T,i)=_{i,1,2,3,..,i-1,i+1,..,n}
										  ! note the output will decompose all the dimension
	end interface
	interface operator(.pi.)!T_{1,2,3,..,i,..,n},permuteInde(T,i)=_{2,3,..,i,1,i+1,..,n}
									! note the output will decompose all the dimension
		module procedure permuteInde
	end interface
		
	interface operator(+)
		module procedure add!the element of one Tensor push the other one
		module procedure add_DTen1!T1+T2,T1 is complex,T2 is real
		module procedure add_DTen2!T1+T2,T2 is complex,T1 is real
		module procedure add_num!one Tensor push an identity Tensor 
		module procedure add_num_!num+ Tensor
		module procedure connectlink!combine two Tensor link
		module procedure connectlist!combine two list
	end interface
	
!**************** contract   ****************
!		combine two index of the Tensor,which is index1 index1+1,..,index2-1,index2
!		if index2 larger than rank,the function will contract the index of index1 -> rank
	interface operator(.c.)
		module procedure contract!combine two index of the Tensor,which is con_index and con_index+1
		module procedure contractsv!combine two index of the Tensor,input vector specify the to be combined index
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
		module procedure compose_decompse
	end interface
!************************************************************		
! decompose the de_index index of the Tensor into n(1),n(2)
!		for example the de_index index is (1,2,3,4,..inde,inde+1,...rank)
!		(1,2,3,4,..inde,inde+1,...rank)-->(1,2,3,4,..inde),(inde+1,...rank)		
!		if inde larger than rank ,the function will return no change	
	interface operator(.dc.)
		module procedure decompose1
		module procedure decomposeAll
		module procedure decomposev!decompose index of the Tensor,input vector specify the de_index and index
	end interface
	
	interface operator(-)
		module procedure minus!the element of one Tensor minue the other one
		module procedure minus_num!one Tensor minue an identity Tensor 
		module procedure minus_num_!num - Tensor
	end interface
!*******************************************************************************
!		productTen is the old version
!		the function productTensor add at 2013.10.28
	interface operator(*)
		!module procedure productTen!Matrixe product,two Tensor intput should be of 1 or 2 dimension
		module procedure productTensor!A*B,product the last index of A and first index of B
		module procedure multiply_number!A Tensor times a number,T*num
		module procedure multiply_number_!num*T
		module procedure multiply_real!A Tensor times a real*8 number
		module procedure multiply_real_!a real*8 number times a Tensor 
		module procedure multiply_real4!A Tensor times a real*4 number
	end interface
	
	interface operator(/)
		module procedure divide!the element of one Tensor divide the element of other one
		module procedure divide_number!A Tensor divide a number of type complex*16
		module procedure divide_numberReal!A Tensor divide a number of type real*8
		module procedure divide_numberReal4!A Tensor divide a number of type real*4
	end interface
	
	interface operator(.x.)
		module procedure dot! dot product conjugating the first vector,The Tensor will be regard as a vector
	end interface
	
	interface operator(.xx.)
		module procedure directProduct! direct Product [m1,n1] * [m2,n2] => [m1,m2,n1,n2]
												!  or direct Product return a matrix [n1] * [n2] => [ n1,n2 ]
	end interface
		
	interface operator(.mxx.)
		module procedure directProductM! direct Product return a matrix [m1,n1] * [m2,n2] => [(m1*m2),(n1*n2)]
												!  or direct Product return a matrix [n1] * [n2] => [ (n1*n2) ]
	end interface
		
	interface operator(.su.)
		module procedure direct_sum! direct sum
	end interface
		
	interface operator(.msu.)
		module procedure direct_sumM! direct sum return a matrix
	end interface
		
	interface operator(.elex.)
		module procedure multiply_Tensor!the element of one Tensor times the element of other one
	end interface
	
	interface operator(.numx.)
		module procedure TensorProduct1Dim!two one-dimension Tensors product,return a number of type complex*16
	end interface
	
	interface operator(.H.)
		module procedure Htranspose! conjugateTranspose! one or two dimension case
											! If input rank=1,the result is the same as conjugate
	end interface
	interface operator(.con.)
		module procedure conjugate! conjugate
	end interface
	interface operator(.equ.)
		module procedure equal_of_Tensor!If two Tensor are equal
	end interface
	
	interface operator(.subDim.)
		module procedure getTenSubDim!output the ith dimension of the Tensor,return type(Dimension)
		module procedure getTenDim!output the Dimension of Tensor,return type(Dimension)
											!It can use to commit the assignment to a vector
											!vector=type(Dimension),if the vector is allocatable
	end interface

	interface operator(.i.)
		module procedure Ti!output the i Tensor in the Tensorlink
		module procedure Element!output the data in the Tensor of Tim=(i,j)
		module procedure Element2!output the data in the Tensor of inde
		module procedure TiI!output the index Tensor in the Tensorlink,index is a vector
		module procedure Liv
		module procedure Liiv
	end interface
	
	
	
	interface operator(.dim.)
		module procedure getTenDim_i!output the i dimension of the Tensor,output an integer
	end interface
		
	interface assignment(=)
		module procedure copy!T1=T2 ,both T1 and T2 are Tensor
		module procedure Ten_DTen!T1=T2 ,T1 is Tensor and T2 is DTensor
		module procedure copy_dim1!Vec=T ,Vec is a vector and T are Tensor,vec=T%Tensor_Data,used for any case
		module procedure copy_dim2!Mat=T ,Mat is a matrix and T are Tensor,Mat=T%Tensor_Data,used only rank=2
		module procedure copy_dim3!Mat=T ,Mat is a rank=3 and T are Tensor,Mat=T%Tensor_Data,used only rank=3
		module procedure assignmentTenReal3!T=vec or matrix
		module procedure assignmentTenReal2
		module procedure assignmentTenReal1
		module procedure assignmentTen3!T=vec or matrix
		module procedure assignmentTen2
		module procedure assignmentTen1
		module procedure com_Ten!val=T,val is complex*16 ,and T is a Tensor,used only there is one element in T
		module procedure  TenLink!Tensor=Tensorlink,if every element(Tensor) in the link is a one-element Tensor
		
		module procedure copylink!link1=link2,both link1 and link2 are TensorLink,and the result ,link1 and link2 will not be the same link
		module procedure copylist!list1=list2,both list1 and list2 are TensorList
		!module procedure copylinkhead!link1=link2,both link1 and link2 will be the same link
	!if one want to write a function return a link,use copylinkhead as (=)
!     when using a subroutine that return a link,do not use result=link to output,because the data in link will not
!	empty after call the function.But one could use 
!		1. call copylinkhead(result,link) 
!		2. result=link
!			call cleanlink(link)
!	the first chose is better,it is faster
	end interface
!	other way to build a Tensor:T=buildTen(dimen,data),dimen is a vector or type(Dimension),data is any type
!**********************************************************
!	Other Function or Subroutine:
!
!		generate:generate a Tensor with random number
!
!		cleanTensor: clean all the data in type(Tensor)
!
!		Print Tensor: 
!			Tprint(Tensor):print all the information of the Tensor
!			TDprint(Tensor):print Dimension
!			TMprint(Tesnor):Print as matrix 
!			Lprint(Tensorlink):print every indice of the Tensor in the link
!			LDprint(Tensorlink):print every dimension and indice of the Tensor in the link
!
!		get rank or totalData or flag
!			getRank
!			getTotalData
!			getflag
!
!		equal_of_Tensor:to judge if two Tensor are equak
!
!		contracts((T1,index1,index2)):combine two index of the Tensor,
!			which is index1 index1+1,..,index2-1,index2.if index2 larger 
!			than rank,the function will contract the index of index1 -> rank
!
!		decompose(T1,de_index,inde):decompose the de_index index of the Tensor 
!			into n(1),n(2) for example the de_index index is (1,2,3,4,..inde,inde+1,...rank).
!			(1,2,3,4,..inde,inde+1,...rank)-->(1,2,3,4,..inde),(inde+1,...rank)		
!			if inde larger than rank ,the function will return no change	
!
!		max or min element:
!			maxElement(T):return the max abs value element of the Tensor
!			maxRealElement(T):return the max real part element of the Tensor
!			minElement(T):return the min abs value element of the Tensor
!			minRealElement(T):return the min real part element of the Tensor
!
!		realTen(T):return the real part of the Tensor
!
!		imagTen(T):return the imag part of the Tensor
!
!		TensorProduct1Dim(T1,T2): Dot product,return a comple*16
!
!		subTen(T,inde):inde=[-1,inde_min,inde_max] output data(:,inde_min:inde_max)
!			or[-2,inde_min,inde_max],data(inde_min:inde_max,:)
!			or [-1,inde_row] [-2,inde_col],output row or col
!
!		SVD :
!			svdright(T,T1,T2,check_discard): T=U*s*V^T,return the (U*s) 
!				in T1 and V^T in T2 
!			svdleft(T,T1,T2,check_discard):	T=U*s*V^T ,return the U 
!				 in T1,and (s*V^T) in T2
!			Both in svdright and svdleft the data s(i) in the matrix s will 
!				be discard if s(i)<s(1)*check_discard_svd.if check_discard_svd<0
!				,all the data in s will be keep
!
!		linear equations
!			linequ:X=linequ(A,B),solve A*X=B,X could be a matrix
!			linequ_routine(A,B):on output,A will change and X is in B	
!
!		diagonal matrix:
!			eye_com(s,m,n):return a m*n diagonal matrix will diag(s(1),s(2),...s(slen),0,0....),
!				s(:) are complex*16
!			eye_real::return a m*n diagonal matrix will diag(s(1),s(2),...s(slen),0,0....),
!				s(:) are real*8
!			one_com(m,n)  :return a m*n identity matrix
!
!		pauli_matrix(Tx,Ty,Tz,num)
!
!		expm(H):
!			return a Tensor e^H,
!
!		eng(H,vec,val):
!			return a eigenvalues matrix of H,H=vec*val*vec^H
!
!		sqrt:
!			return sqrt(A),A is a matrix
!
!		resetdim:
!			reset the dimension of the Tensor
!
!		conbine
!			T1 :a [l,m,n,...] matrix
!			T2 :a [l,m,n,...] matrix
!			conbine(T1,T2):a [2,l,m,n,...] matrix
!			or 
!			T1 :a [l,m,n,...] matrix
!			T2 :a [m,n,...] matrix
!			conbine(T1,T2):a [l+1,m,n,...] matrix
!		conbinerow(T1,T2)
!			T1 :a [...,l,m,n] matrix
!			T2 :a [...,l,m,n] matrix
!			conbine(T1,T2):a [...,l,m,n,2] matrix
!			or 
!			T1 :a [...,m,n,l] matrix
!			T2 :a [...,m,n] matrix
!			conbine(T1,T2):a [...,m,n,l+1] matrix
!			
!		operator value:
!			value(F,Tr):return the value of <phi|F|phi>,where F is an operator,
!				F=<I',s1',s2',..|F|I,s1,s2,..>
!			norm2(Tr_):return  <phi|phi>
!			norm(Tr):return  sqrt(<phi|phi>)
!
!		Tensor1(dimen,num):return a Tensor with all the elememt num,dimen is dimension of
!			the Tensor
!		
!		Ten_Ten(T1_,T2_,i1,i2,ifDecompose):contract The i1 index of T1 and the i2 index of t2
!			T1: (D1,D2,..D_i1,...),T2 :(E1,E2,...,E_i2,....),then the result will be 
!			T:(D1,D2,..D_{i1-1},D_{i1+1},...,D_{rank1},E1,E2,...,E_{i2-1},E_{i2+1},...,E_{rank2})	
!			if if Decompose=1, T will be Decompose
!
!		Tensor index
!			TenIndice(h,inde,indice):return the index of the inde Tensor.T_{1,1}->T_{1,2}->T_{1,3}
!				->T_{2,1}->T_{2,2},if input inde=2 => indice=[1,2].on entry the size of indice should
!				 be equal to the one in h
!			TenIndiceLog(h,inde,indice):if there no index in h,return .false..other while return inde
!				of the Tensor
!			lenOfIndicelist(h,inde,indice):return the index of the inde Tensorlink.T_{1,1}->T_{1,2}->T_{1,3}
!				->T_{2,1}->T_{2,2},if input inde=2 => indice=[1,2].on entry the size of indice should
!				 be equal to the one in h
!			linkIndiceLog(h,inde,indice):if there no index in h,return .false..other while return inde
!				of the Tensorlink
!
!		pointer subroutine
!			ifnextnode:if the next node exict,p point to the next node,or p remain nochange
!			headnode:make the pointer points to the head of the link
!			nextnode:make the pointer points to the next node
!			endnode:make the pointer points to the end of the link
!			addnode:add a empty node to the link(then use a pointer to point to it)
!
!			ifnextlistnode:if the next node exict,p point to the next node,or p remain nochange(Tensorlist subroutine)
!			headlist:make the pointer points to the head of the list(Tensorlist subroutine)
!			nextlistnode:make the pointer points to the next node(Tensorlist subroutine)
!			endlist:make the pointer points to the end of the list(Tensorlist subroutine)
!			addlist:add a empty node to the link(then use a pointer to point to it)(Tensorlist subroutine)
!
!		checklength(link):check and output length of the link,check if the length of the link is equal 
!				to the value in head of the link
!		checklistlength(list):check and output length of the list,check if the length of the list is equal 
!				to the value in head of the list
!
!		linklength(link):output length of the link
!		listlength(list):output length of the list
!
!		copylinkhead:copy the head of a link
!
!		modify(h,T,inde):modify the inde element of the link.The Tensor will be modify
!
!		cleanlink(h):clean the link h
!		cleanlist(h):clean the list h
!
!		deletelink(link_in,inde1,inde2):
!			link   :	T1->T2->...->T_{inde1-1}->T_inde1->T_{inde1+1}->....->T_{inde2-1}->T_inde2->T_{inde2+1}...->TN 
!		  result:		T1->T2->...->T_{inde1-1}->T_{inde1+1}....->T_{inde2-1}->->T_{inde2+1}...->TN 
!		  delete the inde1 to inde2 Tensor in the link,note that the deleted Tensor are include the Tensors of inde1 and inde2
!		deletelist(list_in,inde1,inde2):(Tensorlist subroutine)
!			list   :	T1->T2->...->T_{inde1-1}->T_inde1->T_{inde1+1}->....->T_{inde2-1}->T_inde2->T_{inde2+1}...->TN 
!		  result:		T1->T2->...->T_{inde1-1}->T_{inde1+1}....->T_{inde2-1}->->T_{inde2+1}...->TN 
!		  delete the inde1 to inde2 Tensorlink in the list,note that the deleted Tensorlink are include the Tensorlinks of inde1 and inde2
!
!		init_random_seed
!
!	Other operator ,see the note in file Dimension.f90
!
!			assignment(=):
!				assignment of the Dimension type;
!				Dimension=vector 
!					or
!				Dimension1=Dimension2 
!					or
!				vector=Dimension
!
!			operator(+):
!				Dimension1+Dimension2:[2,3]+[4,5,6]->[2,3,4,5,6] 
!				Dimension+vector 
!					 
!			operator(.sub.)
!				get the ith of dimension renturn dimenison tpye
!
!			operator(.i.)
!				get the ith of dimension data
!
!			operator(.equ.)
!				equal_of_array:If two array of one dimension are equal
!				equal_of_dim:If two type(dimension) are equal
!
!			delta:delta function
!
!			interface assignments:
!					A is a allocatable array and B is array too,when use
!				A=B before allocating storage space for A,A will get
!			 	no data.This happen fortran90 but do not in fortran77.
!				assignments is used for A=B
!
!**********************************************************	
	contains
	
!***************** assignment *********************
	subroutine copy(T,T2)
		type(Tensor),intent(out) ::T
		type(Tensor),intent(in) :: T2
		if(T2%flag) then
			T%rank=T2%rank
			T%TenDim=T2%TenDim
			call storeTenData(T,T2%Tensor_Data)
			T%totalData=T2%totalData
			T%flag=T2%flag
		else
			T%flag=T2%flag
		endif
		return
	end subroutine
	subroutine com_Ten(val,T)
		complex*16,intent(out)::val
		type(Tensor),intent(in)::T
		if(T%totalData.eq.1) then
			val=T.i.1
		else
			write(*,*)"ERROR in assignment for Tensor to complex"
			stop
		end if
		return
	end subroutine
	subroutine Ten_DTen(T,DT)
		type(Tensor),intent(inout)::T
		type(DTensor),intent(in)::DT
		real*8,allocatable::TData(:)
		if(Dgetflag(DT)) then
			TData=DT
			T=TData
			T%rank=dgetRank(DT)
			T%TenDim=.subdim.DT
			T%totalData=DgetTotalData(DT)
			T%flag=Dgetflag(DT)
		else
			T%flag=Dgetflag(DT)
		endif
		return
	end subroutine
	subroutine assignmentTen1(T,Tensor_data)
		complex*16,intent(in)::Tensor_data(:)
		type(Tensor),intent(inout)::T
		call cleanTensor(T)
		T=buildTen(Tensor_data)
		return
	end subroutine
	subroutine assignmentTen2(T,Tensor_data)
		complex*16,intent(in)::Tensor_data(:,:)
		type(Tensor),intent(inout)::T
		call cleanTensor(T)
		T=buildTen(Tensor_data)
		return
	end subroutine
	subroutine assignmentTen3(T,Tensor_data)
		complex*16,intent(in)::Tensor_data(:,:,:)
		type(Tensor),intent(inout)::T
		call cleanTensor(T)
		T=buildTen(Tensor_data)
		return
	end subroutine
	subroutine assignmentTenReal1(T,Tensor_data)
		real*8,intent(in)::Tensor_data(:)
		type(Tensor),intent(inout)::T
		call cleanTensor(T)
		T=buildTen(Tensor_data)
		return
	end subroutine
	subroutine assignmentTenReal2(T,Tensor_data)
		real*8,intent(in)::Tensor_data(:,:)
		type(Tensor),intent(inout)::T
		call cleanTensor(T)
		T=buildTen(Tensor_data)
		return
	end subroutine
	subroutine assignmentTenReal3(T,Tensor_data)
		real*8,intent(in)::Tensor_data(:,:,:)
		type(Tensor),intent(inout)::T
		call cleanTensor(T)
		T=buildTen(Tensor_data)
		return
	end subroutine
		
!*************** cleanTensor  *****************
	subroutine cleanTensor(T)
		Type(Tensor),intent(inout)::T
		T%rank=0
		call cleanDimension(T%TenDim)
		if(allocated(T%Tensor_data)) then
			deallocate(T%Tensor_data)
		end if
		T%totalData=0
		T%flag=.false.
		return
	end subroutine
!copy the Tensor_Data to a vector	
	subroutine copy_dim1(Vec,T)
		complex*16,allocatable,intent(out) ::Vec(:)
		type(Tensor),intent(in) :: T
		integer::length
		length=T%TotalData
		if(allocated(Vec)) then
			deallocate(Vec)
		end if
		allocate(Vec(length))
		call ZCOPY(length,T%Tensor_Data,1,Vec,1)
		return
	end subroutine
!if rank=2,then copy the Tensor_Data to a mat
! mat is a matrix	
	subroutine copy_dim2(Mat,T)
		complex*16,allocatable,intent(out) ::Mat(:,:)
		type(Tensor),intent(in) :: T
		integer::m,n
		if(T%rank.eq.2) then
			m=T.dim.1
			n=T.dim.2
			if(allocated(Mat)) then
				deallocate(Mat)
			end if
			allocate(Mat(m,n))
			call ZCOPY(m*n,T%Tensor_Data,1,Mat,1)
		else
			write(*,*)"The Tensor is not a matrix"
			write(*,*)"Assignment ERROR,stop"
			stop
		end if
		return
	end subroutine
!if rank=3,then copy the Tensor_Data to a mat
!mat is a 3 dimension array	
	subroutine copy_dim3(Mat,T)
		complex*16,allocatable,intent(out) ::Mat(:,:,:)
		type(Tensor),intent(in) :: T
		integer::m,n,l
		if(T%rank.eq.3) then
			m=T.dim.1
			n=T.dim.2
			l=T.dim.3
			if(allocated(Mat)) then
				deallocate(Mat)
			end if
			allocate(Mat(m,n,l))
			call ZCOPY(m*n*l,T%Tensor_Data,1,Mat,1)
		else
			write(*,*)"Assignment ERROR,stop"
			stop
		end if
		return
	end subroutine	
!***********	storeTenData  **************
	subroutine storeTenData_dim1(T,inData)
		complex*16,intent(in)::inData(:)
		type(Tensor),intent(inout)::T
		integer::length
		length=size(inData)
		if(allocated(T%Tensor_Data)) then
			deallocate(T%Tensor_Data)
		end if
		allocate(T%Tensor_Data(length))
		call ZCOPY(length,inData,1,T%Tensor_Data,1)
		return
	end subroutine
	subroutine storeTenData_dim2(T,inData)
		complex*16,intent(in)::inData(:,:)
		type(Tensor),intent(inout)::T
		integer::length
		length=size(inData)
		if(allocated(T%Tensor_Data)) then
			deallocate(T%Tensor_Data)
		end if
		allocate(T%Tensor_Data(length))
		call ZCOPY(length,inData,1,T%Tensor_Data,1)
		return
	end subroutine
	subroutine storeTenData_dim3(T,inData)
		complex*16,intent(in)::inData(:,:,:)
		type(Tensor),intent(inout)::T
		integer::length
		length=size(inData)
		if(allocated(T%Tensor_Data)) then
			deallocate(T%Tensor_Data)
		end if
		allocate(T%Tensor_Data(length))
		call ZCOPY(length,inData,1,T%Tensor_Data,1)
		return
	end subroutine
	
!***************** buildTen *********************
!	build a Tensor
	type(Tensor) function buildTen1(rank,TenDim,Tensor_data)
		integer,intent(in)::rank
		type(Dimension),intent(in)::TenDim
		complex*16,intent(in)::Tensor_data(:)
		buildTen1%rank=rank
		buildTen1%TenDim=TenDim
		call storeTenData(buildTen1,Tensor_Data)
		buildTen1%totalData=size(Tensor_Data)
		buildTen1%flag=.true.
		return
	end function	
	type(Tensor) function buildTen2(rank,TenDim,Tensor_data)
		integer,intent(in)::rank
		integer,intent(in)::TenDim(:)
		complex*16,intent(in)::Tensor_data(:)
		buildTen2%rank=rank
		buildTen2%TenDim=TenDim
		call storeTenData(buildTen2,Tensor_Data)
		buildTen2%totalData=size(Tensor_Data)
		buildTen2%flag=.true.
		return
	end function		
	type(Tensor) function buildTen3(TenDim,Tensor_data)
		integer,intent(in)::TenDim(:)
		complex*16,intent(in)::Tensor_data(:)
		buildTen3%rank=size(TenDim)
		buildTen3%TenDim=TenDim
		call storeTenData(buildTen3,Tensor_Data)
		buildTen3%totalData=size(Tensor_Data)
		buildTen3%flag=.true.
		if(product(TenDim).ne.size(Tensor_Data)) then
			write(*,*)"ERROR in dimension and input data in function buildTen"
			write(*,*)"stop"
		end if
		return
	end function	
	type(Tensor) function buildTen4(Tensor_data)
		complex*16,intent(in)::Tensor_data(:)
		buildTen4%rank=1
		buildTen4%TenDim=(/size(Tensor_Data)/)
		call storeTenData(buildTen4,Tensor_Data)
		buildTen4%totalData=size(Tensor_Data)
		buildTen4%flag=.true.
		return
	end function
	type(Tensor) function buildTenreal4(Tensor_data)
		real*8,intent(in)::Tensor_data(:)
		complex*16,allocatable::TenData(:)
		allocate(TenData(size(Tensor_data)))
		TenData=Tensor_data
		buildTenreal4%rank=1
		buildTenreal4%TenDim=(/size(Tensor_Data)/)
		call storeTenData(buildTenreal4,TenData)
		buildTenreal4%totalData=size(Tensor_Data)
		buildTenreal4%flag=.true.
		return
	end function			
	type(Tensor) function buildTen5(Tensor_data)
		complex*16,intent(in)::Tensor_data(:,:)
		buildTen5%rank=2
		buildTen5%TenDim=(/size(Tensor_Data,1),size(Tensor_Data,2)/)
		call storeTenData(buildTen5,Tensor_Data)
		buildTen5%totalData=size(Tensor_Data)
		buildTen5%flag=.true.
		return
	end function		
	type(Tensor) function buildTenreal5(Tensor_data)
		real*8,intent(in)::Tensor_data(:,:)
		complex*16,allocatable::TenData(:,:)
		allocate(TenData(size(Tensor_data,1),size(Tensor_data,2)))
		TenData=Tensor_data
		buildTenreal5%rank=2
		buildTenreal5%TenDim=(/size(TenData,1),size(TenData,2)/)
		call storeTenData(buildTenreal5,TenData)
		buildTenreal5%totalData=size(TenData)
		buildTenreal5%flag=.true.
		return
	end function		
	type(Tensor) function buildTen6(Tensor_data)
		complex*16,intent(in)::Tensor_data(:,:,:)
		buildTen6%rank=3
		buildTen6%TenDim=(/size(Tensor_Data,1),size(Tensor_Data,2),size(Tensor_Data,3)/)
		call storeTenData(buildTen6,Tensor_Data)
		buildTen6%totalData=size(Tensor_Data)
		buildTen6%flag=.true.
		return
	end function		
	type(Tensor) function buildTenreal6(Tensor_data)
		real*8,intent(in)::Tensor_data(:,:,:)
		complex*16,allocatable::TenData(:,:,:)
		allocate(TenData(size(Tensor_data,1),size(Tensor_data,2),size(Tensor_data,3)))
		TenData=Tensor_data
		buildTenreal6%rank=3
		buildTenreal6%TenDim=(/size(TenData,1),size(TenData,2),size(TenData,3)/)
		call storeTenData(buildTenreal6,TenData)
		buildTenreal6%totalData=size(TenData)
		buildTenreal6%flag=.true.
		return
	end function		
	type(Tensor) function buildTen7(TenDim,Tensor_data)
		integer,intent(in)::TenDim(:)
		real*8,intent(in)::Tensor_data(:)
		complex*16,allocatable::TenData(:)
		allocate(TenData(size(Tensor_data)))
		TenData=Tensor_data
		buildTen7=buildTen3(TenDim,TenData)
		return
	end function	
	type(Tensor) function buildTen8(TenDim,Tensor_data)
		type(Dimension),intent(in)::TenDim
		integer,allocatable::Dimen(:)
		complex*16,intent(in)::Tensor_data(:)
		Dimen=TenDim
		buildTen8%rank=size(Dimen)
		buildTen8%TenDim=TenDim
		call storeTenData(buildTen8,Tensor_Data)
		buildTen8%totalData=size(Tensor_Data)
		buildTen8%flag=.true.
		
		if(product(Dimen).ne.size(Tensor_Data)) then
			write(*,*)"ERROR in dimension and input data in function buildTen"
			write(*,*)"stop"
		end if
		return
	end function	
	type(Tensor) function buildTen9(TenDim,Tensor_data)
		type(Dimension),intent(in)::TenDim
		integer,allocatable::Dimen(:)
		real*8,intent(in)::Tensor_data(:)
		complex*16,allocatable::TenData(:)
		allocate(TenData(size(Tensor_data)))
		TenData=Tensor_data
		buildTen9=buildTen8(TenDim,TenData)
		return
	end function

!********************* set *********************
!   	set the data in a Tensor 	 	c
	function set1(dim_,Ten_data) result (T)
		type(Tensor) :: T
		complex*16,intent(in) :: Ten_data(:)
		integer,intent(in)	:: dim_(:)
		integer :: total,i,rank,test_total
		type(Dimension) :: TenDim
		rank=size(dim_)
		total=size(Ten_data)
		test_total=product(dim_)
		if(test_total.ne.total) then
			write(*,*)"ERROR in set1,stop"
			stop
		end if
		TenDim=dim_
		call storeTenData(T,Ten_data)
		T%rank=rank
		T%totalData=total
		T%TenDim=TenDim
		T%flag=.true.
		return
	end function
		 
	function set2(TenDim,Ten_data) result (T)
		type(Tensor) :: T
		complex*16,intent(in) :: Ten_data(:)
		type(Dimension),intent(in)	:: TenDim
		integer :: total,i,rank
		rank=DimSize(TenDim)
		total=size(Ten_data)
!		T%Tensor_data=Ten_data
!		call assignment_com_dim1(T%Tensor_data,Ten_data)
		call storeTenData(T,Ten_data)
		T%rank=rank
		T%totalData=total
		T%TenDim=TenDim
		T%flag=.true.
		return
	end function
!*********************  getRank	 **********************
	integer function getRank(T)
		type(Tensor),intent(in) :: T
		getRank=T%rank
	end function
!*********************  getTotalData	 **********************
	integer function getTotalData(T)
		type(Tensor),intent(in) :: T
		getTotalData=T%totalData
	end function	
!*********************  getFlag	 **********************
	logical function getFlag(T)
		type(Tensor),intent(in) :: T
		getFlag=T%flag
	end function		
!*********************  generate *********************
!		generate a Tensor with random number
	type(Tensor) function generate(Tdim) result (T)
		integer,intent(in) :: Tdim(:)
		integer :: rank,totalData,i
		complex*16,allocatable :: Tensor_data(:)
		real*8 ::temp_real,temp_imag
		type(Dimension):: TenDim
		call init_random_seed()
		rank=size(Tdim)
		totalData=product(Tdim)
		allocate(Tensor_data(totalData))
		do i=1,totalData
			call random_number(temp_real)
			call random_number(temp_imag)
			Tensor_data(i)=dcmplx((0.5-temp_real),(0.5-temp_imag))
		end do
		TenDim=Tdim
		call storeTenData(T,Tensor_data)
		T%rank=rank
		T%totalData=totalData
		T%TenDim=TenDim
		T%flag=.true.
		return
	end function
!*********************  get Tensor dimenison *********************
	integer function getTenDim_i(T,inde)
		type(Tensor),intent(in) :: T
		integer,intent(in) :: inde
		getTenDim_i=Dim_i(T%TenDim,inde)
		return
	end function
	type(Dimension) function getTenDim(T)
		type(Tensor),intent(in) :: T
		getTenDim=T%TenDim
		return
	end function
	type(Dimension) function getTenSubDim(T,inde)
		type(Tensor),intent(in)::T
		integer,intent(in)::inde
		getTenSubDim=T%TenDim.sub.inde
		return
	end function
!********************* print_Tensor *********************       
	subroutine Tprint1(T)
		type(Tensor),intent(in) :: T
		if(T%flag) then
			write(*,*) "The rank of the Tensor is"
			write(*,*) T%rank
			write(*,*) "The number of  data of the Tensor is"
			write(*,*) T%totalData
			write(*,*) "The data of the Tensor is"
			write(*,*) T%Tensor_data
			call Dprint0(T%TenDim)
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
!********************* print_dim *********************            
	subroutine TDprint(T)
		type(Tensor) :: T
		if(T%flag) then
			write(*,*) "*** START ***"
			write(*,*) "The rank of the Tensor is"
			write(*,*) T%rank
			write(*,*) "The number of  data of the Tensor is"
			write(*,*) T%totalData
			call Dprint0(T%TenDim)
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
	subroutine Tprint2(T,printtype)
		type(Tensor),intent(in) :: T
		integer,intent(in)::printtype
		if(T%flag) then
			write(*,*) "The rank of the Tensor is"
			write(*,*) T%rank
			write(*,*) "The number of  data of the Tensor is"
			write(*,*) T%totalData
			write(*,*) "The data of the Tensor is"
			select case (printtype)
			case (1)
				write(*,*) real(T%Tensor_data)
			case (2)
				write(*,*) aimag(T%Tensor_data)
			case (3)
				write(*,'(99999F12.8)') real(T%Tensor_data)
			case (4)
				write(*,'(99999F12.8)') aimag(T%Tensor_data)
			case (5)
				write(*,'(99999F10.2)') real(T%Tensor_data)
			case (6)
				write(*,'(99999F10.2)') aimag(T%Tensor_data)
			end 	select
			call Dprint0(T%TenDim)
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
!********************* print_Matrix *********************
	subroutine TMprint1(T)
		type(Tensor) T
		complex*16,allocatable :: Tdata(:,:)
		complex*16,allocatable :: Tdata3(:,:,:)
		complex*16,allocatable :: Tdata4(:,:,:,:)
		integer,allocatable ::Tdim(:)
		integer :: i,j,k
		if(T%flag) then!if 1
			write(*,*) "*** START ***"
			select case(T%rank)
				case(1)
					write(*,*) T%Tensor_Data
					write(*,*) "*** END ***"
				case(2)
					allocate(Tdim(2))
					Tdim=T%TenDim
					allocate(Tdata(Tdim(1),Tdim(2)))
					Tdata=reshape(T%Tensor_Data,shape(Tdata))
						do i=1,Tdim(1)
							write (*,*) Tdata(i,:)
							write(*,*) ""
						end do 
						write(*,*) "*** END ***"
				case(3)
					allocate(Tdim(3))
					Tdim=T%TenDim
					allocate(Tdata3(Tdim(1),Tdim(2),Tdim(3)))
					Tdata3=reshape(T%Tensor_Data,shape(Tdata3))
					do j=1,Tdim(3)
						write(*,*) "(*,*,",j,")"
						do i=1,Tdim(1)
							write (*,*) Tdata3(i,:,j)
							write(*,*) ""
						end do
					end do 
					write(*,*) "*** END ***"
				case(4)
					allocate(Tdim(4))
					Tdim=T%TenDim
					allocate(Tdata4(Tdim(1),Tdim(2),Tdim(3),Tdim(4)))
					Tdata4=reshape(T%Tensor_Data,shape(Tdata4))
					do j=1,Tdim(4)
						do k=1,Tdim(3)
							write(*,*) "(*,*,",k,j,")"
							do i=1,Tdim(1)
							write (*,*) Tdata4(i,:,k,j)
							write(*,*) ""
							end do
						end do 
					end do
					write(*,*) "*** END ***"
				case default
					write(*,*) "rank of the Tensor is large than 4"
					write(*,*) "The data of the Tensor is"
					write(*,*) T%Tensor_data
					write(*,*) "***end***"
					write(*,*) "The dimension of the Tensor is"
					call Dprint0(T%TenDim)
					write(*,*) "The rank,total data are"
					write(*,*) T%rank,T%totalData
			end select
		else!if 1
			write(*,*) "There is no data"
		end if!if 1
		return
	end subroutine
	subroutine TMprint2(T,printtype)
		type(Tensor),intent(in) :: T
		integer,intent(in)::printtype
		complex*16,allocatable :: Tdata(:,:)
		complex*16,allocatable :: Tdata3(:,:,:)
		complex*16,allocatable :: Tdata4(:,:,:,:)
		integer,allocatable ::Tdim(:)
		integer :: i,j,k
		if(T%flag) then!if 1
			write(*,*) "*** START ***"
			select case(T%rank)
				case(1)
					select case (printtype)
					case (1)
						write(*,*) real(T%Tensor_data)
					case (2)
						write(*,*) aimag(T%Tensor_data)
					case (3)
						write(*,'(99999F12.8)') real(T%Tensor_data)
					case (4)
						write(*,'(99999F10.2)') real(T%Tensor_data)
					end 	select
					write(*,*) "*** END ***"
				case(2)
					allocate(Tdim(2))
					Tdim=T%TenDim
					allocate(Tdata(Tdim(1),Tdim(2)))
					Tdata=reshape(T%Tensor_Data,shape(Tdata))
					select case (printtype)
					case (1)
						do i=1,Tdim(1)
							write (*,*) real(Tdata(i,:))
							write(*,*) ""
						end do 
					case (2)
						do i=1,Tdim(1)
							write (*,*) aimag(Tdata(i,:))
							write(*,*) ""
						end do 
					case (3)
						do i=1,Tdim(1)
							write (*,'(99999F12.8)') real(Tdata(i,:))
							write(*,*) ""
						end do 
					case (4)
						do i=1,Tdim(1)
							write (*,'(99999F12.8)') aimag(Tdata(i,:))
							write(*,*) ""
						end do 
					case (5)
						do i=1,Tdim(1)
							write (*,'(99999F10.2)') real(Tdata(i,:))
							write(*,*) ""
						end do 
					case (6)
						do i=1,Tdim(1)
							write (*,'(99999F10.2)') aimag(Tdata(i,:))
							write(*,*) ""
						end do 
					end 	select
						write(*,*) "*** END ***"
				case(3)
					allocate(Tdim(3))
					Tdim=T%TenDim
					allocate(Tdata3(Tdim(1),Tdim(2),Tdim(3)))
					Tdata3=reshape(T%Tensor_Data,shape(Tdata3))
					select case (printtype)
					case (1)
						do j=1,Tdim(3)
							write(*,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write (*,*) real(Tdata3(i,:,j))
								write(*,*) ""
							end do
						end do 
					case (2)
						do j=1,Tdim(3)
							write(*,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write (*,*) aimag(Tdata3(i,:,j))
								write(*,*) ""
							end do
						end do 
					case (3)
						do j=1,Tdim(3)
							write(*,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write (*,'(99999F12.8)') real(Tdata3(i,:,j))
								write(*,*) ""
							end do
						end do 
					case (4)
						do j=1,Tdim(3)
							write(*,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write (*,'(99999F12.8)') aimag(Tdata3(i,:,j))
								write(*,*) ""
							end do
						end do 
					case (5)
						do j=1,Tdim(3)
							write(*,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write (*,'(99999F10.2)') real(Tdata3(i,:,j))
								write(*,*) ""
							end do
						end do 
					case (6)
						do j=1,Tdim(3)
							write(*,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write (*,'(99999F10.2)') aimag(Tdata3(i,:,j))
								write(*,*) ""
							end do
						end do 
					end 	select
					write(*,*) "*** END ***"
				case(4)
					allocate(Tdim(4))
					Tdim=T%TenDim
					allocate(Tdata4(Tdim(1),Tdim(2),Tdim(3),Tdim(4)))
					Tdata4=reshape(T%Tensor_Data,shape(Tdata4))
					select case (printtype)
					case (1)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(*,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write (*,*) real(Tdata4(i,:,k,j))
								write(*,*) ""
								end do
							end do 
						end do
					case (2)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(*,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write (*,*) aimag(Tdata4(i,:,k,j))
								write(*,*) ""
								end do
							end do 
						end do
					case (3)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(*,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write (*,'(99999F12.8)') real(Tdata4(i,:,k,j))
								write(*,*) ""
								end do
							end do 
						end do
					case (4)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(*,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write (*,'(99999F12.8)') aimag(Tdata4(i,:,k,j))
								write(*,*) ""
								end do
							end do 
						end do
					case (5)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(*,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write (*,'(99999F10.2)') real(Tdata4(i,:,k,j))
								write(*,*) ""
								end do
							end do 
						end do
					case (6)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(*,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write (*,'(99999F10.2)') aimag(Tdata4(i,:,k,j))
								write(*,*) ""
								end do
							end do 
						end do
					end 	select
					write(*,*) "*** END ***"
				case default
					call Tprint2(T,printtype)
			end select
		else!if 1
			write(*,*) "There is no data"
		end if!if 1
		return
	end subroutine


!cccccccccccccccc add          cccccccccccccccccc
	type(Tensor) function add(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		type(Dimension)::dim1,dim2
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*)"The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		add=set(dim1,T1%Tensor_data+T2%Tensor_data)
		return
	end function
	type(Tensor) function add_DTen1(T1,T2)
		type(Tensor),intent(in) :: T1
		type(DTensor),intent(in) ::T2
		type(Dimension)::dim1,dim2
		real*8,allocatable::realData(:)
		dim1=T1%TenDim
		dim2=.subDim.T2
		if(.not.(dim1.equ.dim2)) then
			write(*,*)"The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		realData=T2
		add_DTen1=set(dim1,T1%Tensor_data+realData)
		return
	end function
	type(Tensor) function add_DTen2(T1,T2)
		type(DTensor),intent(in) :: T1
		type(Tensor),intent(in) ::T2
		type(Dimension)::dim1,dim2
		real*8,allocatable::realData(:)
		dim1=.subDim.T1
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*)"The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		realData=T1
		add_DTen2=set(dim1,realData+T2%Tensor_data)
		return
	end function
	type(Tensor) function add_num(T1,num)
		type(Tensor),intent(in) :: T1
		complex*16,intent(in)::num
		type(Tensor)::T2
		integer,allocatable::dimT(:)
		dimT=T1%TenDim
		T2=eye(dimT)
		add_num=set(T1%TenDim,T1%Tensor_data+T2%Tensor_data)
		return
	end function
type(Tensor) function add_num_(num,T1)
		type(Tensor),intent(in) :: T1
		complex*16,intent(in)::num
		type(Tensor)::T2
		integer,allocatable::dimT(:)
		dimT=T1%TenDim
		T2=eye(dimT)
		add_num_=set(T1%TenDim,T1%Tensor_data+T2%Tensor_data)
		return
	end function
!cccccccccccccccc    minus       cccccccccccccccccc		
	type(Tensor) function minus(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		type(Dimension)::dim1,dim2
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*)"The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		minus=set(dim1,T1%Tensor_data-T2%Tensor_data)
		return
	end function
	type(Tensor) function minus_num(T1,num)
		type(Tensor),intent(in) :: T1
		complex*16,intent(in)::num
		type(Tensor)::T2
		integer,allocatable::dimT(:)
		dimT=T1%TenDim
		T2=eye(dimT)
		minus_num=set(T1%TenDim,T1%Tensor_data-T2%Tensor_data)
		return
	end function
	type(Tensor) function minus_num_(num,T1)
		type(Tensor),intent(in) :: T1
		complex*16,intent(in)::num
		type(Tensor)::T2
		integer,allocatable::dimT(:)
		dimT=T1%TenDim
		T2=eye(dimT)
		minus_num_=set(T1%TenDim,T1%Tensor_data-T2%Tensor_data)
		return
	end function
!cccccccccccccccc divide          cccccccccccccccccc	
	type(Tensor) function divide(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		type(Dimension)::dim1,dim2
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*) "The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		divide=set(dim1,T1%Tensor_data/T2%Tensor_data)
		return
	end function
!cccccccccccccccc divide_number          cccccccccccccccccc	
	type(Tensor) function divide_number(T1,num)
		type(Tensor),intent(in) :: T1
		complex*16,intent(in) :: num
		divide_number%rank=T1%rank
		divide_number%TenDim=T1%TenDim
		call storeTenData(divide_number,T1%Tensor_data/num)
		divide_number%totalData=T1%totalData
		divide_number%flag=.true.
		return
	end function
!ccccccccccccccccc divide_number          cccccccccccccccccc	
	type(Tensor) function divide_numberReal(T1,num)
		type(Tensor),intent(in) :: T1
		real*8,intent(in) :: num
		divide_numberReal%rank=T1%rank
		divide_numberReal%TenDim=T1%TenDim
		call storeTenData(divide_numberReal,T1%Tensor_data/num)
		divide_numberReal%totalData=T1%totalData
		divide_numberReal%flag=.true.
		return
	end function
!ccccccccccccccccc divide_number          cccccccccccccccccc	
	type(Tensor) function divide_numberReal4(T1,num)
		type(Tensor),intent(in) :: T1
		real*4,intent(in) :: num
		divide_numberReal4%rank=T1%rank
		divide_numberReal4%TenDim=T1%TenDim
		call storeTenData(divide_numberReal4,T1%Tensor_data/num)
		divide_numberReal4%totalData=T1%totalData
		divide_numberReal4%flag=.true.
		return
	end function			
!cccccccccccccccc multiply_Tensor          cccccccccccccccccc	
	type(Tensor) function multiply_Tensor(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		type(Dimension)::dim1,dim2
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*) "The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		multiply_Tensor=set(dim1,T1%Tensor_data*T2%Tensor_data)
		return
	end function
!cccccccccccccccc multiply_number          cccccccccccccccccc	
	type(Tensor) function multiply_number(T1,num)
		type(Tensor),intent(in) :: T1
		complex*16,intent(in) ::   num
		multiply_number%rank=T1%rank
		multiply_number%TenDim=T1%TenDim
		call storeTenData(multiply_number,T1%Tensor_data*num)
		multiply_number%totalData=T1%totalData
		multiply_number%flag=.true.
		return
	end function
		type(Tensor) function multiply_number_(num,T1)
		type(Tensor),intent(in) :: T1
		complex*16,intent(in) ::   num
		multiply_number_%rank=T1%rank
		multiply_number_%TenDim=T1%TenDim
		call storeTenData(multiply_number_,T1%Tensor_data*num)
		multiply_number_%totalData=T1%totalData
		multiply_number_%flag=.true.
		return
	end function
!cccccccccccccccc multiply_real          cccccccccccccccccc	
	type(Tensor) function multiply_real(T1,num)
		type(Tensor),intent(in) :: T1
		real*8,intent(in) ::   num
		multiply_real%rank=T1%rank
		multiply_real%TenDim=T1%TenDim
		call storeTenData(multiply_real,T1%Tensor_data*num)
		multiply_real%totalData=T1%totalData
		multiply_real%flag=.true.
		return
	end function
	type(Tensor) function multiply_real_(num,T1)
		type(Tensor),intent(in) :: T1
		real*8,intent(in) ::   num
		multiply_real_%rank=T1%rank
		multiply_real_%TenDim=T1%TenDim
		call storeTenData(multiply_real_,T1%Tensor_data*num)
		multiply_real_%totalData=T1%totalData
		multiply_real_%flag=.true.
		return
	end function	
!ccccccccccccccccc multiply_real4          cccccccccccccccccc	
	type(Tensor) function multiply_real4(T1,num)
		type(Tensor),intent(in) :: T1
		real*4,intent(in) ::   num
		multiply_real4%rank=T1%rank
		multiply_real4%TenDim=T1%TenDim
		call storeTenData(multiply_real4,T1%Tensor_data*num)
		multiply_real4%totalData=T1%totalData
		multiply_real4%flag=.true.
		return
	end function			
!ccccccccccccccccc permute_rank3   cccccccccccccccccc		
	type(Tensor) function permute_rank3(T1,index_not_permute)
		type(Tensor),intent(in) :: T1
		integer,intent(in) ::   index_not_permute
		integer ::i,newdimen(3),totaldata
		integer,allocatable :: dimen(:)
		type(Dimension) ::newTDim
		complex*16,allocatable :: Tdata(:,:,:),newdata(:,:,:)
		dimen=T1%TenDim
		allocate(Tdata(dimen(1),dimen(2),dimen(3)))
		Tdata=T1
		if(index_not_permute.eq.1) then
			newdimen(1)=dimen(1)
			newdimen(2)=dimen(3)
			newdimen(3)=dimen(2)
			newTDim=Dimpermute(T1%TenDim,(/1,3,2/))
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(1)
				newdata(i,:,:)=transpose(Tdata(i,:,:))
			end do
			call storeTenData(permute_rank3,newdata)
			permute_rank3%rank=3
			permute_rank3%totalData=T1%totalData
			permute_rank3%TenDim=newTDim
			permute_rank3%flag=.true.
		end if
		if(index_not_permute.eq.2) then
			newdimen(1)=dimen(3)
			newdimen(2)=dimen(2)
			newdimen(3)=dimen(1)
			newTDim=Dimpermute(T1%TenDim,(/3,2,1/))
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(2)
				newdata(:,i,:)=transpose(Tdata(:,i,:))
			end do
			call storeTenData(permute_rank3,newdata)
			permute_rank3%rank=3
			permute_rank3%totalData=T1%totalData
			permute_rank3%TenDim=newTDim
			permute_rank3%flag=.true.
		end if
		if(index_not_permute.eq.3) then
			newdimen(1)=dimen(2)
			newdimen(2)=dimen(1)
			newdimen(3)=dimen(3)
			newTDim=Dimpermute(T1%TenDim,(/2,1,3/))
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(3)
				newdata(:,:,i)=transpose(Tdata(:,:,i))
			end do
			call storeTenData(permute_rank3,newdata)
			permute_rank3%rank=3
			permute_rank3%totalData=T1%totalData
			permute_rank3%TenDim=newTDim
			permute_rank3%flag=.true.
		end if
	 end function

!cccccccccccccccc permute_rank2   cccccccccccccccccc			 
	type(Tensor) function permute_rank2 (T)
		type(Tensor),intent(in) :: T
		complex*16,allocatable :: Tdata(:,:),newdata(:,:)
		integer,allocatable ::dimen(:)
		type(Dimension) ::newTDim
		dimen=T%TenDim
		newTDim=Dimpermute(T%TenDim,(/2,1/))
		allocate(Tdata(dimen(1),dimen(2)))
		Tdata=T
		allocate(newdata(dimen(2),dimen(1)))
		newdata=transpose(Tdata)
		call storeTenData(permute_rank2,newdata)
		permute_rank2%rank=2
		permute_rank2%totalData=T%totalData
		permute_rank2%TenDim=newTDim
		permute_rank2%flag=.true.
	end function
!cccccccccccccccc permutation   cccccccccccccccccc
	type(Tensor) function permutation(T,newOrder)
		type(Tensor),intent(in) :: T
		integer,intent(in)::newOrder(:)
		integer,allocatable ::inde(:),newinde(:),maxdim(:),newmaxdim(:)
		integer::i,total,j,lendim,newaddress
		logical::fastoutput!in no permutation,output
		maxdim=T%TenDim
		lendim=size(maxdim)
		fastoutput=.true.
		i=1
		do while (fastoutput)
			if(newOrder(i).ne.i) then
				fastoutput=.false.
			end if
			i=i+1
			if(i.gt.lendim) then
				fastoutput=.false.
			end if
		end do
		if(fastoutput) then
			permutation=T
			return
		end if
		total=T%TotalData
		permutation%TenDim=Dimpermute(T%TenDim,newOrder)
		allocate(permutation%Tensor_Data(total))
		maxdim=T%TenDim
		newmaxdim=permutation%TenDim
		allocate(newinde(lendim))
		i=1
		do i=1,total
			call IndesToaddress(maxdim,inde,i)
			do j=1,lendim
				newinde(j)=inde(newOrder(j))
			end do
			newaddress=addressToIndes2(newmaxdim,newinde)
			permutation%Tensor_Data(newaddress)=T%Tensor_Data(i)
		end do
		permutation%rank=T%rank
		permutation%totalData=total
		permutation%flag=T%flag
	end function

!ccccccccccccccccc   productTen   ccccccccccccccccccc
!			the dimension of T1 and T2 should be 2 or 1
!			choose the righ form of
!						TensorProduct
!					TensorProduct1Dim
!					TensorProduct1Dim2
!		this is the old version.ProductTensor can do the same with
!	  higher perfromence.
	type(Tensor)	function productTen(T1,T2)
		type(Tensor),intent(in) :: T1
		type(Tensor),intent(in) :: T2
		integer :: rank1,rank2
		character*50 ::wt
		rank1=T1%rank
		rank2=T2%rank
		if((rank1.le.2).and.(rank2.le.2)) then
			if((rank1.eq.2).and.(rank2.eq.2)) then
				productTen=TensorProduct(T1,T2)
			else if((rank1.eq.1).and.(rank2.eq.1)) then
				productTen=TensorProduct1Dim1(T1,T2)
			else if((rank1.eq.1).or.(rank2.eq.1)) then
				productTen=TensorProduct1Dim2 (T1,T2) 
			else
				write(*,*) "error in the dimension of the Tensor"
			end if
		else
		end if
		return
	end function
				 	
!cccccccccccccccc TensorProduct   cccccccccccccccccc			 
!		both T1 and T2 should be a 2 \times 2 Tensor
	type(Tensor) function TensorProduct(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		integer :: totalData
		integer ::T1m,T1n,T2m,T2n
		type(Dimension) :: TenDim,dim1,dim2,newdim
		if((T1%rank.eq.2).and.(T2%rank.eq.2)) then
			T1m=T1.dim.1
			T1n=T1.dim.2
			T2m=T2.dim.1
			T2n=T2.dim.2
			totalData=T1m*T2n
			allocate(TensorProduct%Tensor_data(totalData))
			call ZGEMM('N','N',T1m,T2n,T1n,dcmplx(1d0,0d0),T1%Tensor_Data,T1m,&
				T2%Tensor_Data,T2m,dcmplx(0d0,0d0),TensorProduct%Tensor_data,T1m)
			dim1=T1%TenDim.sub.1
			dim2=T2%TenDim.sub.2
			newdim=dim1+dim2
			TensorProduct%rank=2
			TensorProduct%totalData=totalData
			TensorProduct%TenDim=newdim
			TensorProduct%flag=.true.
		else
			write(*,*)"The two Tensor input for product Should both be matrix"
		end if
	end function
!cccccccccccccccc  TensorProduct1Dim  cccccccccccccccccc
	complex*16 function TensorProduct1Dim(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		complex*16,allocatable :: T1data(:,:)
		complex*16,allocatable :: T2data(:,:)
		complex*16 :: Ndata
		integer ::T1m,T1n,T2m,T2n
		if((T1%rank.eq.1).and.(T2%rank.eq.1)) then
			T1n=T1.dim.1
			T2m=T1.dim.2
			if(T1n.ne.T2m) then
				write(*,*)"ERROR in TensorProduct1Dim1,stop"
				stop
			end if
			Ndata=ZDOTU(T1n,T1%Tensor_Data,1,T2%Tensor_Data,1)
			TensorProduct1Dim=Ndata
		else
			write(*,*)"The two Tensor input for product Should both be a vector"
		end if
	end function
!cccccccccccccccc  TensorProduct1Dim1  cccccccccccccccccc				 
	type(Tensor) function TensorProduct1Dim1(T1,T2)  
		type(Tensor),intent(in) :: T1,T2
		complex*16,allocatable :: T1data(:,:)
		complex*16,allocatable :: T2data(:,:)
		complex*16 :: Ndata(1,1)
		integer ::T1m,T1n,T2m,T2n
		if((T1%rank.eq.1).and.(T2%rank.eq.1)) then
			T1n=T1.dim.1
			T2m=T2.dim.1
			if(T1n.ne.T2m) then
				write(*,*)"ERROR in TensorProduct1Dim1,stop"
				stop
			end if
			Ndata=ZDOTU(T1n,T1%Tensor_Data,1,T2%Tensor_Data,1)
			TensorProduct1Dim1=set((/1/),(/Ndata/))
		else
			write(*,*)"The two Tensor input for product Should both be a vector"
		end if
		return
	end function
!cccccccccccccccc	 TensorProduct1Dim2   ccccccccccccccccc
	type(Tensor) function TensorProduct1Dim2 (T1,T2) result(ResultT)
		integer :: m,n,tempDim(2),RDim(1)
		type(Tensor),intent(in) :: T1,T2
		type(Tensor) :: temp
		type(Dimension) ::diml,dimr,newdim1,newdim2
		complex*16,allocatable::reData(:)
		if((T1%rank.eq.1).and.(T2%rank.eq.2)) then
			M=T2.dim.1
			N=T2.dim.2
			if((T1.dim.1) .ne. M) then
				write(*,*)"ERROR in TensorProduct1Dim2,stop"
				stop
			end if
			allocate(reData(N))
			call ZGEMV('T',M,N,dcmplx(1d0,0d0),T2%Tensor_Data,M,T1%Tensor_Data,1,dcmplx(0d0,0d0),reData,1)
			ResultT=set(T2%TenDim.sub.2,reData)
		else if((T1%rank.eq.2).and.(T2%rank.eq.1)) then
			M=T1.dim.1
			N=T1.dim.2
			if((T2.dim.1) .ne. N) then
				write(*,*)"ERROR in TensorProduct1Dim2,stop"
				stop
			end if
			allocate(reData(M))
			call ZGEMV('N',M,N,dcmplx(1d0,0d0),T1%Tensor_Data,M,T2%Tensor_Data,1,dcmplx(0d0,0d0),reData,1)
			ResultT=set(T1%TenDim.sub.1,reData)
		else
			write(*,*)"Error in the dimension of the input Tensors"
		end if
		return
	end function
!**************** ProductTensor  ***************************
!		ProductTensor regard the last index of T1 and the first index
!	  of T2 as the dimenion for matrix-product,other index will be see
!	  as another dimenison.T1 and T2 can be any rank,but the last dimenion
!	  of T1 and the first diemsion of T2 should be equal.
	type(Tensor) function ProductTensor (T1,T2)
		type(Tensor),intent(in) :: T1,T2
		integer::rank1,rank2,flag,T1m,T1n,T2m,T2n,T1l,T2l
		type(Dimension)::D1,D2,newD
		complex*16::Ndata
		integer::i,totaldata
		rank1=getRank(T1)
		rank2=getRank(T2)
		D1=.subDim.T1
		D2=.subDim.T2
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
				T1m=T1.dim.1
				T2n=T2.dim.1
				if(T1m.ne.T2n) then
					write(*,*)"ERROR in ProductTensor,case 1,stop"
					stop
				end if
				Ndata=ZDOTU(T1m,T1%Tensor_Data,1,T2%Tensor_Data,1)
				ProductTensor=set((/1/),(/Ndata/))
			case (2)
				if(rank2.ge.2) then
					newD=D2.sub.2
					do i=3,rank2
						newD=newD+(D2.sub.i)
					end do
					D2=DimConstract(D2,2,rank2)
				end if
				T2m=D2.i.1
				T2n=D2.i.2
				if((D1.i.1) .ne. T2m) then
					write(*,*)"ERROR in ProductTensor,case 2,stop"
					stop
				end if
				allocate(ProductTensor%Tensor_Data(T2n))
				call ZGEMV('T',T2m,T2n,dcmplx(1d0,0d0),T2%Tensor_Data,T2m,T1%Tensor_Data,&
							1,dcmplx(0d0,0d0),ProductTensor%Tensor_Data,1)
				ProductTensor%totalData=T2n
				ProductTensor%TenDim=newD
				ProductTensor%rank=DimSize(newD)
				ProductTensor%flag=.true.
		case (3)
				if(rank1.ge.2) then
					newD=D1.sub.1
					do i=2,rank1-1
						newD=newD+(D1.sub.i)
					end do
					D1=DimConstract(D1,1,rank1-2)
				end if
				T1m=D1.i.1
				T1n=D1.i.2
				if((D2.i.1) .ne. T1n) then
					write(*,*)"ERROR in ProductTensor,case 3,stop"
					stop
				end if
				allocate(ProductTensor%Tensor_Data(T1m))
				call ZGEMV('N',T1m,T1n,dcmplx(1d0,0d0),T1%Tensor_Data,T1m,T2%Tensor_Data,&
						1,dcmplx(0d0,0d0),ProductTensor%Tensor_Data,1)
				ProductTensor%totalData=T1m
				ProductTensor%TenDim=newD
				ProductTensor%rank=DimSize(newD)
				ProductTensor%flag=.true.
		case (4)
				if(rank1.ge.2) then
					newD=D1.sub.1
					do i=2,rank1-1
						newD=newD+(D1.sub.i)
					end do
					D1=DimConstract(D1,1,rank1-2)
				end if
				if(rank2.ge.2) then
					do i=2,rank2
						newD=newD+(D2.sub.i)
					end do
					D2=DimConstract(D2,2,rank2)
				end if
				if((DimSize(D1).ne.2).or.(DimSize(D2).ne.2)) then
					write(*,*)"error ProductTensor,dimension"
					stop
				end if
				T1m=D1.i.1
				T1n=D1.i.2
				T2m=D2.i.1
				T2n=D2.i.2
				totaldata=T1m*T2n
				allocate(ProductTensor%Tensor_data(totalData))
				call ZGEMM('N','N',T1m,T2n,T1n,dcmplx(1d0,0d0),T1%Tensor_Data,T1m,&
					T2%Tensor_Data,T2m,dcmplx(0d0,0d0),ProductTensor%Tensor_data,T1m)
				ProductTensor%totalData=totaldata
				ProductTensor%TenDim=newD
				ProductTensor%rank=DimSize(newD)
				ProductTensor%flag=.true.
		case default 
			write(*,*) "ERROR in ProductTensor,no such data"
			stop
		end 	select
		return
	end function
!************************  operateTensor  ****************************************
!	Operar is a matrix,the first element of every row specify
!		0: do nothing
!		1:contracts
!		2:decompose
!		3:decomposeAll
!	other element of every row is input parameter of the function
!		at last A*B .when flag=1,result store in A with B no change when productflag=0,
!	 the result will store in B
!		this subroutine run for big inoutA and B,do as less as possible operation on the
!	 Tensor_Data of inoutA and B.
	subroutine operateTensor2(A,B,oA,oB,productflag)
		type(Tensor),intent(inout)::A
		type(Tensor),intent(inout)::B
		integer,intent(in)::oA(:,:)
		integer,intent(in)::oB(:,:)
		integer,intent(in)::productflag
		type(Dimension)::Dim_backups!store the dimension which is no change
		integer::rank_backups
		if(productflag.eq.1) then
			Dim_backups=B%TenDim
			rank_backups=B%rank
		else if(productflag.eq.2) then
			Dim_backups=A%TenDim
			rank_backups=A%rank
		else
			write(*,*)"ERROR in operateTensor,productflag",productflag
			stop
		end if
		if(oA(1,1).ne.0) then
			call compose_decompse_subroutine2(A%TenDim,A%rank,oA)
		end if
		if(oB(1,1).ne.0) then
			call compose_decompse_subroutine2(B%TenDim,B%rank,oB)
		end if
		if(productflag.eq.1) then
			A=A*B
			B%TenDim=Dim_backups
			B%rank=rank_backups
		else if(productflag.eq.2) then
			B=A*B
			A%TenDim=Dim_backups
			A%rank=rank_backups
		else
			write(*,*)"ERROR in operateTensor"
		end if
		return
	end subroutine
	subroutine operateTensor1(A,B,oA,oB,productflag)
		type(Tensor),intent(inout)::A
		type(Tensor),intent(inout)::B
		integer,intent(in)::oA(:)
		integer,intent(in)::oB(:)
		integer,intent(in)::productflag
		type(Dimension)::Dim_backups!store the dimension which is no change
		integer::rank_backups
		if(productflag.eq.1) then
			Dim_backups=B%TenDim
			rank_backups=B%rank
		else if(productflag.eq.2) then
			Dim_backups=A%TenDim
			rank_backups=A%rank
		else
			write(*,*)"ERROR in operateTensor,productflag",productflag
			stop
		end if
		if(oA(1).ne.0) then
			call compose_decompse_subroutine1(A%TenDim,A%rank,oA)
		end if
		if(oB(1).ne.0) then
			call compose_decompse_subroutine1(B%TenDim,B%rank,oB)
		end if
		if(productflag.eq.1) then
			A=A*B
			B%TenDim=Dim_backups
			B%rank=rank_backups
		else if(productflag.eq.2) then
			B=A*B
			A%TenDim=Dim_backups
			A%rank=rank_backups
		else
			write(*,*)"ERROR in operateTensor"
		end if
		return
	end subroutine
		
		
		
		
		
		
		
		
		
!**************** contract   ****************
!		combine two index of the Tensor,which is con_index and con_index+1
	type(Tensor) function contract(T1,con_index)
		integer,intent(in) :: con_index
		type(Tensor),intent(in) :: T1
		type(Dimension) ::newdim
		newdim=DimConstract(T1%TenDim,con_index,1)
		call storeTenData(contract,T1%Tensor_Data)
		contract%rank=DimSize(newDim)
		contract%totalData=T1%totalData
		contract%TenDim=newdim
		contract%flag=.true.
		return
	end function
!**************** contracts   ****************
!		combine two index of the Tensor,which is index1 index1+1,..,index2-1,index2
!		if index2 larger than rank,the function will contract the index of index1 -> rank
	type(Tensor) function contracts(T1,index1,index2)
		integer,intent(in) :: index1,index2
		type(Tensor),intent(in) :: T1
		type(Dimension) ::newdim
		integer ::num
		num=index2-index1
		newdim=DimConstract(T1%TenDim,index1,num)
		call storeTenData(contracts,T1%Tensor_Data)
		contracts%rank=DimSize(newDim)
		contracts%totalData=T1%totalData
		contracts%TenDim=newdim
		contracts%flag=.true.
		return
	end function			
	type(Tensor) function contractsv(T1,vector)
		integer,intent(in) ::vector(2)
		type(Tensor),intent(in) :: T1
		type(Dimension) ::newdim
		integer ::num
		num=vector(2)-vector(1)
		newdim=DimConstract(T1%TenDim,vector(1),num)
		call storeTenData(contractsv,T1%Tensor_Data)
		contractsv%rank=DimSize(newDim)
		contractsv%totalData=T1%totalData
		contractsv%TenDim=newdim
		contractsv%flag=.true.
		return
	end function		
!*****************  decompose  *****************
! decompose the de_index index of the Tensor into n(1),n(2)
!		for example the de_index index is (1,2,3,4,..inde,inde+1,...rank)
!		(1,2,3,4,..inde,inde+1,...rank)-->(1,2,3,4,..inde),(inde+1,...rank)		
!		if inde larger than rank ,the function will return no change	
	type(Tensor) function decompose(T1,de_index,inde)
		type(Tensor),intent(in) :: T1
		integer,intent(in) :: de_index,inde
		type(Dimension) :: newDim
		newDim=DimDecompose(T1%TenDim,de_index,inde)
		call storeTenData(decompose,T1%Tensor_Data)
		decompose%rank=DimSize(newDim)
		decompose%totalData=T1%totalData
		decompose%TenDim=newdim
		decompose%flag=.true.
		return
	end function
	type(Tensor) function decomposev(T1,vector)
		type(Tensor),intent(in) :: T1
		integer,intent(in) ::vector(2)
		integer:: de_index,inde
		type(Dimension) :: newDim
		de_index=vector(1)
		inde=vector(2)
		newDim=DimDecompose(T1%TenDim,de_index,inde)
		call storeTenData(decomposev,T1%Tensor_Data)
		decomposev%rank=DimSize(newDim)
		decomposev%totalData=T1%totalData
		decomposev%TenDim=newdim
		decomposev%flag=.true.
		return
	end function
			
	type(Tensor) function decompose1(T1,de_index)
		type(Tensor),intent(in) :: T1
		integer,intent(in) :: de_index
		type(Dimension) :: newDim
		newDim=DimDecompose(T1%TenDim,de_index,1)
		call storeTenData(decompose1,T1%Tensor_Data)
		decompose1%rank=DimSize(newDim)
		decompose1%totalData=T1%totalData
		decompose1%TenDim=newdim
		decompose1%flag=.true.
		return
	end function
			
	type(Tensor) function decomposeAll(T1)
		type(Tensor),intent(in) :: T1
		integer:: i
		type(Dimension) :: newDim
		newDim=DimDecomposeAll(T1%TenDim)
		call storeTenData(decomposeAll,T1%Tensor_Data)
		decomposeAll%rank=DimSize(newDim)
		decomposeAll%totalData=T1%totalData
		decomposeAll%TenDim=newdim
		decomposeAll%flag=.true.
		return
	end function
!**********************  compose_decompse  *********************
!		compose_decompse do the contract or decompose on the Tensor
!		Operar is a matrix,the first element of every row specify
!			1:contracts
!			2:decompose
!			3:decomposeAll
!		other element of every row is input parameter of the function
	type(Tensor) function compose_decompse(T1,Operar)
		integer,intent(in)::Operar(:,:)
		type(Tensor),intent(in) :: T1
		integer::i,j,size1
		integer :: index1,index2
		integer ::num
		integer :: de_index,inde
		type(dimension)::dimen
		size1=size(Operar,1)
		dimen=T1%TenDim
		do i=1,size1
			if(Operar(i,1).eq.1) then
				index1=Operar(i,2)
				index2=Operar(i,3)
				num=index2-index1
				dimen=DimConstract(dimen,index1,num)
			else if(Operar(i,1).eq.2) then
				de_index=Operar(i,2)
				inde=Operar(i,3)
				dimen=DimDecompose(dimen,de_index,inde)
			else if(Operar(i,1).eq.3) then
				dimen=DimDecomposeAll(dimen)
			else
				write(*,*) "error in compose_decompse"
				stop
			end if
		end do
		call storeTenData(compose_decompse,T1%Tensor_Data)
		compose_decompse%TenDim=dimen
		compose_decompse%rank=DimSize(dimen)
		compose_decompse%flag=.true.
		compose_decompse%totalData=T1%totalData
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
	subroutine dimOperations(T1,Operar)
		integer,intent(in)::Operar(:,:)
		type(Tensor),intent(inout) :: T1
		integer::i,j,size1
		integer :: index1,index2
		integer ::num
		integer :: de_index,inde
		size1=size(Operar,1)
		do i=1,size1
			if(Operar(i,1).eq.1) then
				index1=Operar(i,2)
				index2=Operar(i,3)
				num=index2-index1
				T1%TenDim=DimConstract(T1%TenDim,index1,num)
				T1%rank=DimSize(T1%TenDim)
			else if(Operar(i,1).eq.2) then
				de_index=Operar(i,2)
				inde=Operar(i,3)
				T1%TenDim=DimDecompose(T1%TenDim,de_index,inde)
				T1%rank=DimSize(T1%TenDim)
			else if(Operar(i,1).eq.3) then
				T1%TenDim=DimDecomposeAll(T1%TenDim)
				T1%rank=DimSize(T1%TenDim)
			else
				write(*,*) "error in compose_decompse"
				stop
			end if
		end do
		return
	end subroutine
	subroutine dimOperation1(T1,Operat)
		integer,intent(in)::Operat(:)
		type(Tensor),intent(inout) :: T1
		integer,allocatable::Operat2(:,:)
		allocate(Operat2(1,size(Operat)))
		Operat2(1,:)=Operat
		call dimOperations(T1,Operat2)
		return
	end subroutine
	subroutine compose_decompse_subroutine2(TenDim,rank,Operar)
		integer,intent(in)::Operar(:,:)
		type(Dimension),intent(inout) :: TenDim
		integer,intent(inout)::rank
		integer::i,j,size1
		integer :: index1,index2
		integer ::num
		integer :: de_index,inde
		size1=size(Operar,1)
		do i=1,size1
			if(Operar(i,1).eq.1) then
				index1=Operar(i,2)
				index2=Operar(i,3)
				num=index2-index1
				TenDim=DimConstract(TenDim,index1,num)
				rank=DimSize(TenDim)
			else if(Operar(i,1).eq.2) then
				de_index=Operar(i,2)
				inde=Operar(i,3)
				TenDim=DimDecompose(TenDim,de_index,inde)
				rank=DimSize(TenDim)
			else if(Operar(i,1).eq.3) then
				TenDim=DimDecomposeAll(TenDim)
				rank=DimSize(TenDim)
			else
				write(*,*) "error in compose_decompse_subroutine"
				stop
			end if
		end do
		return
	end subroutine
	subroutine compose_decompse_subroutine1(TenDim,rank,Operar)
		integer,intent(in)::Operar(:)
		type(Dimension),intent(inout) :: TenDim
		integer,intent(inout)::rank
		integer :: index1,index2
		integer ::num
		integer :: de_index,inde
		if(Operar(1).eq.1) then
			index1=Operar(2)
			index2=Operar(3)
			num=index2-index1
			TenDim=DimConstract(TenDim,index1,num)
			rank=DimSize(TenDim)
		else if(Operar(1).eq.2) then
			de_index=Operar(2)
			inde=Operar(3)
			TenDim=DimDecompose(TenDim,de_index,inde)
			rank=DimSize(TenDim)
		else if(Operar(1).eq.3) then
			TenDim=DimDecomposeAll(TenDim)
			rank=DimSize(TenDim)
		else
			write(*,*) "error in compose_decompse_subroutine"
			stop
		end if
		return
	end subroutine

	
!***************** equal_of_Tensor *****************
	logical function equal_of_Tensor(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		integer :: l
		equal_of_Tensor=.true.
		if(T1%rank.ne.T2%rank) then
			equal_of_Tensor=.false.
			return
		end if
		if(.not.equal_of_dim(T1%TenDim,T2%TenDim)) then
			equal_of_Tensor=.false.
			return
		end if
		l=count(abs(T1%Tensor_data-T2%Tensor_data).gt.1d-10)
		if(l.gt.0) then
			equal_of_Tensor=.false.
			return
		end if
		return
	end function
!*****************  maxElement  *****************
	real*8 function maxElement(T)
		type(Tensor),intent(in) :: T
		maxElement=maxval(abs(T%Tensor_Data))
		return
	end function
!*****************  maxRealElement  *****************
	real*8 function maxRealElement(T)
		type(Tensor),intent(in) :: T
		maxRealElement=maxval(real(T%Tensor_Data))
		return
	end function
!*****************  minElement  *****************
	real*8 function minElement(T)
		type(Tensor),intent(in) :: T
		minElement=minval(abs(T%Tensor_Data))
		return
	end function
!*****************  minRealElement  *****************
	real*8 function minRealElement(T)
		type(Tensor),intent(in) :: T
		minRealElement=minval(real(T%Tensor_Data))
		return
	end function
!ccccccccccccccccc svd cccccccccccccccccc	

!*****************  svdright  *****************
! 			T should be two dimension		
!			T=U*s*V^T
!			return the (U*s) in T1
!			and V^T in T2
!			the data s(i) in the matrix s will be discard if
!			s(i)<s(1)*check_discard_svd.
!			if check_discard_svd<=0,all the data in s will be keep
!			note:I do not write the faster version of this part
	subroutine svdright(T,T1,T2,check_discard)
		type(Tensor),intent(in) :: T
		type(Tensor),intent(out) :: T1,T2
		real*8,intent(in) :: check_discard
		complex*16,allocatable :: T1data(:,:)
		complex*16,allocatable :: T2data(:,:)
		complex*16,allocatable :: newT2data(:,:)
		complex*16,allocatable :: newT1data(:,:)
		complex*16,allocatable :: Tdata(:,:)
		complex*16,allocatable :: work(:)
		integer m,n,info,lw,rwdim,i,j
		integer,allocatable :: iw(:)
		real*8,allocatable :: rw(:),sdata(:)
		integer :: ms_max
		type(Dimension) :: T1Dim,T2Dim
		real*8::cputime(2)
		if(T%rank.eq.2) then
			m=T.dim.1
			n=T.dim.2
			lw=min(M,N)*min(M,N)+2*min(M,N)+max(M,N)
			rwdim=min(M,N)*max(5*min(M,N)+7,2*max(M,N)+2*min(M,N)+1)
			allocate(rw(rwdim))
			allocate(iw(8*min(M,N)))
			allocate(T1data(m,m))
			allocate(T2data(m,n))
			allocate(Tdata(m,n))
			allocate(sdata(m)) 
			if(m.gt.n) then
				write(*,*)"error in svdright"
			end if
			allocate(work(lw))
			Tdata=reshape(T%Tensor_Data,(/m,n/))
			call ZGESDD('S',m,n,Tdata,m,sdata,T1data,m,T2data,m,WORK,lw,rw,iw,INFO)
			deallocate(Tdata)
			if(info.ne.0) then
				write (*,*) "Error in svd ,info=",info
			end if
			if(check_discard.gt.0) then 
				ms_max=1
				do i=2,m
					if(sdata(i).gt.sdata(1)*check_discard) then
						ms_max=ms_max+1
					end if
				end do
			else 
				ms_max=m
			end if
			allocate(newT2data(ms_max,n))
			do j=1,n
				do i=1,ms_max
					newT2data(i,j)=sdata(i)*T2data(i,j)
				end do
			end do
			if(ms_max.eq.m) then
				T1Dim=(/m/)
				T1Dim=(T%TenDim.sub.1)+T1Dim
				T2Dim=(/ms_max/)
				T2Dim=T2Dim+(T%TenDim.sub.2)
				T1=set(T1Dim,reshape(T1data,(/m*m/)))
				T2=set(T2Dim,reshape(newT2data,(/n*ms_max/)))
			else
				allocate(newT1data(m,ms_max))
				newT1data=T1data(:,1:ms_max)
				T1Dim=(/ms_max/)
				T1Dim=(T%TenDim.sub.1)+T1Dim
				T2Dim=(/ms_max/)
				T2Dim=T2Dim+(T%TenDim.sub.2)
				T1=set(T1Dim,reshape(newT1data,(/m*ms_max/)))
				T2=set(T2Dim,reshape(newT2data,(/n*ms_max/)))
			end if
		else
			write(*,*) "Input Tensor should be 2 dimension"
		end if
		return
	end subroutine												
!cccccccccccccccc svdleft cccccccccccccccccc	
! 				T should be two dimension		
!				T=U*s*V^T
!				return the transpose U in T1
!			and (s*V^T) in T2
!				the data s(i) in the matrix s will be discard if
!			s(i)<s(1)*check_discard_svd.
!				if check_discard_svd<0,all the data in s will be keep
!			wornning:I do not check the part check_discard>0
	subroutine svdleft(T,T1,T2,check_discard)
		type(Tensor),intent(in) :: T
		type(Tensor),intent(out) :: T1,T2
		real*8,intent(in) :: check_discard
		complex*16,allocatable :: T1data(:)
		complex*16,allocatable :: T2data(:)
		complex*16,allocatable :: newT2data(:)
		complex*16,allocatable :: Tdata(:)
		complex*16,allocatable :: work(:)
		integer m,n,info,lw,rwdim,i,j
		integer,allocatable :: iw(:)
		real*8,allocatable :: rw(:),sdata(:)
		integer :: ns_max,inde_ij,inde_ij2
		real*8::cputime(2)
		type(Dimension) :: T1Dim,T2Dim
		if(T%rank.eq.2) then
			m=T.dim.1
			n=T.dim.2
			lw=min(M,N)*min(M,N)+2*min(M,N)+max(M,N)
			rwdim=min(M,N)*max(5*min(M,N)+7,2*max(M,N)+2*min(M,N)+1)
			allocate(rw(rwdim))
			allocate(iw(8*min(M,N)))
			allocate(T1data(m*n))
			allocate(T2data(n*n))
			allocate(Tdata(m*n))
			allocate(sdata(n))
			if(m.lt.n) then
				write(*,*)"error in svdleft"
			end if
			allocate(work(lw))
			Tdata=T
			call ZGESDD('S',m,n,Tdata,m,sdata,T1data,m,T2data,n,WORK,lw,rw,iw,INFO)
			deallocate(Tdata)
			if(info.ne.0) then
				write (*,*) "Error in svd ,info=",info
			end if
			if(check_discard.gt.0) then 
				ns_max=1
				do i=2,n
					if(sdata(i).gt.sdata(1)*check_discard) then
						ns_max=ns_max+1
					end if
				end do
			else 
				ns_max=n
			end if
			do j=1,ns_max
				do i=1,m
					inde_ij=addressToIndes2((/m,n/),(/i,j/))
					T1data(inde_ij)=T1data(inde_ij)*sdata(j)
				end do
			end do
			if(ns_max.eq.n) then
				T1Dim=(/n/)
				T1Dim=(T%TenDim.sub.1)+T1Dim
				T2Dim=(/n/)
				T2Dim=T2Dim+(T%TenDim.sub.2)
				T2=set(T2Dim,T2data)
				T1=set(T1Dim,T1data)
			else
				allocate(newT2data(ns_max*n))
				do j=1,n
					do i=1,ns_max
						inde_ij=addressToIndes2((/n,n/),(/i,j/))
						inde_ij2=addressToIndes2((/ns_max,n/),(/i,j/))
						newT2data(inde_ij)=T2data(inde_ij2)
					end do
				end do
				T1Dim=(/ns_max/)
				T1Dim=(T%TenDim.sub.1)+T1Dim
				T2Dim=(/ns_max/)
				T2Dim=T2Dim+(T%TenDim.sub.2)
				T2=set(T2Dim,newT2data)
				T1=set(T1Dim,T1data)
			end if
		else
			write(*,*) "Input Tensor should be 2 dimension"
		end if
		return
	end subroutine
!Solves a general system of linear equations
!A*X=B,find X
!X is a vector or a matrix,the dimension of which is the same as B
	type(Tensor) function linequ(A,B)
		type(Tensor),intent(in)::A,B
		complex*16,allocatable::Xdata(:),Adata(:)
		type(dimension)::Xdim
		integer,allocatable::IPIV(:)
		integer::Na,Nb,INFO
		Na=A.dim.1
		if(Na.ne.(A.dim.2)) then
			write(*,*)"error in linequ,dimension of A"
			stop
		end if
		if(getRank(B).eq.1)then
			Nb=1
		else
			Nb=B.dim.2
		end if
		if(Na.ne.(B.dim.1)) then
			write(*,*)"error in linequ,dimension of A and B"
			stop
		end if
		allocate(IPIV(Na))
		Xdata=B
		Adata=A!the subroutine ZGESV will change A
		call ZGESV(Na,Nb,Adata,Na,IPIV,Xdata,Na,INFO)
		if(INFO.ne.0)then
			write(*,*)"ZGESV is not successful"
			write(*,*)"INFO",INFO
			stop
		end if
		
		Xdim=(A%TenDim.sub.2)
		if(getRank(B).ne.1)then
			Xdim=Xdim+(B%TenDim.sub.2)
		end if
		linequ=set(Xdim,Xdata)
		return
	end function
!Solves a general system of linear equations
!A*X=B,find X
!X is a vector or a matrix,the dimension of which is the same as B
!on output,A will change and X is in B	
	subroutine linequ_routine(A,B)
		type(Tensor),intent(inout)::A,B
		integer,allocatable::IPIV(:)
		integer::Na,Nb,INFO
		type(dimension)::Xdim
		Na=A.dim.1
		if(Na.ne.(A.dim.2)) then
			write(*,*)"error in linequ,dimension of A"
			stop
		end if
		if(getRank(B).eq.1)then
			Nb=1
		else
			Nb=B.dim.2
		end if
		if(Na.ne.(B.dim.1)) then
			write(*,*)"error in linequ,dimension of A and B"
			stop
		end if
		allocate(IPIV(Na))
		call ZGESV(Na,Nb,A%Tensor_Data,Na,IPIV,B%Tensor_Data,Na,INFO)
		if(INFO.ne.0)then
			write(*,*)"ZGESV is not successful"
			write(*,*)"INFO",INFO
			stop
		end if
		Xdim=(A%TenDim.sub.2)
		if(getRank(B).ne.1)then
			Xdim=Xdim+(B%TenDim.sub.2)
		end if
		call resetdim(B,Xdim)
		return
	end subroutine
!********************************************************
!		eigen value near some value, say s0
!		on output,a Tensor link,link:[eigenvalue]->[eigenvector]
!		link.i.2 is a matrix of eigenvector , and the columns of which is the approximate 
!	eigenvectors (Ritz vectors) corresponding eigenvalue
!	ncv number of lonzcos vector
!	if outtype=1 eigenvalues
!	if outtype=2 eigenvalues and eigenvectors
!	if outtype=3 eigenvalues,but using |T-s0|^2 to find the eigenvalue near s0
!	if outtype=4 eigenvalues and eigenvectors,but using |T-s0|^2 to find the eigenvalue near s0
	 subroutine eigen1(resultlink,T,s0,num_of_eig,ncv,tol,outtype)
	 	type(Tensorlink),intent(inout)::resultlink
		type(Tensor),intent(in)::T
		complex*16,intent(in)::s0
		integer,intent(in)::num_of_eig,ncv,outtype
		real*8,intent(in)::tol
		complex*16,allocatable::A(:),v(:,:),resid(:),d(:),Cdata(:),workdata(:)
		complex*16,allocatable::workd(:),workev(:),workl(:)
		logical::outvector
		integer::ido,info,ishfts,maxitr,mode,iparam(11), ipntr(14)
		integer,allocatable::ipiv(:)
		logical::continu
		real*8,allocatable::rwork(:)
		logical,allocatable::selework(:)
		integer,allocatable::IPIV2(:)
		integer::icouter,info2
		complex*16,allocatable::C0(:,:),C(:,:)
		integer::i,j,n,lworkl,ierr
		type(Tensornode),pointer::res_eigvalue,res_eigenvector
		type(Tensor)::T2
		if(outtype.eq.1)then
			outvector=.false.
		else if(outtype.gt.4) then
			write(*,*)"ERROR in eigen,on such out type"
			write(*,*)outtype
			stop
		else
			outvector=.true.
		end if
		if(getRank(T).ne.2) then
			write(*,*)"ERROR in dimension of input Tensor,eigen"
			write(*,*)"stop"
			stop
		end if
		n=T.dim.1
		if(n.ne.(T.dim.2)) then
		   write(*,*)"ERROR in dimension of input Tensor,eigen"
		   write(*,*)"Should be a matrix of n * n"
			write(*,*)"stop"
			stop
		end if
		call cleanlink(resultlink)
		if((outtype.eq.1).or.(outtype.eq.2))then
		   allocate(C0(n,n))
		   C0=0.
		   do i=1,n
				C0(i,i)=s0
		   end do
		   C=T
		   C0=C-C0
		else
			T2=T-(s0*eye(n,n))
			T2=(.h.T2)*T2
		end if
      
      lworkl = 3*ncv**2+5*ncv
		allocate(workd(3*n))
		allocate(workev(2*ncv))
		allocate(workl(lworkl))
		allocate(ipiv(n))
		allocate(rwork(ncv))
		allocate(v(n,ncv))
		allocate(resid(n))
		allocate(d(num_of_eig+1))
		allocate(IPIV2(n))
		allocate(workdata(n))
		allocate(selework(ncv))
		!tol    = 0.0 
		ido    = 0
		info   = 0
		info2=0
				
				
		ishfts = 1
		maxitr = 300
		if((outtype.eq.1).or.(outtype.eq.2))then
			mode   = 3
		else
			mode   = 1
		end if
		iparam(1) = ishfts 
		iparam(3) = maxitr 
		iparam(7) = mode 	 
		continu=.true.
		if((outtype.eq.1).or.(outtype.eq.2))then!if[1],znaupd to find eigenvalue near s0
			do while(continu)
				call znaupd(ido,'I',n,'LM',num_of_eig,tol,resid,ncv,v,n,iparam,ipntr,workd,workl,lworkl,rwork,info)
				if (ido .eq. -1 .or. ido .eq. 1 ) then
					call zcopy ( n, workd(ipntr(1)),1, workd(ipntr(2)), 1)
					C=C0
					call ZGESV(n,1,C,n,IPIV2,workd(ipntr(2)),n,info2)
				else
					continu=.false.
				end if
			end do
			if(ido.ne.99)then
				write(*,*)"znaupd do not success"
				write(*,*)info
				stop
			end if
			call zneupd(outvector,'A',selework,d,v,n,s0,workev,'I',n,'LM',num_of_eig,tol,resid,ncv,v,n,iparam,&
						ipntr,workd,workl,lworkl,rwork,ierr)
			if(outvector) then
				allocate(res_eigvalue)
				allocate(res_eigenvector)
				res_eigvalue%Ten=d(1:num_of_eig)
				res_eigenvector%Ten=v(:,1:num_of_eig)
				call push_back(resultlink,res_eigvalue)
				call push_back(resultlink,res_eigenvector)
			else
				allocate(res_eigvalue)
				res_eigvalue%Ten=d(1:num_of_eig)
				call push_back(resultlink,res_eigvalue)
			end if
		else!if[1],znaupd to find the min eigenvalue of |T-s0|^2
			do while(continu)
				call znaupd(ido,'I',n,'SM',num_of_eig,tol,resid,ncv,v,n,iparam,ipntr,workd,workl,lworkl,rwork,info)
				if (ido .eq. -1 .or. ido .eq. 1 ) then
					call ZGEMV('N',n,n,dcmplx(1d0,0d0),T2%Tensor_Data,n,workd(ipntr(1):ipntr(1)+n-1),&
						1,dcmplx(0d0,0d0),workd(ipntr(2):ipntr(2)+n-1),1)
				else
					continu=.false.
				end if
			end do
			if(ido.ne.99)then
				write(*,*)"znaupd do not success"
				write(*,*)info
				stop
			end if
			call zneupd(outvector,'A',selework,d,v,n,s0,workev,'I',n,'SM',num_of_eig,tol,resid,ncv,v,n,iparam,&
						ipntr,workd,workl,lworkl,rwork,ierr)
			allocate(res_eigvalue)
			allocate(res_eigenvector)
			res_eigenvector%Ten=v(:,1:num_of_eig)
			do i=1,num_of_eig
				T2=subTen(res_eigenvector%Ten,(/-2,i/))
				T2=T*T2
				d(i)=subTen(res_eigenvector%Ten,(/-2,i/)) .x. T2
			end do
			res_eigvalue%Ten=d(1:num_of_eig)
			call push_back(resultlink,res_eigvalue)
			if(outtype.eq.4)then
				call push_back(resultlink,res_eigenvector)
			end if
		end if!if[1]
		
		return
	end subroutine
	 subroutine eigen2(resultlink,T,s0,num_of_eig,outtype)
	 	type(Tensorlink),intent(inout)::resultlink
		type(Tensor),intent(in)::T
		complex*16,intent(in)::s0
		integer,intent(in)::num_of_eig,outtype
		complex*16,allocatable::A(:),v(:,:),resid(:),d(:),Cdata(:),workdata(:)
		complex*16,allocatable::workd(:),workev(:),workl(:)
		logical::outvector
		real*8::tol
		integer::ido,info,ishfts,maxitr,mode,iparam(11), ipntr(14)
		integer,allocatable::ipiv(:)
		logical::continu
		real*8,allocatable::rwork(:)
		logical,allocatable::selework(:)
		integer,allocatable::IPIV2(:)
		integer::icouter,info2,ncv
		complex*16,allocatable::C0(:,:),C(:,:)
		integer::i,j,n,lworkl,ierr
		type(Tensornode),pointer::res_eigvalue,res_eigenvector
		type(Tensor)::T2
		if(outtype.eq.1)then
			outvector=.false.
		else if(outtype.gt.4) then
			write(*,*)"ERROR in eigen,on such out type"
			write(*,*)outtype
			stop
		else
			outvector=.true.
		end if
		if(getRank(T).ne.2) then
			write(*,*)"ERROR in dimension of input Tensor,eigen"
			write(*,*)"stop"
			stop
		end if
		n=T.dim.1
		ncv=(2+num_of_eig+n)/2
		if(ncv.gt.max_ncv) then
			ncv=max_ncv
		end if
		if(n.ne.(T.dim.2)) then
		   write(*,*)"ERROR in dimension of input Tensor,eigen"
		   write(*,*)"Should be a matrix of n * n"
			write(*,*)"stop"
			stop
		end if
		call cleanlink(resultlink)
		if((outtype.eq.1).or.(outtype.eq.2))then
		   allocate(C0(n,n))
		   C0=0.
		   do i=1,n
				C0(i,i)=s0
		   end do
		   C=T
		   C0=C-C0
		else
			T2=T-(s0*eye(n,n))
			T2=(.h.T2)*T2
		end if
      
      lworkl = 3*ncv**2+5*ncv
		allocate(workd(3*n))
		allocate(workev(2*ncv))
		allocate(workl(lworkl))
		allocate(ipiv(n))
		allocate(rwork(ncv))
		allocate(v(n,ncv))
		allocate(resid(n))
		allocate(d(num_of_eig+1))
		allocate(IPIV2(n))
		allocate(workdata(n))
		allocate(selework(ncv))
		tol    = 0.0 
		ido    = 0
		info   = 0
		info2=0
				
				
		ishfts = 1
		maxitr = 300
		if((outtype.eq.1).or.(outtype.eq.2))then
			mode   = 3
		else
			mode   = 1
		end if
		iparam(1) = ishfts 
		iparam(3) = maxitr 
		iparam(7) = mode 	 
		continu=.true.
		if((outtype.eq.1).or.(outtype.eq.2))then!if[1],znaupd to find eigenvalue near s0
			do while(continu)
				call znaupd(ido,'I',n,'LM',num_of_eig,tol,resid,ncv,v,n,iparam,ipntr,workd,workl,lworkl,rwork,info)
				if (ido .eq. -1 .or. ido .eq. 1 ) then
					call zcopy ( n, workd(ipntr(1)),1, workd(ipntr(2)), 1)
					C=C0
					call ZGESV(n,1,C,n,IPIV2,workd(ipntr(2)),n,info2)
				else
					continu=.false.
				end if
			end do
			if(ido.ne.99)then
				write(*,*)"znaupd do not success"
				write(*,*)info
				stop
			end if
			call zneupd(outvector,'A',selework,d,v,n,s0,workev,'I',n,'LM',num_of_eig,tol,resid,ncv,v,n,iparam,&
						ipntr,workd,workl,lworkl,rwork,ierr)
			if(outvector) then
				allocate(res_eigvalue)
				allocate(res_eigenvector)
				res_eigvalue%Ten=d(1:num_of_eig)
				res_eigenvector%Ten=v(:,1:num_of_eig)
				call push_back(resultlink,res_eigvalue)
				call push_back(resultlink,res_eigenvector)
			else
				allocate(res_eigvalue)
				res_eigvalue%Ten=d(1:num_of_eig)
				call push_back(resultlink,res_eigvalue)
			end if
		else!if[1],znaupd to find the min eigenvalue of |T-s0|^2
			do while(continu)
				call znaupd(ido,'I',n,'SM',num_of_eig,tol,resid,ncv,v,n,iparam,ipntr,workd,workl,lworkl,rwork,info)
				if (ido .eq. -1 .or. ido .eq. 1 ) then
					call ZGEMV('N',n,n,dcmplx(1d0,0d0),T2%Tensor_Data,n,workd(ipntr(1):ipntr(1)+n-1),&
						1,dcmplx(0d0,0d0),workd(ipntr(2):ipntr(2)+n-1),1)
				else
					continu=.false.
				end if
			end do
			if(ido.ne.99)then
				write(*,*)"znaupd do not success"
				write(*,*)info
				stop
			end if
			call zneupd(outvector,'A',selework,d,v,n,s0,workev,'I',n,'SM',num_of_eig,tol,resid,ncv,v,n,iparam,&
						ipntr,workd,workl,lworkl,rwork,ierr)
			allocate(res_eigvalue)
			allocate(res_eigenvector)
			res_eigenvector%Ten=v(:,1:num_of_eig)
			do i=1,num_of_eig
				T2=subTen(res_eigenvector%Ten,(/-2,i/))
				T2=T*T2
				d(i)=subTen(res_eigenvector%Ten,(/-2,i/)) .x. T2
			end do
			res_eigvalue%Ten=d(1:num_of_eig)
			call push_back(resultlink,res_eigvalue)
			if(outtype.eq.4)then
				call push_back(resultlink,res_eigenvector)
			end if
		end if!if[1]
	end subroutine
	subroutine eigen3(resultlink,T,s0,num_of_eig)!outvector=.ture.
		type(Tensorlink),intent(inout)::resultlink
		type(Tensor),intent(in)::T
		complex*16,intent(in)::s0
		integer,intent(in)::num_of_eig
		call eigen2(resultlink,T,s0,num_of_eig,1)
		return
	end subroutine
!********************************************************
!		the max of min eigen value
!		on output,a Tensor link,link:[eigenvalue]->[eigenvector]
!		link.i.2 is a matrix of eigenvector , and the columns of which is the approximate 
!	eigenvectors (Ritz vectors) corresponding eigenvalue
!	ncv number of lonzcos vector	
!	if outvector=.false. no eigenvectors
!				WHICH   Character*2.  (INPUT)
!          'LM' -> want the NEV eigenvalues of largest magnitude.
!          'SM' -> want the NEV eigenvalues of smallest magnitude.
!          'LR' -> want the NEV eigenvalues of largest real part.
!          'SR' -> want the NEV eigenvalues of smallest real part.
!          'LI' -> want the NEV eigenvalues of largest imaginary part.
!          'SI' -> want the NEV eigenvalues of smallest imaginary part.
	subroutine eigen_max_min1(resultlink,T,WHICH,num_of_eig,ncv,tol,outvector)
		type(Tensorlink),intent(inout)::resultlink
		type(Tensor),intent(in)::T
		Character*2,intent(in)::WHICH
		logical,intent(in)::outvector
		integer,intent(in)::num_of_eig,ncv
		complex*16,allocatable::A(:),v(:,:),resid(:),d(:),Cdata(:),workdata(:)
		complex*16,allocatable::workd(:),workev(:),workl(:)
		complex*16::s0
		real*8,intent(in)::tol
		integer::ido,info,ishfts,maxitr,mode,iparam(11), ipntr(14)
		integer,allocatable::ipiv(:)
		logical::continu
		real*8,allocatable::rwork(:)
		logical,allocatable::selework(:)
		integer,allocatable::IPIV2(:)
		integer::icouter,info2
	!	complex*16,allocatable::C0(:,:),C(:,:)
		integer::i,j,n,lworkl,ierr
		type(Tensornode),pointer::res_eigvalue,res_eigenvector
		
		if(getRank(T).ne.2) then
			write(*,*)"ERROR in dimension of input Tensor,eigen"
			write(*,*)"stop"
			stop
		end if
		n=T.dim.1
		if(n.ne.(T.dim.2)) then
		   write(*,*)"ERROR in dimension of input Tensor,eigen"
		   write(*,*)"Should be a matrix of n * n"
			write(*,*)"stop"
			stop
		end if
		
     call cleanlink(resultlink)
      
      lworkl = 3*ncv**2+5*ncv
		allocate(workd(3*n))
		allocate(workev(2*ncv))
		allocate(workl(lworkl))
		allocate(ipiv(n))
		allocate(rwork(ncv))
		allocate(v(n,ncv))
		allocate(resid(n))
		allocate(d(num_of_eig+1))
		allocate(IPIV2(n))
		allocate(workdata(n))
		allocate(selework(ncv))
		!tol    = 0.0 
		ido    = 0
		info   = 0
		info2=0
				
				
		ishfts = 1
		maxitr = 300
		mode   = 1
		iparam(1) = ishfts 
		iparam(3) = maxitr 
		iparam(7) = mode 	 
		continu=.true.
		
		do while(continu)
			call znaupd(ido,'I',n,WHICH,num_of_eig,tol,resid,ncv,v,n,iparam,ipntr,workd,workl,lworkl,rwork,info)
			if (ido .eq. -1 .or. ido .eq. 1 ) then
				call ZGEMV('N',n,n,dcmplx(1d0,0d0),T%Tensor_Data,n,workd(ipntr(1):ipntr(1)+n-1),&
						1,dcmplx(0d0,0d0),workd(ipntr(2):ipntr(2)+n-1),1)
				!call zcopy ( n, workd(ipntr(1)),1, workd(ipntr(2)), 1)
			else
				continu=.false.
			end if
		end do
		if(ido.ne.99)then
			write(*,*)"znaupd do not success"
			write(*,*)info
			stop
		end if
		call zneupd(outvector,'A',selework,d,v,n,s0,workev,'I',n,WHICH,num_of_eig,tol,resid,ncv,v,n,iparam,&
						ipntr,workd,workl,lworkl,rwork,ierr)
		if(ierr.ne.0)then
			write(*,*)"zneupd do not success"
			write(*,*)ierr
			stop
		end if
		if(outvector) then
			allocate(res_eigvalue)
			allocate(res_eigenvector)
			res_eigvalue%Ten=d(1:num_of_eig)
			res_eigenvector%Ten=v(:,1:num_of_eig)
			call push_back(resultlink,res_eigvalue)
			call push_back(resultlink,res_eigenvector)
		else
			allocate(res_eigvalue)
			res_eigvalue%Ten=d(1:num_of_eig)
			call push_back(resultlink,res_eigvalue)
		end if
		return
	end subroutine
	subroutine eigen_max_min2(resultlink,T,WHICH,num_of_eig,outvector)
		type(Tensorlink),intent(inout)::resultlink
		type(Tensor),intent(in)::T
		Character*2,intent(in)::WHICH
		logical,intent(in)::outvector
		integer,intent(in)::num_of_eig
		complex*16,allocatable::A(:),v(:,:),resid(:),d(:),Cdata(:),workdata(:)
		complex*16,allocatable::workd(:),workev(:),workl(:)
		complex*16::s0
		real*8::tol
		integer::ido,info,ishfts,maxitr,mode,iparam(11), ipntr(14)
		integer,allocatable::ipiv(:)
		logical::continu
		real*8,allocatable::rwork(:)
		logical,allocatable::selework(:)
		integer,allocatable::IPIV2(:)
		integer::icouter,info2
	!	complex*16,allocatable::C0(:,:),C(:,:)
		integer::i,j,n,lworkl,ierr,ncv
		type(Tensornode),pointer::res_eigvalue,res_eigenvector
		
		if(getRank(T).ne.2) then
			write(*,*)"ERROR in dimension of input Tensor,eigen"
			write(*,*)"stop"
			stop
		end if
		n=T.dim.1
		ncv=(2+num_of_eig+n)/2
		if(ncv.gt.max_ncv) then
			ncv=max_ncv
		end if
		if(n.ne.(T.dim.2)) then
		   write(*,*)"ERROR in dimension of input Tensor,eigen"
		   write(*,*)"Should be a matrix of n * n"
			write(*,*)"stop"
			stop
		end if
		call cleanlink(resultlink)
     ! allocate(C0(n,n))
      !C0=0.
     ! do i=1,n
	!		C0(i,i)=s0
   !   end do
   !   C=T
   !   C0=C-C0
      
      lworkl = 3*ncv**2+5*ncv
		allocate(workd(3*n))
		allocate(workev(2*ncv))
		allocate(workl(lworkl))
		allocate(ipiv(n))
		allocate(rwork(ncv))
		allocate(v(n,ncv))
		allocate(resid(n))
		allocate(d(num_of_eig+1))
		allocate(IPIV2(n))
		allocate(workdata(n))
		allocate(selework(ncv))
		tol    = 0.0 
		ido    = 0
		info   = 0
		info2=0
				
				
		ishfts = 1
		maxitr = 300
		mode   = 1
		iparam(1) = ishfts 
		iparam(3) = maxitr 
		iparam(7) = mode 	 
		continu=.true.
		
		do while(continu)
			call znaupd(ido,'I',n,WHICH,num_of_eig,tol,resid,ncv,v,n,iparam,ipntr,workd,workl,lworkl,rwork,info)
			if (ido .eq. -1 .or. ido .eq. 1 ) then
				call ZGEMV('N',n,n,dcmplx(1d0,0d0),T%Tensor_Data,n,workd(ipntr(1):ipntr(1)+n-1),&
						1,dcmplx(0d0,0d0),workd(ipntr(2):ipntr(2)+n-1),1)
				!call zcopy ( n, workd(ipntr(1)),1, workd(ipntr(2)), 1)
			else
				continu=.false.
			end if
		end do
		if(ido.ne.99)then
			write(*,*)"znaupd do not success"
			write(*,*)info
			stop
		end if
		call zneupd(outvector,'A',selework,d,v,n,s0,workev,'I',n,WHICH,num_of_eig,tol,resid,ncv,v,n,iparam,&
						ipntr,workd,workl,lworkl,rwork,ierr)
		if(ierr.ne.0)then
			write(*,*)"zneupd do not success"
			write(*,*)ierr
			stop
		end if
		if(outvector) then
			allocate(res_eigvalue)
			allocate(res_eigenvector)
			res_eigvalue%Ten=d(1:num_of_eig)
			res_eigenvector%Ten=v(:,1:num_of_eig)
			call push_back(resultlink,res_eigvalue)
			call push_back(resultlink,res_eigenvector)
		else
			allocate(res_eigvalue)
			res_eigvalue%Ten=d(1:num_of_eig)
			call push_back(resultlink,res_eigvalue)
		end if
		return
	end subroutine
	subroutine eigen_max_min3(resultlink,T,WHICH,num_of_eig)!outvector=.ture.
		type(Tensorlink),intent(inout)::resultlink
		type(Tensor),intent(in)::T
		Character*2,intent(in)::WHICH
		integer,intent(in)::num_of_eig
		call eigen_max_min2(resultlink,T,WHICH,num_of_eig,.false.)
		return
	end subroutine
!*****************  eye  *****************
	type(Tensor) function eye_com(s,m,n)
		complex*16,intent(in) :: s(:)
		integer,intent(in) :: m,n
		complex*16,allocatable :: sdata(:,:)
		integer :: i,j,lens
		lens=size(s)
		allocate(sdata(m,n))
		do j=1,n
			do i=1,m
				if(i.eq.j) then
					if(lens.lt.i) then
						sdata(i,j)=dcmplx(0d0,0d0)
					else
						sdata(i,j)=s(i)
					end if
				else
					sdata(i,j)=0
				end if
			end do
		end do
		call storeTenData(eye_com,sdata)
		eye_com%rank=2
		eye_com%totalData=m*n
		eye_com%TenDim=(/m,n/)
		eye_com%flag=.true.
		return
	end function
	type(Tensor) function eye_Ten1(T)
		type(Tensor),intent(in) :: T
		complex*16,allocatable :: sdata(:,:)
		integer :: i,j,lens,n
		if(T%rank.ne.1) then
			write(*,*)"ERROR in eye_Ten,input should be a vector"
			stop
		end if
		n=T%TotalData
		allocate(sdata(n,n))
		sdata=dcmplx(0d0,0d0)
		do i=1,n
			sdata(i,i)=T%Tensor_Data(i)
		end do
		call storeTenData(eye_Ten1,sdata)
		eye_Ten1%rank=2
		eye_Ten1%totalData=n*n
		eye_Ten1%TenDim=(/n,n/)
		eye_Ten1%flag=.true.
		return
	end function
	type(Tensor) function eye_Ten2(T,m,n)
		type(Tensor),intent(in) :: T
		integer,intent(in) :: m,n
		complex*16,allocatable :: sdata(:,:)
		integer :: i,j,lens,k
		if(T%rank.ne.1) then
			write(*,*)"ERROR in eye_Ten,input should be a vector"
			stop
		end if
		lens=T%TotalData
		allocate(sdata(m,n))
		sdata=dcmplx(0d0,0d0)
		do i=1,min(m,n)
			if(lens.gt.i) then
				sdata(i,i)=T%Tensor_Data(i)
			end if
		end do
		call storeTenData(eye_Ten2,sdata)
		eye_Ten2%rank=2
		eye_Ten2%totalData=m*n
		eye_Ten2%TenDim=(/m,n/)
		eye_Ten2%flag=.true.
		return
	end function
!*****************  eye_com  *****************
	type(Tensor) function eye_real(s,m,n)
		real*8,intent(in) :: s(:)
		integer,intent(in) :: m,n
		complex*16,allocatable :: sdata(:,:)
		integer :: i,j,lens
		lens=size(s)
		allocate(sdata(m,n))
		do j=1,n
			do i=1,m
				if(i.eq.j) then
					if(lens.lt.i) then
						sdata(i,j)=dcmplx(0d0,0d0)
					else
						sdata(i,j)=dcmplx(s(i),0)
					end if
				else
				sdata(i,j)=dcmplx(0,0)
				end if
			end do
		end do
		call storeTenData(eye_real,sdata)
		eye_real%rank=2
		eye_real%totalData=m*n
		eye_real%TenDim=(/m,n/)
		eye_real%flag=.true.
		return
	end function
!*****************  one_com  *****************
	type(Tensor) function one_com(m,n)
		integer,intent(in) :: m,n
		complex*16,allocatable :: sdata(:,:)
		integer :: i,j
		allocate(sdata(m,n))
		do j=1,n
			do i=1,m
				if(i.eq.j) then
					sdata(i,j)=EE
				else
					sdata(i,j)=dcmplx(0,0)
				end if
			end do
		end do
		call storeTenData(one_com,sdata)
		one_com%rank=2
		one_com%totalData=m*n
		one_com%TenDim=(/m,n/)
		one_com%flag=.true.
		return
	end function	
!****************  one_Ten  *************
	type(Tensor) function one_Ten(dimen)
		integer,intent(in) :: dimen(:)
		integer :: i,total,lenDim,min_dim,address,dim_i
		complex*16,allocatable::TenData(:)
		total=product(dimen)
		lenDim=size(dimen)
		min_dim=minval(dimen)
		allocate(TenData(total))
		TenData=0
		do dim_i=1,min_dim
			address=0
			do i=lenDim,2,-1
				address=address+(dim_i-1)*product(dimen(1:(i-1)))
			end do
			address=address+dim_i
			TenData(address)=dcmplx(1d0,0d0)
		end do
		one_Ten=buildTen(dimen,TenData)
		return
	end function	
! return exp(H)
	type(Tensor) function expm(H)
		type(Tensor),intent(in) ::H
		complex*16,allocatable :: Hdata(:,:),work(:)
		complex*16,allocatable :: expHdata(:,:)
		complex*16,allocatable :: engVal(:),engVec(:,:)
		integer :: hdim(2),i,j,N,lwork,info,SDIM
		real*8,allocatable::rwork(:)
		logical,allocatable::bwork(:)
		type(dimension)::newHdim
		hdim(1)=H.dim.1
		hdim(2)=H.dim.2
		newHdim=H%TenDim
		if(getRank(H).ne.2) then
			write(*,*)"ERROR in expm"
			write(*,*)"input Tensor should be a matrix"
			stop
		end if
		if(hdim(1).ne.hdim(2)) then
			write(*,*)"ERROR in expm"
			stop
		end if
		N=hdim(1)
		lwork=2*N
		
		allocate(Hdata(N,N))
		allocate(engVal(N))
		allocate(engVec(N,N))
		allocate(rwork(N))
		allocate(work(lwork))
		allocate(bwork(N))
		allocate(expHdata(N,N))
		Hdata=H
		call ZGEES('V','N',1,N,Hdata,N,SDIM,engVal,engVec,N,WORK,LWORK,RWORK,BWORK,INFO)
		do i=1,N
			engVal(i)=exp(engVal(i))
		end do
		do i=1,N
			do j=1,N
				expHdata(i,j)=engVec(i,j)*engVal(j)
			end do
		end do
		engVec=dconjg(transpose(engVec))
		expHdata=matmul(expHdata,engVec)
		expm=set(newHdim,reshape(expHdata,(/N*N/)))
		return
	end function
!return eigenvalue and eigenvectoers		
	subroutine eng(H,vec,val)
		type(Tensor),intent(in) ::H
		type(Tensor),intent(inout) ::val,vec
		complex*16,allocatable :: Hdata(:,:),work(:)
		complex*16,allocatable :: engVal(:),engVec(:,:)
		integer :: hdim(2),i,j,N,lwork,info,SDIM
		real*8,allocatable::rwork(:)
		logical,allocatable::bwork(:)
		type(dimension)::newHdim
		hdim(1)=H.dim.1
		hdim(2)=H.dim.2
		newHdim=H%TenDim
		if(getRank(H).ne.2) then
			write(*,*)"ERROR in eng"
			write(*,*)"input Tensor should be a matrix"
			stop
		end if
		if(hdim(1).ne.hdim(2)) then
			write(*,*)"ERROR in eng"
			stop
		end if
		N=hdim(1)
		lwork=2*N
		
		allocate(Hdata(N,N))
		allocate(engVal(N))
		allocate(engVec(N,N))
		allocate(rwork(N))
		allocate(work(lwork))
		allocate(bwork(N))
		Hdata=H
		call ZGEES('V','N',1,N,Hdata,N,SDIM,engVal,engVec,N,WORK,LWORK,RWORK,BWORK,INFO)
		val=eye(engVal,N,N)
		vec=set((/N,N/),reshape(engVec,(/N*N/)))
		return
	end subroutine	
!return eigenvalue,eigenvalue store in a vector
	subroutine eng2(H,val)
		type(Tensor),intent(in) ::H
		type(Tensor),intent(inout) ::val
		complex*16,allocatable :: Hdata(:,:),work(:)
		complex*16,allocatable :: engVal(:),engVec(:,:)
		integer :: hdim(2),i,j,N,lwork,info,SDIM
		real*8,allocatable::rwork(:)
		logical,allocatable::bwork(:)
		type(dimension)::newHdim
		hdim(1)=H.dim.1
		hdim(2)=H.dim.2
		newHdim=H%TenDim
		if(getRank(H).ne.2) then
			write(*,*)"ERROR in eng"
			write(*,*)"input Tensor should be a matrix"
			stop
		end if
		if(hdim(1).ne.hdim(2)) then
			write(*,*)"ERROR in eng"
			stop
		end if
		N=hdim(1)
		lwork=2*N
		
		allocate(Hdata(N,N))
		allocate(engVal(N))
		allocate(engVec(N,N))
		allocate(rwork(N))
		allocate(work(lwork))
		allocate(bwork(N))
		Hdata=H
		call ZGEES('V','N',1,N,Hdata,N,SDIM,engVal,engVec,N,WORK,LWORK,RWORK,BWORK,INFO)
		val=engVal
		return
	end subroutine	
!*****************  sqrt  *********************
! return sqrt(A),A is a matrix

	type(Tensor) function sqrtTen(T)
		type(Tensor),intent(in) :: T
		type(Tensor)::val,vec
		complex*16,allocatable :: engVal(:),engVec(:,:)
		complex*16,allocatable :: Tdata(:,:),work(:)
		integer :: Tdim(2),i,j,N,lwork,info,SDIM
		real*8,allocatable::rwork(:)
		logical,allocatable::bwork(:)
		type(dimension)::newHdim
		Tdim(1)=T.dim.1
		Tdim(2)=T.dim.2
		newHdim=T%TenDim
		if(getRank(T).ne.2) then
			write(*,*)"ERROR in sqrt of Tensor"
			write(*,*)"input Tensor should be a matrix"
			stop
		end if
		if(Tdim(1).ne.Tdim(2)) then
			write(*,*)"ERROR in sqrt of Tensor"
			write(*,*)"input Tensor should be a matrix"
			write(*,*)Tdim(1),Tdim(2)
			stop
		end if
		N=Tdim(1)
		lwork=2*N
		
		allocate(Tdata(N,N))
		allocate(engVal(N))
		allocate(engVec(N,N))
		allocate(rwork(N))
		allocate(work(lwork))
		allocate(bwork(N))
		Tdata=T
		call ZGEES('V','N',1,N,Tdata,N,SDIM,engVal,engVec,N,WORK,LWORK,RWORK,BWORK,INFO)
		do i=1,N
			if((dreal(engVal(i)).ge.0).and.(dabs(aimag(engVal(i))).lt.1d-8)) then
				engVal(i)=dsqrt(real(engVal(i)))
			else
				engVal(i)=0d0
				!write(*,*)"ERROR in sqrt of Tensor"
				write(*,*)engVal(i)
				!write(*,*)Tdata
				!stop
			end if
		end do
		val=eye(engVal,N,N)
		vec=set((/N,N/),reshape(engVec,(/N*N/)))
		sqrtTen=vec * val * (.H. vec)
		return
	end function	
!*****************  sqrt  *********************
! return trace(A),A is a matrix

	complex*16 function trace(T)
		type(Tensor),intent(in) :: T
		integer::rank,i
		rank=getRank(T)
		if(rank.ne.2) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			stop
		end if
		if((T.dim.1).ne.(T.dim.2)) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			write(*,*)(T.dim.1),(T.dim.2)
			stop
		end if
		trace=dcmplx(0d0,0d0)
		do i=1,(T.dim.1)
			trace=trace+(T.i.(/i,i/))
		end do
		return
	end function	
		
		
!*****************  Htranspose  *****************
!		note: update at 2013.10.29
	type(Tensor) function Htranspose(T)
		type(Tensor),intent(in) :: T
		complex*16,allocatable :: Tdata(:,:),newdata(:,:)
		integer::rank,m,n
		type(dimension)::dimen
		rank=T%rank
		if(rank.eq.1) then
			m=T.dim. 1
			allocate(Htranspose%Tensor_Data(m))
			Htranspose%Tensor_Data=dconjg(T%Tensor_Data)
			Htranspose%rank=T%rank
			Htranspose%totalData=T%totalData
			Htranspose%TenDim=T%TenDim
			Htranspose%flag=.true.
		else if(rank.eq.2) then
			m=T.dim.1
			n=T.dim.2
			Dimen=(T%TenDim.sub.2)+(T%TenDim.sub.1)
			allocate(Tdata(m,n))
			allocate(newdata(n,m))
			Tdata=T
			newdata=dconjg(transpose(Tdata))
			Htranspose%TenDim=Dimen
			Htranspose%rank=2
			Htranspose%flag=.true.
			Htranspose%TotalData=m*n
			call storeTenData(Htranspose,newdata)	
		else
			write(*,*) "Tensor should be 1 or 2 dimension"
		end if
		return
	end function
!*****************  conjugate  *****************
!		note: update at 2013.10.29
	type(Tensor) function conjugate(T)
		type(Tensor),intent(in) :: T
		complex*16,allocatable :: newdata(:)
		integer :: m
		m=T%totalData
		allocate(conjugate%Tensor_Data(m))
		conjugate%Tensor_Data=dconjg(T%Tensor_Data)
		conjugate%rank=T%rank
		conjugate%totalData=T%totalData
		conjugate%TenDim=T%TenDim
		conjugate%flag=.true.
		return
	end function
!*****************  conjugateTranspose  *****************
	type(Tensor) function conjugateTranspose(T)
		type(Tensor),intent(in) :: T
		complex*16,allocatable :: Tdata(:,:),newdata(:,:)
		integer :: m,n
		type(Dimension) :: Dimen
		if(T%rank.eq.2) then
			m=T.dim.1
			n=T.dim.2
			Dimen=(T%TenDim.sub.2)+(T%TenDim.sub.1)
			allocate(Tdata(m,n))
			allocate(newdata(n,m))
			Tdata=T
			newdata=dconjg(transpose(Tdata))
			conjugateTranspose%TenDim=Dimen
			conjugateTranspose%rank=2
			conjugateTranspose%flag=.true.
			conjugateTranspose%TotalData=m*n
			call storeTenData(conjugateTranspose,newdata)			
		else
			write(*,*) "The input Tensor should 1 or 2 dimension"
		end if
		return
	end function
!***************   realTen  **********************
	type(DTensor) function realTen(T)
		type(Tensor),intent(in)::T
		realTen=DbuildTen8(T%TenDim,dreal(T%Tensor_Data))
		return
	end function	
!***************   imagTen  **********************
	type(DTensor) function imagTen(T)
		type(Tensor),intent(in)::T
		imagTen=DbuildTen8(T%TenDim,dimag(T%Tensor_Data))
		return
	end function		
!*****************  value  *****************
!			return the value of <phi|F|phi>,where F is an operator
!			F=<I',s1',s2',..|F|I,s1,s2,..>
!			update at 2013.10.29
	complex*16 function value(F_,Tr_)
		type(Tensor),intent(in) :: F_,Tr_
		type(Tensor) :: Tr,Tl,F,midResult
		integer :: i,rank,oper(2,3)
		rank=Tr_%rank
		Tr=Tr_.c.(/1,rank/)
		oper(1,:)=(/1,1,rank/)
		oper(2,:)=(/1,2,rank+1/)
		F=F_.cd.oper
		Tl=.h.Tr
		value=Tl*F*Tr
		return
	end function			
!***************   dot   ******************
!		return <phi1|phi2>
!	
	complex*16 function dot(phi1,phi2)
		Type(Tensor),intent(in)::phi1,phi2
		integer::N1,N2
		N1=getTotalData(phi1)
		N2=getTotalData(phi2)
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			stop
		end if
		dot=ZDOTC(N1,phi1%Tensor_Data,1,phi2%Tensor_Data,1)
		RETURN
	end function

!*****************  norm   *****************
!			return  <phi|phi>
!			update at 2013.10.29
	real*8 function norm2(Tr)
		type(Tensor),intent(in) :: Tr
		norm2=DZNRM2(Tr%TotalData,Tr%Tensor_Data,1)
		norm2=norm2*norm2
		return
	end function	
!		return  sqrt(<phi|phi>	)
	real*8 function norm(Tr)
		type(Tensor),intent(in) :: Tr
		norm=DZNRM2(Tr%TotalData,Tr%Tensor_Data,1)
		return
	end function		
!*****************  addressToIndes   *****************
	integer function addressToIndes(T,Adim)
		type(Tensor),intent(in) :: T
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
		
! 		modify the code below in 2013.5.18		
!		addressToIndes=Adim(1)
!		Tdim=T%TenDim
!		do i=2,size(Adim)
!			addressToIndes=addressToIndes+(Adim(i)-1)*Tdim(i-1)
!		end do
!		return
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
! 		modify the code below in 2013.5.18
!		addressToIndes2=Adim(1)
!		do i=2,size(Adim)  
!			addressToIndes2=addressToIndes2+(Adim(i)-1)*Tdim(i-1)
!		end do
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
!*****************  Element   *****************
	complex*16 function Element(T,Tdim)
		type(Tensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		integer::Inde
		inde=addressToIndes(T,Tdim)
		Element=T%Tensor_Data(inde)
		return
	end function		  
	complex*16 function Element2(T,inde)
		type(Tensor),intent(in) ::T
		integer,intent(in)::inde
		Element2=T%Tensor_Data(inde)
		return
	end function
!		inde=[-1,inde_min,inde_max] output data(:,inde_min:inde_max)
!	or[-2,inde_min,inde_max],data(inde_min:inde_max,:)
!	or [-1,inde_row] [-2,inde_col],output row or col

	type(Tensor) function subTen(T,inde)	
		type(Tensor),intent(in) ::T
		integer,intent(in)::inde(:)
		complex*16,allocatable::Tdata(:,:),subTen1(:),subTen2(:,:)
		integer::num
		if(getRank(T).ne.2)then
			write(*,*)"error in subTen,only matrix is allowed"
			stop
		end if
		Tdata=T
		if(size(inde).eq.2) then
			select case (inde(1))
				case (-1)!output row
					allocate(subTen1(T.dim.2))
					subTen1=Tdata(inde(2),:)
				case (-2)!
					allocate(subTen1(T.dim.1))
					subTen1=Tdata(:,inde(2))
				case default 
					write(*,*) "no such case in subTen"
					write(*,*)inde
					stop
			end 	select
			subTen=subTen1
			return
		end if
		if(size(inde).ne.3) then
			write(*,*) "no such case in subTen"
			write(*,*) "length of inde is ",size(inde)
			write(*,*)inde
			stop
		end if
		num=inde(3)-inde(2)+1
		if(num.le.0) then
			write(*,*)"error subTen"
			stop
		end if
		select case (inde(1))
			case (-1)!output row
				allocate(subTen2(num,T.dim.2))
				subTen2=Tdata(inde(2):inde(3),:)
			case (-2)!
				allocate(subTen2(T.dim.1,num))
				subTen2=Tdata(:,inde(2):inde(3))
			case default 
				write(*,*) "no such case in subTen"
				write(*,*)inde
				stop
		end 	select
		subTen=subTen2
		return
	end function
!***************  directProduct  *********************
	type(Tensor) function directProduct(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		complex*16,allocatable :: T1data(:,:),T2data(:,:)
		complex*16,allocatable :: newdata(:,:,:,:)
		integer :: m1,n1,m2,n2,i,j,k,l,rank(2)
		type(Dimension):: D1,Dtemp
		if((getRank(T1).eq.1).and.(getRank(T2).eq.1)) then
			directProduct=directProduct1(T1,T2)
			return
		end if
		if((getRank(T1).eq.2).and.(getRank(T2).eq.2)) then
			D1=T1%TenDim.sub.1
			Dtemp=T2%TenDim.sub.1
			D1=D1+Dtemp
			Dtemp=T1%TenDim.sub.2
			D1=D1+Dtemp
			Dtemp=T2%TenDim.sub.2
			D1=D1+Dtemp
			m1=T1.dim.1
			n1=T1.dim.2
			m2=T2.dim.1
			n2=T2.dim.2
			allocate(T1data(m1,n1))
			allocate(T2data(m2,n2))
			allocate(newdata(m1,m2,n1,n2))
			T1data=reshape(T1%Tensor_Data,(/m1,n1/))
			T2data=reshape(T2%Tensor_Data,(/m2,n2/))
			do l=1,n2
				do k=1,n1
					 do j=1,m2
						 do i=1,m1
							newdata(i,j,k,l)=T1data(i,k)*T2data(j,l)
						end do
					end do
				end do
			end do
			directProduct=set(D1,reshape(newdata,(/m1*m2*n1*n2/)))
		else
			write(*,*)"The directProduct is just for matrix"
			stop
		end if
	return
	end function
!********************  directProductM  ************************
	type(Tensor) function directProductM(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		complex*16,allocatable :: T1data(:,:),T2data(:,:)
		complex*16,allocatable :: newdata(:,:,:,:)
		integer :: m1,n1,m2,n2,i,j,k,l
		type(Dimension):: D1,Dtemp
		if((getRank(T1).eq.1).and.(getRank(T2).eq.1)) then
			directProductM=directProduct1V(T1,T2)
			return
		end if
		if((getRank(T1).eq.2).and.(getRank(T2).eq.2)) then
			D1=T1%TenDim.sub.1
			Dtemp=T2%TenDim.sub.1
			D1=D1+Dtemp
			Dtemp=T1%TenDim.sub.2
			D1=D1+Dtemp
			Dtemp=T2%TenDim.sub.2
			D1=D1+Dtemp
			D1=DimConstract(D1,1,1)
			D1=DimConstract(D1,2,1)
			m1=T1.dim.1
			n1=T1.dim.2
			m2=T2.dim.1
			n2=T2.dim.2
			allocate(T1data(m1,n1))
			allocate(T2data(m2,n2))
			allocate(newdata(m1,m2,n1,n2))
			T1data=reshape(T1%Tensor_Data,(/m1,n1/))
			T2data=reshape(T2%Tensor_Data,(/m2,n2/))
			do l=1,n2
				do k=1,n1
					do j=1,m2
						do i=1,m1
							newdata(i,j,k,l)=T1data(i,k)*T2data(j,l)
						end do
					end do
				end do
			end do
			directProductM=set(D1,reshape(newdata,(/m1*m2*n1*n2/)))
		else
			write(*,*)"The directProductM is just for matrix"
			stop
		end if
		return
	end function
!************* directProduct1 **************************
!
!       for two Tensor whose ranks are 1		
!
	type(Tensor) function directProduct1V(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		complex*16,allocatable :: newdata(:,:)
		integer :: m1,n1,m2,n2,i,j,k,l
		type(Dimension):: D1,Dtemp
		if((getRank(T1).eq.1).and.(getRank(T2).eq.1)) then
			D1=T1%TenDim
			D1=D1+T2%TenDim
			D1=DimConstract(D1,1,1)
			n1=T1.dim.1
			n2=T2.dim.1
			allocate(newdata(n1,n2))
			do l=1,n2
				do k=1,n1
					newdata(k,l)=T1%Tensor_Data(k)*T2%Tensor_Data(l)
				end do
			end do
			directProduct1V=set(D1,reshape(newdata,(/n1*n2/)))
		else
			write(*,*)"The directProduct1 is just for vector"
			stop
		end if
		return
	end function	
	type(Tensor) function directProduct1(T1,T2)
		type(Tensor),intent(in) :: T1,T2
		complex*16,allocatable :: newdata(:,:)
		integer :: m1,n1,m2,n2,i,j,k,l
		type(Dimension):: D1,Dtemp
		if((getRank(T1).eq.1).and.(getRank(T2).eq.1)) then
			D1=T1%TenDim
			D1=D1+T2%TenDim
			n1=T1.dim.1
			n2=T2.dim.1
			allocate(newdata(n1,n2))
			do l=1,n2
				do k=1,n1
					newdata(k,l)=T1%Tensor_Data(k)*T2%Tensor_Data(l)
				end do
			end do
			directProduct1=set(D1,reshape(newdata,(/n1*n2/)))
		else
			write(*,*)"The directProduct1 is just for vector"
			stop
		end if
		return
	end function
!*******************  direct sum  **********************
	type(Tensor) function  direct_sum(T1,T2) result(T)
		type(Tensor),intent(in) :: T1,T2
		type(Tensor) ::temp1,temp2
		type(Dimension) ::dimT1,dimT2
		integer :: m,n
		dimT1=T1%TenDim
		dimT2=T2%TenDim
		if(dimT1.equ.dimT2) then
			temp1=T1.mxx.one_com(m,n)
			temp2=one_com(m,n).mxx.T2
			T=temp1+temp2
			T=.dc.T
		else
			write(*,*) "error in direct_sumM"
		end if
		return 
	end function
!****************  direct sum returning matrix *******************
	type(Tensor) function  direct_sumM(T1,T2) result(T)
		type(Tensor),intent(in) :: T1,T2
		type(Tensor) ::temp1,temp2
		type(Dimension) ::dimT1,dimT2
		integer :: m,n
		dimT1=T1%TenDim
		dimT2=T2%TenDim
		if(dimT1.equ.dimT2) then
			m=T1.dim.1
			n=T1.dim.2
			temp1=T1.mxx.one_com(m,n)
			temp2=one_com(m,n).mxx.T2
			T=temp1+temp2
		else
			write(*,*) "error in direct_sumM"
		end if
		return 
	end function
!**************   conbine   ********************
!
!                    / T1 \
!	conbine(T1,T2)	=  |----|
!		               \ T2 /
!	T1 :a [...,m,n,l] matrix
!	T2 :a [...,m,n,l] matrix
!	conbine(T1,T2):a [...,m,n,l,2] matrix
!or 
!	T1 :a [...,m,n,l] matrix
!	T2 :a [...,m,n] matrix
!	conbine(T1,T2):a [...,m,n,l+1] matrix
	type(Tensor) function conbine(T1,T2)
		type(Tensor),intent(in)::T1,T2
		integer,allocatable::dim1(:),dim2(:)
		integer::total,total1,i,dim_n
		complex*16,allocatable::newdata(:),olddata(:,:)
		type(Dimension)::newDim,dimen1,dimen2
		if(T1%rank.eq.T2%rank) then
			dim1=.subDim.T1
			dim2=.subDim.T2
			if(.not.(dim1.equ.dim2)) then
				write(*,*)"can not conbie two Tensor"
				write(*,*)dim1
				write(*,*)dim2
				write(*,*)"stop"
				stop
				return
			end if
			total1=gettotalData(T1)
			total=total1+total1
			allocate(newdata(total1*2))
			newdata(:total1)=T1%Tensor_Data
			newdata(total1+1:)=T2%Tensor_Data
			call storeTenData(conbine,newdata)
			conbine%Rank=getRank(T1)+1
			newDim=.subDim.T1
			newDim=newDim+(/2/)
			conbine%TenDim=newDim
			conbine%totalData=total
			conbine%flag=.true.
			return
		end if
		dimen1=.subDim.T1
		dimen1=DimConstract(dimen1,1,T1%rank-2)
		dimen2=.subDim.T2
		dimen2=DimConstract(dimen2,1,T2%rank)
		if(.not.((dimen1.sub.1).equ.dimen2)) then
			write(*,*)"can not conbie two Tensor"
		!	write(*,*)T1%rank,T2%rank
			write(*,*)"stop2"
			stop
			return
		end if
		dim_n=dimen1.i.2
		total1=gettotalData(T1)
		total=gettotalData(T2)
		allocate(newdata(total+total1))
		!allocate(olddata(dimen1.i.1,dimen1.i.2))
		!call ZCOPY(gettotalData(T1),T1%Tensor_Data,1,olddata,1)
		newdata(1:total1)=T1%Tensor_Data!olddata
		newdata(total1+1:)=T2%Tensor_Data
		call storeTenData(conbine,newdata)
		conbine%Rank=getRank(T1)
		newDim=(.subDim.T2)
		newDim=newDim+(/dim_n+1/)
		conbine%TenDim=newDim
		conbine%totalData=total+total1
		conbine%flag=.true.
	end function
	type(Tensor) function conbinerow(T1,T2)
		type(Tensor),intent(in)::T1,T2
		integer,allocatable::dim1(:),dim2(:)
		integer::total,total1,i,dim_n
		complex*16,allocatable::newdata(:,:),olddata(:,:)
		type(Dimension)::newDim,dimen1,dimen2
		if(T1%rank.eq.T2%rank) then
			dim1=.subDim.T1
			dim2=.subDim.T2
			if(.not.(dim1.equ.dim2)) then
				write(*,*)"can not conbinerow two Tensor"
				write(*,*)dim1
				write(*,*)dim2
				write(*,*)"stop"
				stop
				return
			end if
			total1=gettotalData(T1)
			total=total1+total1
			allocate(newdata(2,total1))
			newdata(1,:)=T1%Tensor_Data
			newdata(2,:)=T2%Tensor_Data
			call storeTenData(conbinerow,newdata)
			conbinerow%Rank=getRank(T1)+1
			newDim=.subDim.T1
			newDim=(/2/)+newDim
			conbinerow%TenDim=newDim
			conbinerow%totalData=total
			conbinerow%flag=.true.
			return
		end if
		dimen1=.subDim.T1
		dimen1=DimConstract(dimen1,2,T1%rank)
		dimen2=.subDim.T2
		dimen2=DimConstract(dimen2,1,T2%rank)
		if(.not.((dimen1.sub.2).equ.dimen2)) then
			write(*,*)"can not conbinerow two Tensor"
		!	write(*,*)T1%rank,T2%rank
			write(*,*)"stop2"
			stop
			return
		end if
		dim_n=dimen1.i.1
		total1=gettotalData(T1)
		total=gettotalData(T2)
		allocate(newdata(dim_n+1,total))
		allocate(olddata(dimen1.i.1,dimen1.i.2))
		call ZCOPY(gettotalData(T1),T1%Tensor_Data,1,olddata,1)
		newdata(1:dim_n,:)=olddata
		newdata(dim_n+1,:)=T2%Tensor_Data
		call storeTenData(conbinerow,newdata)
		conbinerow%Rank=getRank(T1)
		newDim=(.subDim.T2)
		newDim=(/dim_n+1/)+newDim
		conbinerow%TenDim=newDim
		conbinerow%totalData=total+total1
		conbinerow%flag=.true.
	end function			
				
!****************  Tensor1  ************************
!	dimen is the dimension of the Tensor
!	return the Tensor with all the data are num
!	num is a complex*16 type
	type(Tensor)function Tensor1(dimen,num)
		type(Dimension),intent(in) :: dimen
		complex*16,intent(in) :: num
		integer::total,i
		total=1
		do i=1,DimSize(dimen)
			total=total*(dimen.i.i)
		end do
		allocate(Tensor1%Tensor_data(total))
		do i=1,total
			Tensor1%Tensor_data(i)=num
		end do
		Tensor1%rank=DimSize(dimen)
		Tensor1%totalData=total
		Tensor1%TenDim=dimen
		Tensor1%flag=.true.
		return
	end function
!******************  modify the data in the Tensor  **************
! modify the data in dimen of the Tensor
	subroutine modifyTen_val(Ten,dimen,val)
		type(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		complex*16,intent(in)::val
		integer::addre
		addre=addressToIndes(Ten,dimen)
		Ten%Tensor_Data(addre)=val
		return
	end subroutine
! modify all the data in the Tensor
	subroutine modifyTen_dat(Ten,val)
		type(Tensor),intent(inout)::Ten
		complex*16,intent(in)::val(:)
		if(size(val).eq.Ten%TotalData) then
			write(*,*)"ERROR in modify Tensor data"
			write(*,*)"stop"
			stop
		end if
		!Ten%Tensor_Data=val
		call assignments(Ten%Tensor_Data,val)
		return
	end subroutine
! modify all the data in the Tensor
	subroutine modifyTen_dat2(Ten,val)
		type(Tensor),intent(inout)::Ten
		complex*16,intent(in)::val(:,:)
		if(Ten%rank.ne.2) then
			write(*,*)"ERROR in modify Tensor data"
			write(*,*)"The rank is 2,stop"
			stop
		end if
		if(product(val).eq.Ten%TotalData) then
			write(*,*)"ERROR in modify Tensor data"
			write(*,*)"stop"
			stop
		end if
		call storeTenData(Ten,val)
		return
	end subroutine	
! modify all the data in the Tensor
	subroutine modifyTen_dat3(Ten,val)
		type(Tensor),intent(inout)::Ten
		complex*16,intent(in)::val(:,:,:)
		if(Ten%rank.ne.3) then
			write(*,*)"ERROR in modify Tensor data"
			write(*,*)"The rank is 3,stop"
			stop
		end if
		if(product(val).eq.Ten%TotalData) then
			write(*,*)"ERROR in modify Tensor data"
			write(*,*)"stop"
			stop
		end if
		call storeTenData(Ten,val)
		return
	end subroutine	
! modify all the dimension in the Tensor
	subroutine resetdim1(Ten,dimen)
		type(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		type(Dimension)::dimensi
		if(product(dimen).ne.Ten%TotalData) then
			write(*,*)"ERROR in resetdim"
			write(*,*)product(dimen),Ten%TotalData
			write(*,*)"stop"
			stop
		end if
		dimensi=dimen
		Ten%TenDim=dimensi
		Ten%rank=size(dimen)
		return
	end subroutine	
	subroutine resetdim2(Ten,dimen)
		type(Tensor),intent(inout)::Ten
		type(dimension),intent(in)::dimen
		integer,allocatable::dimensi(:)
		dimensi=dimen
		if(product(dimensi).ne.Ten%TotalData) then
			write(*,*)"ERROR in resetdim"
			write(*,*)product(dimensi),Ten%TotalData
			write(*,*)"stop"
			stop
		end if
		Ten%TenDim=dimen
		Ten%rank=size(dimensi)
		return
	end subroutine	
			
!******************  Ten_Ten  *********************
!		contract The i1 index of T1 and the i2 index of t2
!		T1: (D1,D2,..D_i1,...),T2 :(E1,E2,...,E_i2,....),then the result
!	  will be T:(D1,D2,..D_{i1-1},D_{i1+1},...,D_{rank1},E1,E2,...,E_{i2-1},E_{i2+1},...,E_{rank2})	
!		if if Decompose=1, T will be Decompose
	type(Tensor) function Ten_Ten(T1_,T2_,i1,i2,ifDecompose) result(T)
		type(Tensor),intent(in) :: T1_,T2_
		integer,intent(in) :: i1,i2,ifDecompose
		integer :: i,j,k,D2,decompoD1,decompoD2!Di use for decompese
		type(Tensor) :: T1,T2
		integer :: decompoD1keep,decompoD2keep,rank1,&
			Tdim1(T1_%rank),Tdim2(T2_%rank),rank2,oper(2,3)
		rank1=T1_%rank
		rank2=T2_%rank
		T1=T1_
		T2=T2_
		!T1
		if((rank1.ne.1).and.(rank1.ne.2)) then
			if(i1.eq.1) then
				T1=contracts(T1,2,rank1)
				T1=.p.T1
			else if (i1.eq.rank1) then
				T1=contracts(T1,1,rank1-1)
			else
				oper(1,:)=(/1,1,i1-1/)
				oper(2,:)=(/1,3,rank1/)
				!T1=contracts(T1,1,i1-1)
				!T1=contracts(T1,3,rank1)
				T1=T1.cd.oper
				T1=T1.p.1
				call dimOperation(T1,(/1,1,2/))
				!T1=T1.c.1
			end if
		end if
		if(rank1.eq.2) then
			if(i1.eq.1) then
				T1=.p.T1
			end if
		end if
		!T2
		if((rank2.ne.1).and.(rank2.ne.2)) then
			if(i2.eq.1) then
				T2=contracts(T2,2,rank2)
			else if (i2.eq.rank2) then
				T2=contracts(T2,1,rank2-1)
				T2=.p.T2
			else
				oper(1,:)=(/1,1,i2-1/)
				oper(2,:)=(/1,3,rank2/)
				T2=T2.cd.oper
				!T2=contracts(T2,1,i2-1)
				!T2=contracts(T2,3,rank2)
				T2=T2.p.3
				call dimOperation(T2,(/1,2,3/))
				!T2=T2.c.2
			end if
		end if
		if(rank2.eq.2) then
			if(i2.eq.2) then
				T2=.p.T2
			end if
		end if
		T=T1 * T2
		if(ifDecompose.eq.1) then
			call dimOperation(T,(/3/))
			!T=.dc.T
		end if
		return
	end function
!*****************permutefo******************************
!		T_{1,2,3,..,i,..,n},permutefo(T,i)=_{i,1,2,3,..,i-1,i+1,..,n}
!
	type(Tensor) function permutefo(T,inde)
		type(Tensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		integer::oper(2,3)
		rank=getRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permutefo"
			write(*,*)"stop"
			stop
		end if
		if(inde.eq.1) then
			permutefo=T
			return
		end if
		if(inde.eq.rank) then
			permutefo=contracts(T,1,rank-1)
			permutefo=.p.permutefo
			permutefo%TenDim=DimDecomposeAll(permutefo%TenDim)
			permutefo%rank=DimSize(permutefo%TenDim)
			return
		end if
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,3,rank/)
		permutefo=T.cd.oper
		permutefo=permutefo.p.3
		permutefo%TenDim=DimDecomposeAll(permutefo%TenDim)
		permutefo%rank=DimSize(permutefo%TenDim)
		return
	end function
!****************** permuteInde**************************
!		T_{1,2,3,..,i,..,n},permuteInde(T,i)=_{2,3,..,i,1,i+1,..,n}
!
	type(Tensor) function permuteInde(T,inde)
		type(Tensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		integer::oper(2,3)
		rank=getRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permuteInde",inde,rank
			write(*,*)"stop"
			stop
		end if
		if(inde.eq.1) then
			permuteInde=T
			return
		end if
		if(inde.eq.rank) then
			permuteInde=contracts(T,2,rank)
			permuteInde=.p.permuteInde
			permuteInde%TenDim=DimDecomposeAll(permuteInde%TenDim)
			permuteInde%rank=DimSize(permuteInde%TenDim)
			return
		end if
		oper(1,:)=(/1,2,inde/)
		oper(2,:)=(/1,3,rank/)
		permuteInde=T.cd.oper
		permuteInde=permuteInde.p.3
		permuteInde%TenDim=DimDecomposeAll(permuteInde%TenDim)
		permuteInde%rank=DimSize(permuteInde%TenDim)
		return
	end function

!cccccccccccccccc	pauli_matrix  ccccccccccccccccccccccc
!	the output is sigma_i*num,where i=x,y,z
	subroutine pauli_matrix(Tx,Ty,Tz,num)
		type(Tensor) :: Tx,Ty,Tz
		real*8,intent(in) :: num
		complex*16 :: x(2,2)
		complex*16 :: y(2,2)
		complex*16 :: z(2,2)
			x(1,1)=dcmplx(0,0)
			x(1,2)=EE
			x(2,1)=EE
			x(2,2)=dcmplx(0,0)

			y(1,1)=dcmplx(0,0)
			y(1,2)=-II
			y(2,1)=II
			y(2,2)=dcmplx(0,0)

			z(1,1)=EE
			z(1,2)=dcmplx(0,0)
			z(2,1)=dcmplx(0,0)
			z(2,2)=-EE
			Tx=set((/2,2/),reshape(x,(/4/)))*num
			Ty=set((/2,2/),reshape(y,(/4/)))*num
			Tz=set((/2,2/),reshape(z,(/4/)))*num
		return
	end subroutine

   subroutine check_zero_Ten(T)!in case of the element of rho is -1*10^-70,this is zero
	   type(Tensor),intent(inout)::T
	   integer::i
	   do i=1,getTotalData(T)
         if(abs(T%Tensor_Data(i)).lt.zero_num) then
            T%Tensor_Data(i)=dcmplx(0d0,0d0)
	      end if
	   end do
	end subroutine

!****************************************************************************************
	 
! 	Tensorlink and Tensornode 	  





!	add a Tensor to the end of the link
!	there is no vector index in the  input Tensor
	subroutine push_backTen(h,T)
		type(Tensorlink),intent(inout) ::h
		type(Tensor),intent(in):: T
		type(Tensornode),pointer ::node
		allocate(node)
		node%Ten=T
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Tend=>node
		else
			h%length=h%length+1
			h%Tend%next=>node
			h%Tend=>node
		end if
		return
	end subroutine
!	add a Tensor to the begin of the link	
!	there is no vector index in the  input Tensor
	subroutine push_fo(h,T)
		type(Tensorlink),intent(inout) ::h
		type(Tensor),intent(in):: T
		type(Tensornode),pointer ::node,p
		allocate(node)
		node%Ten=T
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Tend=>node
		else
			p=>h%head
			h%head=>node
			node%next=>p
			h%length=h%length+1
		end if
		return
	end subroutine
!	add a Tensor to the end of the link	
!	inde is the index of the Tensor
	subroutine push_backI(h,T,inde)
		type(Tensorlink),intent(inout) ::h
		type(Tensor),intent(in):: T
		integer,intent(in) :: inde(:)
		type(Tensornode),pointer ::node
		allocate(node)
		node%Ten=T
		call assignment_int_dim1(node%indice,inde)
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Tend=>node
		else
			h%length=h%length+1
			h%Tend%next=>node
			h%Tend=>node
		end if
		return
	end subroutine
!	add a Tensor to the begin of the link	
!	inde is the index of the Tensor	
	subroutine push_foI(h,T,inde)
		type(Tensorlink),intent(inout) ::h
		type(Tensor),intent(in):: T
		integer,intent(in) :: inde(:)
		type(Tensornode),pointer ::node,p
		allocate(node)
		node%Ten=T
		call assignment_int_dim1(node%indice,inde)
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Tend=>node
		else
			p=>h%head
			h%head=>node
			node%next=>p
			h%length=h%length+1
		end if
		return
	end subroutine
!	add a Tensornode to the end of the link
!	there is no vector index in the  input Tensor
	subroutine push_backnode(h,node)
		type(Tensorlink),intent(inout) ::h
		type(Tensornode),pointer,intent(inout) ::node
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Tend=>node
		else
			h%length=h%length+1
			h%Tend%next=>node
			h%Tend=>node
		end if
		return
	end subroutine
!	add a Tensornode to the begin of the link	
!	there is no vector index in the  input Tensor
	subroutine push_fonode(h,node)
		type(Tensorlink),intent(inout) ::h
		type(Tensornode),pointer,intent(inout)  ::node
		type(Tensornode),pointer ::p
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Tend=>node
		else
			p=>h%head
			h%head=>node
			node%next=>p
			h%length=h%length+1
		end if
		return
	end subroutine
!	add a compelx munber to the end of the link
!	there is no vector index in the  input Tensor
	subroutine push_backnum(h,num)
		type(Tensorlink),intent(inout) ::h
		complex*16,intent(in) ::num
		type(Tensor)::Ten
		Ten=(/num/)
		call push_backTen(h,Ten)
		return
	end subroutine
!	add a compelx munber to the begin of the link	
!	there is no vector index in the  input Tensor
	subroutine push_fonum(h,num)
		type(Tensorlink),intent(inout) ::h
		complex*16,intent(in) ::num
		type(Tensor)::Ten
		Ten=(/num/)
		call push_fo(h,Ten)
		return
	end subroutine
!	add a compelx munber to the end of the link
!	there is vector index in the  input Tensor
	subroutine push_backnumI(h,num,inde)
		type(Tensorlink),intent(inout) ::h
		complex*16,intent(in) ::num
		integer,intent(in) :: inde(:)
		type(Tensor)::Ten
		Ten=(/num/)
		call push_backI(h,Ten,inde)
		return
	end subroutine
!	add a compelx munber to the begin of the link	
!	there is vector index in the  input Tensor
	subroutine push_fonumI(h,num,inde)
		type(Tensorlink),intent(inout) ::h
		complex*16,intent(in) ::num
		integer,intent(in) :: inde(:)
		type(Tensor)::Ten
		Ten=(/num/)
		call push_foI(h,Ten,inde)
		return
	end subroutine
!	return the inde Tensor in the link
	type(Tensor) function  Ti(h,inde)		
		type(Tensorlink),intent(in) :: h
		integer,intent(in) :: inde
		integer ::i
		type(Tensornode),pointer ::p
		if(h%length.lt.inde) then
			write(*,*)"ERROR in Ti"
		else
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			Ti=p%Ten
		end if
		return
	end  function
!	return the inde Tensor in the link
!	inde is vector	
	type(Tensor) function  TiI(h,inde)		
		type(Tensorlink),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i
		type(Tensornode),pointer ::p
		logical::continu
		continu=.true.
		if(lenOfIndice(h).ne.size(inde)) then
			write(*,*)"ERROR in TiI"
		else
			p=>h%head
			i=1
			do while (continu)
				if(inde.equ.p%indice) then
					TiI=p%Ten
					continu=.false.
				end if
				p=>p%next
				i=i+1
				if(i.eq.h%length+1) then
					continu=.false.
				end if
			end do
		end if
		return
	end  function
!	return the inde Tensornode's address in the link
!	on return,output is a pointer
	 subroutine  node_i_int(h,p,inde)
	 	type(Tensornode),pointer,intent(inout)::p
		type(Tensorlink),intent(in) :: h
		integer,intent(in) :: inde
		integer ::i
		if(h%length.lt.inde) then
			write(*,*)"ERROR in Ti"
		else
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
		end if
		return
	end  subroutine
!	return the inde Tensornode's address in the link
!	inde is vector	
!	on return,output is a pointer
	subroutine  node_i_vec(h,p,inde) 
		type(Tensorlink),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i
		type(Tensornode),pointer,intent(inout) ::p
		logical::continu
		continu=.true.
		if(lenOfIndice(h).ne.size(inde)) then
			write(*,*)"ERROR in TiI"
		else
		
			i=1
			p=>h%head
			if(inde.equ.p%indice) then
				return
			end if
			
			do while (continu)
				i=i+1
				p=>p%next
				if(inde.equ.p%indice) then
					return
				end if
				if(i.gt.h%length+1) then
					write(*,*)"ERROR in node_i_vec"
					stop
				end if
			end do
			
		end if
		return
	end  subroutine	
! return the length of the index of the first Tesnor
	integer function lenOfIndice(h)
		type(Tensorlink),intent(in) :: h
		type(Tensornode),pointer ::p
		p=>h%head
		lenOfIndice=size(p%indice)
		return
	end function
!	return the index of the inde Tensor
!	T_{1,1}->T_{1,2}->T_{1,3}->T_{2,1}->T_{2,2}
!	inde=2 => indice=[1,2]
! on entry the size of indice should be equal to the one in h
	subroutine TenIndice(h,inde,indice)
		type(Tensorlink),intent(in) :: h
		integer,intent(in) :: inde
		integer,intent(out) :: indice(:)
		integer :: i,lenOfind
		type(Tensornode),pointer ::p
		lenOfind=size(indice)
		if(lenOfind.ne.lenOfIndice(h)) then
			write(*,*) "Error in TenIndice"
			write(*,*) lenOfind,lenOfIndice(h)
			p=>h%head
			p=>p%Next
			write(*,*)p%indice
			write(*,*)"stop"
			stop
		end if
		if(h%length.lt.inde) then
			write(*,*)"ERROR in TenIndice"
		else
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			indice=p%indice
		end if
		return
	end subroutine

! if the next node exict,p point to the next node,or p remain nochange
	logical function ifnextnode(p)
		type(Tensornode),pointer,intent(inout) ::p
		if(associated(p%next)) then
			ifnextnode=.true.
			p=>p%next
		else
			ifnextnode=.false.
		end if
	end function
!	p is a pointer of Tensornode,p points to next
	subroutine nextnode(p)
		type(Tensornode),pointer,intent(inout) ::p
		if(associated(p%next)) then
			p=>p%next
		else
			write(*,*)"Error in nextnode,p is pointing to the end of the link"
			stop
		end if
	end subroutine
!	p is a pointer of Tensornode,p points to head of the h	
	subroutine headnode(h,p)
		type(Tensorlink),intent(in)::h
		type(Tensornode),pointer,intent(inout) ::p
		p=>h%head
	end subroutine
!	p is a pointer of Tensornode,p points to end of the h	
	subroutine endnode(h,p)
		type(Tensorlink),intent(in)::h
		type(Tensornode),pointer,intent(inout) ::p
		p=>h%Tend
	end subroutine
!	add a empty node to the link
	subroutine addnode(h)
		type(Tensorlink),intent(inout)::h
		type(Tensornode),pointer ::p
		allocate(p)
		call push_back(h,p)
		return
	end subroutine
		
		
!	get the index of the inde Tensor,and return .true.
!	T_{1,1}->T_{1,2}->T_{1,3}->T_{2,1}->T_{2,2}
!	inde=2 => indice=[1,2]
!	if there no index in h,return .false.
	logical function TenIndiceLog(h,inde,indice)
		type(Tensorlink),intent(in) :: h
		integer,intent(in) :: inde
		integer,allocatable,intent(out) :: indice(:)
		integer :: i
		type(Tensornode),pointer ::p
		if(h%length.lt.inde) then!if[2]
			write(*,*)"ERROR in TenIndiceLog"
		else!if[2]
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			if(allocated(p%indice)) then!if[3]
				TenIndiceLog=.true.
				allocate(indice(size(p%indice)))
				indice=p%indice
			else!if[3]
				TenIndiceLog=.false.
			end if!if[3]
		end if!if[2]
		return
	end function
!check and output length of the link		
!check if the length of the link is equal to the value in head of the link
	integer function Checklength(link)
		type(Tensorlink),intent(in) :: link
		type(Tensornode),pointer ::p
		integer ::num
		Checklength=link%length
		p=>link%head
		num=0
		do while(associated(p))
			p=>p%next
			num=num+1
		end do
		if(num.ne.Checklength)then
			write(*,*)"The length of link is",num
			write(*,*)"The value of length in link is" ,Checklength
			write(*,*)"Error in length of link"
			write(*,*)"stop !"
			stop
		end if
		return
	end function
!	output length of the link		
	integer function linklength(link)
		type(Tensorlink),intent(in) :: link
		linklength=link%length
		return
	end function	
	
!	modify the inde element of the link
!	The Tensor will be modify
	subroutine modifylink(h,T,inde)
		type(Tensorlink),intent(inout):: h
		integer,intent(in) :: inde
		type(Tensor),intent(in):: T
		integer ::i
		type(Tensornode),pointer ::p
		if(h%length.lt.inde) then
			write(*,*)"ERROR in modify,length"
			else
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			p%Ten=T
		end if
		return
	end subroutine
!	clean the link
	subroutine cleanlink(h)
		type(Tensorlink),intent(inout) :: h
		type(Tensornode),pointer ::p,p1
		if(h%length.eq.0) then!if[2]
			nullify(h%head)
			nullify(h%Tend)
			return
		end if!if[2]
		p=>h%head
		do while(associated(p))
			p1=>p%next
			call cleanTensor(p%Ten)
			deallocate(p)
			p=>p1
		end do
		h%length=0
		nullify(h%head)
		nullify(h%Tend)
		return
	end subroutine
!	copy the link hin to hout
!	hin will not change
	subroutine copylink(hout,hin)
		type(Tensorlink),intent(inout) :: hout
		type(Tensorlink),intent(in) :: hin
		integer :: i
		integer,allocatable::indes(:)
		logical :: logi
		call cleanlink(hout)
		do i=1,hin%length
			logi=TenIndiceLog(hin,i,indes)
			if(logi) then!if[2]
				call push_backI(hout,hin.i.i,indes)
			else
				call push_back(hout,hin.i.i)
			end if!if[2]
		end do
		return
	end subroutine
	subroutine copylinkhead(hout,hin)
		type(Tensorlink),intent(inout) :: hout
		type(Tensorlink),intent(in) :: hin
		integer :: i
		integer,allocatable::indes(:)
		logical :: logi
		call cleanlink(hout)
		hout%head=>hin%head
		hout%Tend=>hin%Tend
		hout%length=hin%length
		return
	end subroutine
!	connect two links:
!			link1 :T1->T2->T4->T5...->TN
!			link2: TT1->TT2->..->TTM
!			result link will b2 T1->T2->..->TN->TT1->TT2->...->TTM
	type(Tensorlink) function connectlink(link1,link2)
		type(Tensorlink),intent(in) :: link1,link2
		type(Tensorlink) :: Templink
		type(Tensornode),pointer ::p
		connectlink=link1
		Templink=link2
		p=>connectlink%Tend
		p%next=>Templink%head
		connectlink%Tend=>Templink%Tend
		connectlink%length=connectlink%length+Templink%length
		return
	end function
!***************  deletelink  **********************
!	link   :	T1->T2->...->T_{inde1-1}->T_inde1->T_{inde1+1}->....->T_{inde2-1}->T_inde2->T_{inde2+1}...->TN 
!  result:		T1->T2->...->T_{inde1-1}->T_{inde1+1}....->T_{inde2-1}->->T_{inde2+1}...->TN 
!  delete the inde1 to inde2 Tensor in the link,note that the deleted Tensor are include the Tensors of inde1 and inde2
	type(Tensorlink) function deletelink(link_in,inde1,inde2) result(link)			
		type(Tensorlink),intent(in) :: link_in
		integer,intent(in)::inde1,inde2
		type(Tensorlink) :: Templink
		type(Tensornode),pointer ::p1,p2,p
		integer::i
		if(inde1.eq.link%length) then!if[2]
			write(*,*) "error"
			write(*,*) "inde1 larger than the length of link"
			write(*,*)"program will stop"
			stop
		end if!if[2]
		if(inde1.gt.inde2) then!if[2]
			write(*,*)"error"
			write(*,*) "inde1 should not larger than inde2"
			write(*,*)"program will stop"
			stop
		end if!if[3]
		link=link_in
		if(inde1.eq.1) then!if[4]
			if(inde2.eq.link%length) then!if[5]
				call cleanlink(link)
				return
			end if!if[5]
			p2=>link%head
			do i=1,inde2
				p=>p2
				p2=>p2%next
				call cleanTensor(p%Ten)
				deallocate(p)
			end do
			link%head=>p2
			link%length=link%length-(inde2-inde1)-1
			return
		end if!if[4]
		if(inde2.eq.link%length) then!if[6]
			p1=>link%head
			do i=2,inde1-1
				p1=>p1%next
			end do
			p2=>p1%next
			link%Tend=>p1
			p=>p1%next
			nullify(p1)
			do while(associated(p))
				p2=>p%next
				call cleanTensor(p%Ten)
				deallocate(p)
				p=>p2
			end do
			link%length=link%length-(inde2-inde1)-1
			return
		end if!if[6]
		p1=>link%head
		do i=2,inde1-1
			p1=>p1%next
		end do
		p2=>p1%next
		do i=inde1,inde2
			p=>p2
			p2=>p2%next
			call cleanTensor(p%Ten)
			deallocate(p)
		end do
		p1%next=>p2
		link%length=link%length-(inde2-inde1)-1
		return
	end function		

!Tensor=Tensorlink,if every element(Tensor) in the link is a one-element Tensor
	subroutine TenLink(Ten,link)
		type(Tensor),intent(inout)::Ten
		type(Tensorlink),intent(in)::link
		complex*16,allocatable::Tensordata(:)
		integer::Tenlen,i
		logical::goon
		type(Tensornode),pointer::p
		Tenlen=linklength(link)
		allocate(Tensordata(Tenlen))
		call headnode(link,p)
		do i=1,Tenlen
			if(getTotaldata(p%Ten).eq.1)then
				Tensordata(i)=p%Ten.i.1
			else
				write(*,*)"error in assignment of a link to Tensor"
			end if
			goon=ifnextnode(p)
		end do
		Ten=Tensordata
		nullify(p)
		return
	end subroutine
	
	subroutine Lprint(link)
		type(Tensorlink),intent(in) :: link
		integer::lenlink,i
		logical::goon,goon2
		integer,allocatable::indice(:)
		lenlink=Checklength(link)
		if(lenlink.eq.0) then
			write(*,*)"There is no Tensor in the link"
			return
		end if
		write(*,*) "length of the link is",lenlink
		goon=TenIndiceLog(link,1,indice)
		if(goon) then
			write(*,*)"indice of every Tensor in the link are"
			do i=1,lenlink
				goon2=TenIndiceLog(link,i,indice)
				if(goon2) then
					write(*,*)indice
				else
					write(*,*)"no indece"
				end if
			end do
		else
			goon=.true.
			do i=1,lenlink
				goon2=TenIndiceLog(link,i,indice)
				if(goon2) then
					write(*,*)"error in the link,there are some Tensor with indice,some are not"
					goon=.false.
				end if
			end do
		end if
		return
	end subroutine
	
	subroutine LDprint(link)
		type(Tensorlink),intent(in) :: link
		type(Tensornode),pointer::p
		integer::lenlink,i
		logical::goon
		integer,allocatable::indice(:)
		lenlink=Checklength(link)
		if(lenlink.eq.0) then
			write(*,*)"There is no Tensor in the link"
			return
		end if
		p=>link%head
		call Dprint0(p%Ten%TenDim)
		if(allocated(p%indice)) then
			write(*,*)"indice",p%indice
		end if
		i=1
		do while(associated(p%next))
			p=>p%next
			write(*,*)"----------------------------------",i
			i=i+1
			call Dprint0(p%Ten%TenDim)
			if(allocated(p%indice)) then
				write(*,*)"indice",p%indice
			end if
		end do
		write(*,*)"--------- End of the link -----------",i
		return
	end subroutine





!****************************************************************************************
	 
! 	Tensorlist and listnode 	  





!	add a Tensorlink to the end of the list
!	there is no vector index in the  input Tensorlink
	subroutine push_backlink(h,L)
		type(Tensorlist),intent(inout) ::h
		type(Tensorlink),intent(in):: L
		type(listnode),pointer ::node
		allocate(node)
		node%link=L
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Lend=>node
		else
			h%length=h%length+1
			h%Lend%next=>node
			h%Lend=>node
		end if
		return
	end subroutine
!	add a Tensorlink to the begin of the list	
!	there is no vector index in the  input Tensorlink
	subroutine push_folink(h,L)
		type(Tensorlist),intent(inout) ::h
		type(Tensorlink),intent(in):: L
		type(listnode),pointer ::node,p
		allocate(node)
		node%link=L
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Lend=>node
		else
			p=>h%head
			h%head=>node
			node%next=>p
			h%length=h%length+1
		end if
		return
	end subroutine
!	add a Tensorlink to the end of the list	
!	inde is the index of the Tensorlink
	subroutine push_backlinkI(h,L,inde)
		type(Tensorlist),intent(inout) ::h
		type(Tensorlink),intent(in):: L
		integer,intent(in) :: inde(:)
		type(listnode),pointer ::node
		allocate(node)
		node%link=L
		call assignment_int_dim1(node%indice,inde)
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Lend=>node
		else
			h%length=h%length+1
			h%Lend%next=>node
			h%Lend=>node
		end if
		return
	end subroutine
!	add a Tensorlink to the begin of the list	
!	inde is the index of the Tensorlink	
	subroutine push_folinkI(h,L,inde)
		type(Tensorlist),intent(inout) ::h
		type(Tensorlink),intent(in):: L
		integer,intent(in) :: inde(:)
		type(listnode),pointer ::node,p
		allocate(node)
		node%link=L
		call assignment_int_dim1(node%indice,inde)
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Lend=>node
		else
			p=>h%head
			h%head=>node
			node%next=>p
			h%length=h%length+1
		end if
		return
	end subroutine
!	add a listnode to the end of the list
!	there is no vector index in the  input Tensorlink
	subroutine push_backlinknode(h,node)
		type(Tensorlist),intent(inout) ::h
		type(listnode),pointer,intent(inout) ::node
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Lend=>node
		else
			h%length=h%length+1
			h%Lend%next=>node
			h%Lend=>node
		end if
		return
	end subroutine
!	add a listnode to the end of the list
!	there is no vector index in the  input Tensor
	subroutine push_folinknode(h,node)
		type(Tensorlist),intent(inout) ::h
		type(listnode),pointer,intent(inout)  ::node
		type(listnode),pointer ::p
		nullify(node%next)
		if(h%length.eq.0) then
			h%head=>node
			h%length=1
			h%Lend=>node
		else
			p=>h%head
			h%head=>node
			node%next=>p
			h%length=h%length+1
		end if
		return
	end subroutine
!	return the inde Tensorlink in the list
	subroutine  Li(link,h,inde)	
		type(Tensorlink),intent(inout) :: link	
		type(Tensorlist),intent(in) :: h
		integer,intent(in) :: inde
		integer ::i
		type(listnode),pointer ::p
		if(h%length.lt.inde) then
			write(*,*)"ERROR in Ti"
		else
			call cleanlink(link)
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			link=p%link
		end if 
		return
	end  subroutine
!	return the inde Tensorlink in the list
!	inde is vector	
	subroutine  LiI(link,h,inde)		
		type(Tensorlink),intent(inout) :: link
		type(Tensorlist),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i
		type(listnode),pointer ::p
		logical::continu
		continu=.true.
		call cleanlink(link)
		p=>h%head
		i=1
		do while (continu)
			if(inde.equ.p%indice) then
				link=p%link
				continu=.false.
			end if
			p=>p%next
			i=i+1
			if(i.eq.h%length+1) then
				continu=.false.
			end if
		end do
		return
	end  subroutine
!	return the inde Tensorlink in the list
	type(Tensor) function  Liv(h,inde)	
		type(Tensorlist),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i,j
		type(listnode),pointer ::p
		type(Tensornode),pointer ::Tp
		if(h%length.lt.inde(1)) then
			write(*,*)"ERROR in Ti"
		else
			p=>h%head
			do i=1,inde(1)-1
				p=>p%next
			end do
			Tp=>p%link%head
			do j=1,inde(2)-1
				Tp=>Tp%next
			end do
			Liv=Tp%Ten
		end if 
		return
	end  function
!	return the inde Tensorlink in the list
!	inde is vector	
	type(Tensor) function  LiIv(h,inde)		
		type(Tensorlist),intent(in) :: h
		integer,intent(in) :: inde(:,:)
		integer ::i,j,size1,size2
		type(listnode),pointer ::p
		type(Tensornode),pointer ::Tp
		logical::continu
		integer,allocatable::indetemp(:)
		continu=.true.
		p=>h%head
		i=1
!search link		
		if(inde(1,1).eq.0)then!if[1],no inde
			p=>h%head
			do i=1,inde(1,2)-1
				p=>p%next
			end do
			Tp=>p%link%head
		else!if[1]
			p=>h%head
			i=1
			size1=inde(1,1)+1
			if(allocated(indetemp))then
				deallocate(indetemp)
			end if
			allocate(indetemp(inde(1,1)))
			indetemp=inde(1,2:size1)
			do while (continu)
				if(indetemp.equ.p%indice) then
					Tp=>p%link%head
					continu=.false.
				end if
				p=>p%next
				i=i+1
				if(i.eq.h%length+1) then
					continu=.false.
				end if
			end do
		end if!if[1]
!search Tensor		
		if(inde(2,1).eq.0)then!if[2]
				do j=1,inde(2,2)-1
					Tp=>Tp%next
				end do
				LiIv=Tp%Ten
				return
			else!if[2]
				j=1
				continu=.true.
				size2=inde(2,1)+1
				if(allocated(indetemp))then
				deallocate(indetemp)
				end if
				allocate(indetemp(inde(2,1)))
				indetemp=inde(2,2:size2)
				do while (continu)
					if(indetemp.equ.Tp%indice) then!if[3]
						LiIv=Tp%Ten
						continu=.false.
						return
					end if!if[3]
					Tp=>Tp%next
					j=j+1
					if(j.eq.p%link%length+1) then!if[3]
						continu=.false.
					end if!if[4]
				end do
				return
			end if!if[2]
	end  function
!	return the inde listnode's address in the list
!	on return,output is a pointer
	 subroutine  listnode_i_int(h,p,inde)
	 	type(listnode),pointer,intent(inout)::p
		type(Tensorlist),intent(in) :: h
		integer,intent(in) :: inde
		integer ::i
		if(h%length.lt.inde) then
			write(*,*)"ERROR in listnode_i_int"
		else
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
		end if
		return
	end  subroutine
!	return the inde listnode's address in the link
!	inde is vector	
!	on return,output is a pointer
	subroutine  listnode_i_vec(h,p,inde) 
		type(Tensorlist),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i
		type(listnode),pointer,intent(inout) ::p
		logical::continu
		continu=.true.
		if(lenOfIndicelist(h).ne.size(inde)) then
			write(*,*)"ERROR in listnode_i_vec"
		else
		
			i=1
			p=>h%head
			if(inde.equ.p%indice) then
				return
			end if
			
			do while (continu)
				i=i+1
				p=>p%next
				if(inde.equ.p%indice) then
					return
				end if
				if(i.gt.h%length+1) then
					write(*,*)"ERROR in listnode_i_vec"
					stop
				end if
			end do
			
		end if
		return
	end  subroutine	
! return the length of the index of the first Tesnorlink
	integer function lenOfIndicelist(h)
		type(Tensorlist),intent(in) :: h
		type(listnode),pointer ::p
		p=>h%head
		lenOfIndicelist=size(p%indice)
		return
	end function
!	return the index of the inde Tensorlink
!	T_{1,1}->T_{1,2}->T_{1,3}->T_{2,1}->T_{2,2}
!	inde=2 => indice=[1,2]
! on entry the size of indice should be equal to the one in h
	subroutine linkIndice(h,inde,indice)
		type(Tensorlist),intent(in) :: h
		integer,intent(in) :: inde
		integer,intent(out) :: indice(:)
		integer :: i,lenOfind
		type(listnode),pointer ::p
		lenOfind=size(indice)
		if(lenOfind.ne.lenOfIndicelist(h)) then
			write(*,*) "Error in linkIndice"
			write(*,*) lenOfind,lenOfIndicelist(h)
			p=>h%head
			p=>p%Next
			write(*,*)p%indice
			write(*,*)"stop"
			stop
		end if
		if(h%length.lt.inde) then
			write(*,*)"ERROR in linkIndice"
		else
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			indice=p%indice
		end if
		return
	end subroutine

! if the next node exict,p point to the next node,or p remain nochange
	logical function ifnextlistnode(p)
		type(listnode),pointer,intent(inout) ::p
		if(associated(p%next)) then
			ifnextlistnode=.true.
			p=>p%next
		else
			ifnextlistnode=.false.
		end if
	end function
!	p is a pointer of listnode,p points to next
	subroutine nextlistnode(p)
		type(listnode),pointer,intent(inout) ::p
		if(associated(p%next)) then
			p=>p%next
		else
			write(*,*)"Error in next;istnode,p is pointing to the end of the list"
			stop
		end if
	end subroutine
!	p is a pointer of listnode,p points to head of the h	
	subroutine headlist(h,p)
		type(Tensorlist),intent(in)::h
		type(listnode),pointer,intent(inout) ::p
		p=>h%head
		return
	end subroutine
!	p is a pointer of listnode,p points to end of the h	
	subroutine endlist(h,p)
		type(Tensorlist),intent(in)::h
		type(listnode),pointer,intent(inout) ::p
		p=>h%Lend
	end subroutine
!	add a empty listnode to the list
	subroutine addlist(h)
		type(Tensorlist),intent(inout)::h
		type(listnode),pointer ::p
		allocate(p)
		call push_backlinknode(h,p)
		return
	end subroutine
		
		
!	get the index of the inde Tensorlink,and return .true.
!	T_{1,1}->T_{1,2}->T_{1,3}->T_{2,1}->T_{2,2}
!	inde=2 => indice=[1,2]
!	if there no index in h,return .false.
	logical function linkIndiceLog(h,inde,indice)
		type(Tensorlist),intent(in) :: h
		integer,intent(in) :: inde
		integer,allocatable,intent(out) :: indice(:)
		integer :: i
		type(listnode),pointer ::p
		if(h%length.lt.inde) then!if[2]
			write(*,*)"ERROR in TenIndiceLog"
		else!if[2]
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			if(allocated(p%indice)) then!if[3]
				linkIndiceLog=.true.
				allocate(indice(size(p%indice)))
				indice=p%indice
			else!if[3]
				linkIndiceLog=.false.
			end if!if[3]
		end if!if[2]
		return
	end function
!check and output length of the link		
!check if the length of the link is equal to the value in head of the link
	integer function Checklistlength(list)
		type(Tensorlist),intent(in) :: list
		type(listnode),pointer ::p
		integer ::num
		Checklistlength=list%length
		p=>list%head
		num=0
		do while(associated(p))
			p=>p%next
			num=num+1
		end do
		if(num.ne.Checklistlength)then
			write(*,*)"The length of link is",num
			write(*,*)"The value of length in link is" ,Checklistlength
			write(*,*)"Error in length of link"
			write(*,*)"stop !"
			stop
		end if
		return
	end function
!	output length of the link		
	integer function listlength(list)
		type(Tensorlist),intent(in) :: list
		listlength=list%length
		return
	end function	
	
!	modify the inde element of the list
!	The Tensorlink will be modify
	subroutine modifylist(h,L,inde)
		type(Tensorlist),intent(inout):: h
		integer,intent(in) :: inde
		type(Tensorlink),intent(in):: L
		integer ::i
		type(listnode),pointer ::p
		if(h%length.lt.inde) then
			write(*,*)"ERROR in modifylink,length"
		else
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			p%link=L
		end if
		return
	end subroutine
!	clean the list
	subroutine cleanlist(h)
		type(Tensorlist),intent(inout) :: h
		type(listnode),pointer ::p,p1
		if(h%length.eq.0) then!if[2]
			nullify(h%head)
			nullify(h%Lend)
			return
		end if!if[2]
		p=>h%head
		do while(associated(p))
			p1=>p%next
			call cleanlink(p%link)
			deallocate(p)
			p=>p1
		end do
		h%length=0
		nullify(h%head)
		nullify(h%Lend)
		return
	end subroutine
!	copy the list hin to hout
!	hin will not change
	subroutine copylist(hout,hin)
		type(Tensorlist),intent(inout) :: hout
		type(Tensorlist),intent(in) :: hin
		type(Tensorlink)::tempL
		integer :: i
		integer,allocatable::indes(:)
		logical :: logi
		call cleanlist(hout)
		do i=1,hin%length
			logi=linkIndiceLog(hin,i,indes)
			call link_i(tempL,hin,indes)
			if(logi) then!if[2]
				call push_backlinkI(hout,tempL,indes)
			else
				
				call push_backlink(hout,tempL)
			end if!if[2]
		end do
		call cleanlink(tempL)
		return
	end subroutine
!	connect two lists:
!			list1 :T1->T2->T4->T5...->TN
!			list2: TT1->TT2->..->TTM
!			result list will b2 T1->T2->..->TN->TT1->TT2->...->TTM
	type(Tensorlist) function connectlist(list1,list2)
		type(Tensorlist),intent(in) :: list1,list2
		type(Tensorlist) :: Templist
		type(listnode),pointer ::p
		connectlist=list1
		Templist=list2
		p=>connectlist%Lend
		p%next=>Templist%head
		connectlist%Lend=>Templist%Lend
		connectlist%length=connectlist%length+Templist%length
		return
	end function
!***************  deletelink  **********************
!	list   :	T1->T2->...->T_{inde1-1}->T_inde1->T_{inde1+1}->....->T_{inde2-1}->T_inde2->T_{inde2+1}...->TN 
!  result:		T1->T2->...->T_{inde1-1}->T_{inde1+1}....->T_{inde2-1}->->T_{inde2+1}...->TN 
!  delete the inde1 to inde2 Tensor in the link,note that the deleted Tensor are include the Tensors of inde1 and inde2
	type(Tensorlist) function deletelist(list_in,inde1,inde2) result(list)			
		type(Tensorlist),intent(in) :: list_in
		integer,intent(in)::inde1,inde2
		type(Tensorlist) :: Templink
		type(listnode),pointer ::p1,p2,p
		integer::i
		if(inde1.gt.list_in%length) then!if[2]
			write(*,*) "error"
			write(*,*) "inde1 larger than the length of list"
			write(*,*)"program will stop"
			stop
		end if!if[2]
		if(inde1.gt.inde2) then!if[2]
			write(*,*)"error"
			write(*,*) "inde1 should not larger than inde2"
			write(*,*)"program will stop"
			stop
		end if!if[3]
		list=list_in
		if(inde1.eq.1) then!if[4]
			if(inde2.eq.list%length) then!if[5]
				call cleanlist(list)
				return
			end if!if[5]
			p2=>list%head
			do i=1,inde2
				p=>p2
				p2=>p2%next
				call cleanlink(p%link)
				deallocate(p)
			end do
			list%head=>p2
			list%length=list%length-(inde2-inde1)-1
			return
		end if!if[4]
		if(inde2.eq.list%length) then!if[6]
			p1=>list%head
			do i=2,inde1-1
				p1=>p1%next
			end do
			p2=>p1%next
			list%Lend=>p1
			p=>p1%next
			nullify(p1)
			do while(associated(p))
				p2=>p%next
				call cleanlink(p%link)
				deallocate(p)
				p=>p2
			end do
			list%length=list%length-(inde2-inde1)-1
			return
		end if!if[6]
		p1=>list%head
		do i=2,inde1-1
			p1=>p1%next
		end do
		p2=>p1%next
		do i=inde1,inde2
			p=>p2
			p2=>p2%next
			call cleanlink(p%link)
			deallocate(p)
		end do
		p1%next=>p2
		list%length=list%length-(inde2-inde1)-1
		return
	end function		
	
	subroutine Listprint(list)
		type(Tensorlist),intent(in) :: list
		integer::lenlist,i
		logical::goon,goon2
		integer,allocatable::indice(:)
		lenlist=Checklistlength(list)
		if(lenlist.eq.0) then
			write(*,*)"There is no Tensorlink in the list"
			return
		end if
		write(*,*) "length of the list is",lenlist
		goon=linkIndiceLog(list,1,indice)
		if(goon) then
			write(*,*)"indice of every Tensorlink in the list are"
			do i=1,lenlist
				goon2=linkIndiceLog(list,i,indice)
				if(goon2) then
					write(*,*)indice
				else
					write(*,*)"no indece"
				end if
			end do
		else
			goon=.true.
			do i=1,lenlist
				goon2=linkIndiceLog(list,i,indice)
				if(goon2) then
					write(*,*)"error in the list,there are some Tensorlink with indice,some are not"
					goon=.false.
				end if
			end do
		end if
		return
	end subroutine
	
	subroutine Listlinkprint(list)
		type(Tensorlist),intent(in) :: list
		type(listnode),pointer::p
		integer::lenlist,i
		logical::goon
		integer,allocatable::indice(:)
		lenlist=Checklistlength(list)
		if(lenlist.eq.0) then
			write(*,*)"There is no Tensor in the link"
			return
		end if
		p=>list%head
		call LDprint(p%link)
		if(allocated(p%indice)) then
			write(*,*)"link indice",p%indice
		end if
		i=1
		do while(associated(p%next))
			p=>p%next
			write(*,*)"======================================"
			write(*,*)"----------------------------------",i
			write(*,*)"======================================"
			i=i+1
			call LDprint(p%link)
			if(allocated(p%indice)) then
				write(*,*)"link indice",p%indice
			end if
		end do
		write(*,*)"*********** End of the list *************",i
		return
	end subroutine
	



!******************  init_random_seed *****************
	subroutine init_random_seed()
		integer :: i,n,clock
		integer,dimension(:),allocatable :: seed
			call random_seed(size=n)
			allocate(seed(n))
			call system_clock(count=clock)
			seed=clock+37*(/(i-1,i=1,n)/)
			call random_seed(PUT=seed) 
			deallocate(seed)
	end subroutine
	

!*****************  warning   *****************
!			if the largest data in T is too large,
!			warning	
	subroutine warning(T,chara)
		character*50,intent(in) :: chara
		type(Tensor),intent(in) :: T
		real*8	::maxV
		maxV=maxelement(T)
		if(maxV.ge.large_number_warning) then
			write (*,*) "***   warning   ***"
			write (*,*) chara
			write (*,*) "data in the Tensor is large,the abs value largest element is"
			write(*,*)		maxV
		end if
	end subroutine
	subroutine warning2(T)
		type(Tensor),intent(in) :: T
		real*8	::maxV
		maxV=maxelement(T)
		if(maxV.ge.large_number_warning) then
			write (*,*) "***   warning   ***"
			write (*,*) "data in the Tensor is large,the abs value largest element is"
			write(*,*)		maxV
		end if
	end subroutine	
!*****************  NANwarning   *****************
!			warning		NAN
	subroutine NANwarning(T,chara)
		character*50,intent(in) :: chara
		type(Tensor),intent(in) :: T
		complex*16::val
		integer ::i,total
		total=T%TotalData
		do i=1,total
			val=T.i. i
			if(isnan(real(val)).or.isnan(aimag(val))) then
				write (*,*) "***   NANwarning   ***"
				write (*,*) chara
				return
			end if
		end do
		return
	end subroutine
	subroutine NANwarning2(T)
		type(Tensor),intent(in) :: T
		complex*16::val
		integer ::i,total
		total=T%TotalData
		do i=1,total
			val=T.i. i
			if(isnan(real(val)).or.isnan(aimag(val))) then
				write (*,*) "***   NANwarning   ***"
				write (*,*) "stop"
				stop
			end if
		end do
		return
	end subroutine	
!*****************  NANjudge   *****************
!			warning	NAN
	logical function NANjudge(T)
		type(Tensor),intent(in) :: T
		integer ::i,total
		complex*16::val
		NANjudge=.false.
		total=T%TotalData
		do i=1,total
			val=T.i. i
			if(isnan(real(val)).or.isnan(aimag(val))) then
				NANjudge=.true.
			return
			end if
		end do
		return
	end function
	logical function NANjudgev(val)
		complex*16,intent(in) :: val(:)
		integer ::i,total
		NANjudgev=.false.
		total=size(val)
		do i=1,total
			if(isnan(real(val(i))).or.isnan(aimag(val(i)))) then
				NANjudgev=.true.
			return
			end if
		end do
		return
	end function			




end module
!***************************************************
!***************** END OF TNESOR *************
!***************************************************


			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
