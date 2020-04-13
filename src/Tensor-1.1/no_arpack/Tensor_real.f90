
!		 This file is modify from Tensor.f90 which is used for complex*16 type
!	no test on this code

!************************************************************
!************* START OF DTensor **********************
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
!			3.A+B=C,where A is real*8,B is real*8,and C is real*8.It goes
!		no wrong in fortran90.This code is in Dadd_DTen1 and Dadd_DTen2
!			4.when link1=link2,and want link1 unchange while operate on link2,use
!		DcopyLink but not (=)
!			5.The code is fortran90 version.
!*****************     note      *******************
!			1.The product of DTensor only allow rank<=2,if larger than 2,
!		use .p. and Dcontract(or .c.) to reshape it.
!			2.Dpermutation(.p.)  for rank<=3 is faster then the one with input a vector
!		,if larger then 3,use Dcontract(or .c.) to reshape it to rank=3 or rank=2,and
!		then do the Dpermutation,and at last Ddecompose the DTensor back to the original
!		rank.One could reshape to any dimension through (.p.vectot),but it is slow
!		than (.p.number).
!			3.to compile the code one should link the files to lapack and blas.
!			4.svd ,Deng and expm could be more faster by define the input matrix of the 
!		lapack subroutine as a vector of the output DTensor,becase the data is already
!		store as a vector,when copy to a matrix,it waste time.This save the time of 
!		assignment of the DTensor
!			5.These files call lapack and blas .mpif90 name.f90 -o name -llapack -lrefblas 
!			6.Send email to tyrants@qq.com to report any bugs.
!***********************************************************
!		bugs do not fix yet:
!			Dpermutation:most sitiation,it gose right
!			expm:if the input matrix is not hermitian,gose wrong			
!
!		Dmodify on 2014.4.13
!			difine DcopyLinkhead as (=) ,not DcopyLink . So when using link or list,
!		test the function if it will create rubbish . The test include some like 
!		link=function()

!  could be more faster after Dmodify:DcopyLink

include "/lib/Tensor/MPI/Dimension.f90"
module Tensor_real
	use Dimension_typede
	implicit none
!****************************************************	
!*********		define data type	***************
	type  DTensor
		integer,private :: rank!length of Dimension.TenDim%Dimsize
		type(Dimension),private	:: TenDim!dimenison
		real*8,allocatable,private :: DTensor_data(:)!The data store in one dimenson array
		integer,private :: totalData!length of DTensor_data
		logical*4,private :: flag=.false. !if flag=true ,it means there are data in DTensor
	end type DTensor
! 	DTensorlink and DTensornode used as vector
!	DTensorlink store the head and the last DElement of DTensornode
!	DTensornode is used to store DTensor
!	indice(:) is used as index such as H_{1,2},H_{1,2,3} and so on
	type DTensornode
		integer,allocatable:: indice(:)
		type(DTensor) ::Ten
		type(DTensornode),private,pointer :: next
	end type DTensornode
		
	type DTensorlink
		integer,private :: length=0
		type(DTensornode),pointer,private  :: head
		type(DTensornode),pointer,private  :: Tend
	end type DTensorlink	
!	-----------
	type DTensorlist
		integer,private :: length=0
		type(Dlistnode),private,pointer :: head
		type(Dlistnode),pointer,private  :: Lend
	end type DTensorlist
	
	type Dlistnode
		integer,allocatable:: indice(:)
		type(DTensorlink) ::link
		type(Dlistnode),private,pointer :: next
	end type Dlistnode
	
!*********		End define		*********************


!****************************************************
!*********		define private data	***************	
	real*8,private :: cputime(2)
	real*8,parameter,private :: II=(0,1),EE=(1,0)! II**2=-1
	real*8,parameter,private :: large_number_warning=1d80!use for checking
	real*8,parameter,private :: zero_num=1d-20!
	logical,private,save::check_flag=.false.!use for checking
	real*8,private::Dnrm2!Dnorm :sqrt(A**H * A)
	real*8,private::ddot
	integer,save,private::max_ncv=500!use in funciton of Deigen2 and Deigen_max_min2 ncv<=max_ncv
	External::Dnrm2,ddot!blas's function
	private::DstoreTenData
	private::addressToIndes
	private::addressToIndes2
	private::IndesToaddress
	private::Dcompose_decompse_subroutine1,Dcompose_decompse_subroutine2
	private::Dset,init_random_seed
!*********		End define	***************		
		
	interface DTMprint
		module procedure DTMprint1
		module procedure DTMprint2
	end interface
	interface Tprint
		module procedure DTprint1
		module procedure DTprint2
	end interface
	
	!interface Dsqrt
	!	module procedure DsqrtTen
	!end interface
	
	!interface Deig!eigenDvalue
	!	module procedure Deng
	!	module procedure Deng2
	!end interface
	
	
	interface Dset!create a DTensor,giving dimension and DTensor_Data
					  !DbuildTen can do all the thing of Dset so I define it as private
		module procedure Dset1 
		module procedure Dset2 
	end interface
	
!	copy the array to the T%DTensor_Data,array can be 1,2 or 3 dimension
	interface DstoreTenData
		module procedure DstoreTenData_dim1
		module procedure DstoreTenData_dim2
		module procedure DstoreTenData_dim3
	end interface
	
	interface	DbuildTen!build a DTensor
! 	note when input a vector in to the DTensor of rank=2
! for example DbuildTen((/3,2/),(/a,b,c,d,e,f/)),then the DTensor return is
!	/ a d \
!	| b e |
!	\ c f /
		module procedure DbuildTen1
		module procedure DbuildTen2
		module procedure DbuildTen3
		module procedure DbuildTen4
		module procedure DbuildTen5
		module procedure DbuildTen6
		module procedure DbuildTen8
	end interface
	
	interface Dmodify
		module procedure Dmodifylink
		module procedure Dmodifylist
		module procedure DmodifyTen_val
		module procedure DmodifyTen_dat
		module procedure DmodifyTen_dat2
		module procedure DmodifyTen_dat3
	end interface
!************************************************************	
!operateDTensor(A,B,oA,oB,productflag):
!	oA,oB is a matrix or vector,the first DElement of every row specify
!		1:Dcontracts
!		2:Ddecompose
!		3:DdecomposeAll
!	other DElement of every row is input parameter of the function
!		at last A*B .when flag=1,result store in A with B no change when productflag=2,
!	 the result will store in B
!		this subroutine run for big inoutA and B,do as less as possible operation on the
!	 DTensor_Data of inoutA and B.
!************************************************************	
	interface operateDTensor
		module procedure operateDTensor1
		module procedure operateDTensor2
	end interface
!DdimOperation(DTensor,Operation):
!	 DdimOperations do the Dcontract or Ddecompose on the DTensor with
!	  no change on the Tensot_data in the DTensor.
!	Operation is a matrix or vector,the first DElement of every row specify
!		1:Dcontracts
!		2:Ddecompose
!		3:DdecomposeAll
!	other DElement of every row is input parameter of the function
!		this subroutine run for big DTensor,do as less as possible operation on the
!	 DTensor_Data of DTensor.
	interface DdimOperation
		module procedure DdimOperations
		module procedure DdimOperation1
	end interface
	
	
	interface Dresetdim
		module procedure Dresetdim1
		module procedure Dresetdim2
	end interface
!************************************************************		
	
	interface Deye
!  return a matrix of diag[s1,s2,s3,...s_{max},0,0...]	,a m*n matrix
		module procedure Deye_real!input s(:) is real*8
		module procedure Done
		module procedure Done_Ten!return a identity DTensor 
		module procedure Deye_Ten1!input DTensor,output diag[T.i.1 , T.i.2 , ... , , T.i.n],rank of T should be 1
		module procedure Deye_Ten2!input DTensor,output a matrix of m*n: diag[T.i.1 , T.i.2 , ... , , T.i.n],rank of T should be 1
	end interface
!***********   interface for DTensorlink  and Tesorlist *************************	
	interface Dpush_back
		module procedure Dpush_backTen
		module procedure Dpush_backnode
		module procedure Dpush_backI
		module procedure Dpush_backnum
		module procedure Dpush_backnumI
		
		module procedure Dpush_backlink
		module procedure Dpush_backlinknode
		module procedure Dpush_backlinkI
!(h,L,pointtoL) if pointtoL=true,than the new node will point to L
!  otherwhile ,create a new link	and that add to the list
		
	end interface
	
	interface Dpush_forward
		module procedure Dpush_fo
		module procedure Dpush_fonode
		module procedure Dpush_foI
		module procedure Dpush_fonum
		module procedure Dpush_fonumI
		
		module procedure Dpush_folink
		module procedure Dpush_folinknode
		module procedure Dpush_folinkI
	end interface
	
	interface Dnode_i!(h,p,inde),return p,points to inde DTensor in the link or inde link in the list
		module procedure Dnode_i_int
		module procedure Dnode_i_vec
		
		module procedure DlistDnode_i_int
		module procedure DlistDnode_i_vec
	end interface
	!if define  DcopyLinkhead as (=),it gose no wrong,but Link and list.i.1 are the same data
	!Dlink_i will copy the data to aother link,the output is a new link
	interface Dlink_i!Do not use Link=list.i.1,because the nodes in Link are pointer,doing so will create rubbish data in memory
		module procedure DcopyLi!output the i DTensorlink in the DTensorlist
		module procedure DcopyLiI!output the index DTensorlink in the DTensorlist,index is a vector
	end interface
!*******************   Dpermutation   ****************************	
	interface operator(.p.)
		module procedure Dpermute_rank2 !Dpermute the DTensor whose rank is 2
		module procedure Dpermute_rank3 !Dpermute the DTensor whose rank is 3
		module procedure Dpermutation	 !Dpermute the DTensor of any rank,give the new order of the dimension
												 !		If operate on a big DTensor,use Dpermute_rank2 or Dpermute_rank3,
												 !	 they are faster.if rank>3,use Dcontract to reshape
	end interface
	interface operator(.pf.)
		module procedure Dpermutefo!Dpermute the inde index to the first
										  !T_{1,2,3,..,i,..,n},Dpermutefo(T,i)=_{i,1,2,3,..,i-1,i+1,..,n}
										  ! note the output will Ddecompose all the dimension
	end interface
	interface operator(.pi.)!T_{1,2,3,..,i,..,n},DpermuteInde(T,i)=_{2,3,..,i,1,i+1,..,n}
									! note the output will Ddecompose all the dimension
		module procedure DpermuteInde
	end interface
		
	interface operator(+)
		module procedure Dadd!the DElement of one DTensor push the other one
		module procedure Dadd_num!one DTensor push an identity DTensor 
		module procedure Dadd_num_!num+ DTensor
		module procedure dconnectlink!combine two DTensor link,the output will be a new link
		module procedure Dconnectlist!combine two list,the output will be a new list
	end interface
	
!**************** Dcontract   ****************
!		combine two index of the DTensor,which is index1 index1+1,..,index2-1,index2
!		if index2 larger than rank,the function will Dcontract the index of index1 -> rank
	interface operator(.c.)
		module procedure Dcontract!combine two index of the DTensor,which is con_index and con_index+1
		module procedure Dcontractsv!combine two index of the DTensor,input vector specify the to be combined index
	end interface
!************************************************************	
!	 Dcompose_decompse do the Dcontract or Ddecompose on the DTensor 
!	(T.cd.Operation)Operation is a matrix,the first DElement of every row specify
!		1:Dcontracts
!		2:Ddecompose
!		3:DdecomposeAll
!	other DElement of every row is input parameter of the function
!		this function run for big T,do as less as possible operation on the
!	 DTensor_Data of T.
	interface operator(.cd.)
		module procedure Dcompose_decompse
	end interface
!************************************************************		
! Ddecompose the de_index index of the DTensor into n(1),n(2)
!		for example the de_index index is (1,2,3,4,..inde,inde+1,...rank)
!		(1,2,3,4,..inde,inde+1,...rank)-->(1,2,3,4,..inde),(inde+1,...rank)		
!		if inde larger than rank ,the function will return no change	
	interface operator(.dc.)
		module procedure Ddecompose1
		module procedure DdecomposeAll
		module procedure Ddecomposev!Ddecompose index of the DTensor,input vector specify the de_index and index
	end interface
	
	interface operator(-)
		module procedure Dminus!the DElement of one DTensor minue the other one
		module procedure Dminus_num!one DTensor minue an identity DTensor 
		module procedure Dminus_num_!num - DTensor
	end interface
!*******************************************************************************
!		productTen is the old version
!		the function productDTensor add at 2013.10.28
	interface operator(*)
		!module procedure productTen!Matrixe product,two DTensor intput should be of 1 or 2 dimension
		module procedure productDTensor!A*B,product the last index of A and first index of B
		module procedure Dmultiply_number!A DTensor times a number,T*num
		module procedure Dmultiply_number_!num*T
		module procedure Dmultiply_real4!A DTensor times a real*4 number
	end interface
	
	interface operator(/)
		module procedure Ddivide!the DElement of one DTensor Ddivide the DElement of other one
		module procedure Ddivide_number!A DTensor Ddivide a number of type real*8
		module procedure Ddivide_numberReal4!A DTensor Ddivide a number of type real*4
	end interface
	
	
	interface operator(.x.)
		module procedure Doudot! dot product ,The Tensor will be regard as a vector
	end interface	
	
	interface operator(.xx.)
		module procedure DdirectProduct! direct Product [m1,n1] * [m2,n2] => [m1,m2,n1,n2]
												!  or direct Product return a matrix [n1] * [n2] => [ n1,n2 ]
	end interface
		
	interface operator(.mxx.)
		module procedure DdirectProductM! direct Product return a matrix [m1,n1] * [m2,n2] => [(m1*m2),(n1*n2)]
												!  or direct Product return a matrix [n1] * [n2] => [ (n1*n2) ]
	end interface
		
	interface operator(.su.)
		module procedure Ddirect_sum! direct sum
	end interface
		
	interface operator(.msu.)
		module procedure Ddirect_sumM! direct sum return a matrix
	end interface
		
	interface operator(.elex.)
		module procedure Dmultiply_DTensor!the DElement of one DTensor times the DElement of other one
	end interface
	
	interface operator(.numx.)
		module procedure DTensorProduct1Dim!two one-dimension DTensors product,return a number of type real*8
	end interface
	
	interface operator(.equ.)
		module procedure equal_of_DTensor!If two DTensor are equal
	end interface
	
	interface operator(.subDim.)
		module procedure DgetTenSubDim!output the ith dimension of the DTensor,return type(Dimension)
		module procedure DgetTenDim!output the Dimension of DTensor,return type(Dimension)
											!It can use to commit the assignment to a vector
											!vector=type(Dimension),if the vector is allocatable
	end interface

	interface operator(.i.)
		module procedure DTi!output the i DTensor in the DTensorlink
		module procedure DElement!output the data in the DTensor of Tim=(i,j)
		module procedure DElement2!output the data in the DTensor of inde
		module procedure DTiI!output the index DTensor in the DTensorlink,index is a vector
		
		module procedure DLi!It gose no wrong when define  DcopyLinkhead as (=),otherwhile it create useless memory
		module procedure DLiI!It gose no wrong when define  DcopyLinkhead as (=),otherwhile it create useless memory
	end interface
	interface operator(.ii.)!the inde(2)th DTensor in inde(1) link
		module procedure DLiv
		module procedure DLiiv
	end interface
	
	interface operator(.dim.)
		module procedure DgetTenDim_i!output the i dimension of the DTensor,output an integer
	end interface
		
	interface assignment(=)
		module procedure DTencopy!T1=T2 ,both T1 and T2 are DTensor
		module procedure Dcopy_dim1!Vec=T ,Vec is a vector and T are DTensor,vec=T%DTensor_Data,used for any case
		module procedure Dcopy_dim2!Mat=T ,Mat is a matrix and T are DTensor,Mat=T%DTensor_Data,used only rank=2
		module procedure Dcopy_dim3!Mat=T ,Mat is a rank=3 and T are DTensor,Mat=T%DTensor_Data,used only rank=3
		module procedure DassignmentTen3!T=vec or matrix
		module procedure DassignmentTen2
		module procedure DassignmentTen1
		module procedure Dcom_Ten!val=T,val is real*8 ,and T is a DTensor,used only there is one DElement in T
		module procedure  DTenLink!DTensor=DTensorlink,if every DElement(DTensor) in the link is a one-DElement DTensor
		
		!module procedure DcopyLink!link1=link2,both link1 and link2 are DTensorLink,and the result ,link1 and link2 will not be the same link
		!module procedure DcopyList!list1=list2,both list1 and list2 are DTensorList
		module procedure DcopyListhead
		module procedure DcopyLinkhead!link1=link2,both link1 and link2 will be the same link
	!if one want to write a function return a link,use DcopyLinkhead as (=)
!     when using a subroutine that return a link,do not use result=link to output,because the data in link will not
!	empty after call the function.But one could use 
!		1. call DcopyLinkhead(result,link) 
!		2. result=link
!			call Dcleanlink(link)
!	the first chose is better,it is faster
	end interface
!	other way to build a DTensor:T=DbuildTen(dimen,data),dimen is a vector or type(Dimension),data is any type
!**********************************************************
!	Other Function or Subroutine:
!
!		generate:generate a DTensor with random number
!
!		cleanDTensor: clean all the data in type(DTensor)
!
!		Print DTensor: 
!			Tprint(DTensor):print all the information of the DTensor
!			TDprint(DTensor):print Dimension
!			TMprint(Tesnor):Print as matrix 
!			DLprint(DTensorlink):print every indice of the DTensor in the link
!			DLDprint(DTensorlink):print every dimension and indice of the DTensor in the link
!
!		get rank or totalData or flag
!			DgetRank
!			DgetTotalData
!			Dgetflag
!
!		equal_of_DTensor:to judge if two DTensor are equak
!
!		Dcontracts((T1,index1,index2)):combine two index of the DTensor,
!			which is index1 index1+1,..,index2-1,index2.if index2 larger 
!			than rank,the function will Dcontract the index of index1 -> rank
!
!		Ddecompose(T1,de_index,inde):Ddecompose the de_index index of the DTensor 
!			into n(1),n(2) for example the de_index index is (1,2,3,4,..inde,inde+1,...rank).
!			(1,2,3,4,..inde,inde+1,...rank)-->(1,2,3,4,..inde),(inde+1,...rank)		
!			if inde larger than rank ,the function will return no change	
!
!		max or min DElement:
!			maxDElement(T):return the max abs Dvalue DElement of the DTensor
!			maxRealDElement(T):return the max real part DElement of the DTensor
!			minDElement(T):return the min abs Dvalue DElement of the DTensor
!			minRealDElement(T):return the min real part DElement of the DTensor
!
!		realTen(T):return the real part of the DTensor
!
!		imagTen(T):return the imag part of the DTensor
!
!		DTensorProduct1Dim(T1,T2): Ddot product,return a comple*16
!
!		DsubTen(T,inde):inde=[-1,inde_min,inde_max] output data(:,inde_min:inde_max)
!			or[-2,inde_min,inde_max],data(inde_min:inde_max,:)
!			or [-1,inde_row] [-2,inde_col],output row or col
!
!		SVD :
!			Dsvdright(T,T1,T2,check_discard): T=U*s*V^T,return the (U*s) 
!				in T1 and V^T in T2 
!			Dsvdleft(T,T1,T2,check_discard):	T=U*s*V^T ,return the U 
!				 in T1,and (s*V^T) in T2
!			Both in Dsvdright and Dsvdleft the data s(i) in the matrix s will 
!				be discard if s(i)<s(1)*check_discard_svd.if check_discard_svd<0
!				,all the data in s will be keep
!
!		linear equations
!			Dlinequ:X=Dlinequ(A,B),solve A*X=B,X could be a matrix
!			Dlinequ_routine(A,B):on output,A will change and X is in B	
!
!		diagonal matrix:
!			Deye_com(s,m,n):return a m*n diagonal matrix will diag(s(1),s(2),...s(slen),0,0....),
!				s(:) are real*8
!			Deye_real::return a m*n diagonal matrix will diag(s(1),s(2),...s(slen),0,0....),
!				s(:) are real*8
!
!		pauli_matrix(Tx,Ty,Tz,num)
!
!		expm(H):
!			return a DTensor e^H,
!
!		Deng(H,vec,val):
!			return a eigenDvalues matrix of H,H=vec*val*vec^H
!
!		sqrt:
!			return sqrt(A),A is a matrix
!
!		Dresetdim:
!			reDset the dimension of the DTensor
!
!		Dconbine
!			T1 :a [l,m,n,...] matrix
!			T2 :a [l,m,n,...] matrix
!			Dconbine(T1,T2):a [2,l,m,n,...] matrix
!			or 
!			T1 :a [l,m,n,...] matrix
!			T2 :a [m,n,...] matrix
!			Dconbine(T1,T2):a [l+1,m,n,...] matrix
!		Dconbinerow(T1,T2)
!			T1 :a [...,l,m,n] matrix
!			T2 :a [...,l,m,n] matrix
!			Dconbine(T1,T2):a [...,l,m,n,2] matrix
!			or 
!			T1 :a [...,m,n,l] matrix
!			T2 :a [...,m,n] matrix
!			Dconbine(T1,T2):a [...,m,n,l+1] matrix
!			
!		operator Dvalue:
!			Dvalue(F,Tr):return the Dvalue of <phi|F|phi>,where F is an operator,
!				F=<I',s1',s2',..|F|I,s1,s2,..>
!			Dnorm2(Tr_):return  <phi|phi>
!			Dnorm(Tr):return  sqrt(<phi|phi>)
!
!		DTensor1(dimen,num):return a DTensor with all the elememt num,dimen is dimension of
!			the DTensor
!		
!		DTen_Ten(T1_,T2_,i1,i2,ifDdecompose):Dcontract The i1 index of T1 and the i2 index of t2
!			T1: (D1,D2,..D_i1,...),T2 :(E1,E2,...,E_i2,....),then the result will be 
!			T:(D1,D2,..D_{i1-1},D_{i1+1},...,D_{rank1},E1,E2,...,E_{i2-1},E_{i2+1},...,E_{rank2})	
!			if if Ddecompose=1, T will be Ddecompose
!
!		DTensor index
!			DTenIndice(h,inde,indice):return the index of the inde DTensor.T_{1,1}->T_{1,2}->T_{1,3}
!				->T_{2,1}->T_{2,2},if input inde=2 => indice=[1,2].on entry the size of indice should
!				 be equal to the one in h
!			DTenIndiceLog(h,inde,indice):if there no index in h,return .false..other while return inde
!				of the DTensor
!			DlenOfIndicelist(h,inde,indice):return the index of the inde DTensorlink.T_{1,1}->T_{1,2}->T_{1,3}
!				->T_{2,1}->T_{2,2},if input inde=2 => indice=[1,2].on entry the size of indice should
!				 be equal to the one in h
!			DlinkIndiceLog(h,inde,indice):if there no index in h,return .false..other while return inde
!				of the DTensorlink
!
!		pointer subroutine
!			DifDnextnode:if the next node exict,p point to the next node,or p remain nochange
!			Dheadnode:make the pointer points to the head of the link
!			Dnextnode:make the pointer points to the next node
!			Dendnode:make the pointer points to the end of the link
!			Daddnode:add a empty node to the link(then use a pointer to point to it)
!
!			Difnextlistnode:if the next node exict,p point to the next node,or p remain nochange(DTensorlist subroutine)
!			headlist:make the pointer points to the head of the list(DTensorlist subroutine)
!			nextDlistnode:make the pointer points to the next node(DTensorlist subroutine)
!			endlist:make the pointer points to the end of the list(DTensorlist subroutine)
!			addlist:add a empty node to the link(then use a pointer to point to it)(DTensorlist subroutine)
!
!		DChecklength(link):check and output length of the link,check if the length of the link is equal 
!				to the Dvalue in head of the link
!		DChecklistlength(list):check and output length of the list,check if the length of the list is equal 
!				to the Dvalue in head of the list
!
!		dlinklength(link):output length of the link
!		listlength(list):output length of the list
!
!		DcopyLinkhead:copy the head of a link
!
!		Dmodify(h,T,inde):Dmodify the inde DElement of the link.The DTensor will be Dmodify
!
!		Dcleanlink(h):clean the link h
!		Dcleanlist(h):clean the list h
!
!		Ddeletelink(Dlink_in,inde1,inde2):
!			link   :	T1->T2->...->T_{inde1-1}->T_inde1->T_{inde1+1}->....->T_{inde2-1}->T_inde2->T_{inde2+1}...->TN 
!		  result:		T1->T2->...->T_{inde1-1}->T_{inde1+1}....->T_{inde2-1}->->T_{inde2+1}...->TN 
!		  delete the inde1 to inde2 DTensor in the link,note that the deleted DTensor are include the DTensors of inde1 and inde2
!		Ddeletelist(list_in,inde1,inde2):(DTensorlist subroutine)
!			list   :	T1->T2->...->T_{inde1-1}->T_inde1->T_{inde1+1}->....->T_{inde2-1}->T_inde2->T_{inde2+1}...->TN 
!		  result:		T1->T2->...->T_{inde1-1}->T_{inde1+1}....->T_{inde2-1}->->T_{inde2+1}...->TN 
!		  delete the inde1 to inde2 DTensorlink in the list,note that the deleted DTensorlink are include the DTensorlinks of inde1 and inde2
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
	subroutine DTencopy(T,T2)
		type(DTensor),intent(out) ::T
		type(DTensor),intent(in) :: T2
		if(T2%flag) then
			T%rank=T2%rank
			T%TenDim=T2%TenDim
			call DstoreTenData(T,T2%DTensor_Data)
			T%totalData=T2%totalData
			T%flag=T2%flag
		else
			T%flag=T2%flag
		endif
		return
	end subroutine
	subroutine Dcom_Ten(val,T)
		real*8,intent(out)::val
		type(DTensor),intent(in)::T
		if(T%totalData.eq.1) then
			val=T.i.1
		else
			write(*,*)"ERROR in assignment for DTensor to complex"
			stop
		end if
		return
	end subroutine
	subroutine DassignmentTen1(T,DTensor_data)
		real*8,intent(in)::DTensor_data(:)
		type(DTensor),intent(inout)::T
		call cleanDTensor(T)
		T=DbuildTen(DTensor_data)
		return
	end subroutine
	subroutine DassignmentTen2(T,DTensor_data)
		real*8,intent(in)::DTensor_data(:,:)
		type(DTensor),intent(inout)::T
		call cleanDTensor(T)
		T=DbuildTen(DTensor_data)
		return
	end subroutine
	subroutine DassignmentTen3(T,DTensor_data)
		real*8,intent(in)::DTensor_data(:,:,:)
		type(DTensor),intent(inout)::T
		call cleanDTensor(T)
		T=DbuildTen(DTensor_data)
		return
	end subroutine
		
!*************** cleanDTensor  *****************
	subroutine cleanDTensor(T)
		Type(DTensor),intent(inout)::T
		T%rank=0
		call cleanDimension(T%TenDim)
		if(allocated(T%DTensor_data)) then
			deallocate(T%DTensor_data)
		end if
		T%totalData=0
		T%flag=.false.
		return
	end subroutine
!copy the DTensor_Data to a vector	
	subroutine Dcopy_dim1(Vec,T)
		real*8,allocatable,intent(out) ::Vec(:)
		type(DTensor),intent(in) :: T
		integer::length
		length=T%TotalData
		if(allocated(Vec)) then
			deallocate(Vec)
		end if
		allocate(Vec(length))
		call dcopy(length,T%DTensor_Data,1,Vec,1)
		return
	end subroutine
!if rank=2,then copy the DTensor_Data to a mat
! mat is a matrix	
	subroutine Dcopy_dim2(Mat,T)
		real*8,allocatable,intent(out) ::Mat(:,:)
		type(DTensor),intent(in) :: T
		integer::m,n
		if(T%rank.eq.2) then
			m=T.dim.1
			n=T.dim.2
			if(allocated(Mat)) then
				deallocate(Mat)
			end if
			allocate(Mat(m,n))
			call dcopy(m*n,T%DTensor_Data,1,Mat,1)
		else
			write(*,*)"The DTensor is not a matrix"
			write(*,*)"Assignment ERROR,stop"
			stop
		end if
		return
	end subroutine
!if rank=3,then copy the DTensor_Data to a mat
!mat is a 3 dimension array	
	subroutine Dcopy_dim3(Mat,T)
		real*8,allocatable,intent(out) ::Mat(:,:,:)
		type(DTensor),intent(in) :: T
		integer::m,n,l
		if(T%rank.eq.3) then
			m=T.dim.1
			n=T.dim.2
			l=T.dim.3
			if(allocated(Mat)) then
				deallocate(Mat)
			end if
			allocate(Mat(m,n,l))
			call dcopy(m*n*l,T%DTensor_Data,1,Mat,1)
		else
			write(*,*)"Assignment ERROR,stop"
			stop
		end if
		return
	end subroutine	
!***********	DstoreTenData  **************
	subroutine DstoreTenData_dim1(T,inData)
		real*8,intent(in)::inData(:)
		type(DTensor),intent(inout)::T
		integer::length
		length=size(inData)
		if(allocated(T%DTensor_Data)) then
			deallocate(T%DTensor_Data)
		end if
		allocate(T%DTensor_Data(length))
		call dcopy(length,inData,1,T%DTensor_Data,1)
		return
	end subroutine
	subroutine DstoreTenData_dim2(T,inData)
		real*8,intent(in)::inData(:,:)
		type(DTensor),intent(inout)::T
		integer::length
		length=size(inData)
		if(allocated(T%DTensor_Data)) then
			deallocate(T%DTensor_Data)
		end if
		allocate(T%DTensor_Data(length))
		call dcopy(length,inData,1,T%DTensor_Data,1)
		return
	end subroutine
	subroutine DstoreTenData_dim3(T,inData)
		real*8,intent(in)::inData(:,:,:)
		type(DTensor),intent(inout)::T
		integer::length
		length=size(inData)
		if(allocated(T%DTensor_Data)) then
			deallocate(T%DTensor_Data)
		end if
		allocate(T%DTensor_Data(length))
		call dcopy(length,inData,1,T%DTensor_Data,1)
		return
	end subroutine
	
!***************** DbuildTen *********************
!	build a DTensor
	type(DTensor) function DbuildTen1(rank,TenDim,DTensor_data)
		integer,intent(in)::rank
		type(Dimension),intent(in)::TenDim
		real*8,intent(in)::DTensor_data(:)
		DbuildTen1%rank=rank
		DbuildTen1%TenDim=TenDim
		call DstoreTenData(DbuildTen1,DTensor_Data)
		DbuildTen1%totalData=size(DTensor_Data)
		DbuildTen1%flag=.true.
		return
	end function	
	type(DTensor) function DbuildTen2(rank,TenDim,DTensor_data)
		integer,intent(in)::rank
		integer,intent(in)::TenDim(:)
		real*8,intent(in)::DTensor_data(:)
		DbuildTen2%rank=rank
		DbuildTen2%TenDim=TenDim
		call DstoreTenData(DbuildTen2,DTensor_Data)
		DbuildTen2%totalData=size(DTensor_Data)
		DbuildTen2%flag=.true.
		return
	end function		
	type(DTensor) function DbuildTen3(TenDim,DTensor_data)
		integer,intent(in)::TenDim(:)
		real*8,intent(in)::DTensor_data(:)
		DbuildTen3%rank=size(TenDim)
		DbuildTen3%TenDim=TenDim
		call DstoreTenData(DbuildTen3,DTensor_Data)
		DbuildTen3%totalData=size(DTensor_Data)
		DbuildTen3%flag=.true.
		if(product(TenDim).ne.size(DTensor_Data)) then
			write(*,*)"ERROR in dimension and input data in function DbuildTen"
			write(*,*)"stop"
		end if
		return
	end function	
	type(DTensor) function DbuildTen4(DTensor_data)
		real*8,intent(in)::DTensor_data(:)
		DbuildTen4%rank=1
		DbuildTen4%TenDim=(/size(DTensor_Data)/)
		call DstoreTenData(DbuildTen4,DTensor_Data)
		DbuildTen4%totalData=size(DTensor_Data)
		DbuildTen4%flag=.true.
		return
	end function
	type(DTensor) function DbuildTen5(DTensor_data)
		real*8,intent(in)::DTensor_data(:,:)
		DbuildTen5%rank=2
		DbuildTen5%TenDim=(/size(DTensor_Data,1),size(DTensor_Data,2)/)
		call DstoreTenData(DbuildTen5,DTensor_Data)
		DbuildTen5%totalData=size(DTensor_Data)
		DbuildTen5%flag=.true.
		return
	end function		
	type(DTensor) function DbuildTen6(DTensor_data)
		real*8,intent(in)::DTensor_data(:,:,:)
		DbuildTen6%rank=3
		DbuildTen6%TenDim=(/size(DTensor_Data,1),size(DTensor_Data,2),size(DTensor_Data,3)/)
		call DstoreTenData(DbuildTen6,DTensor_Data)
		DbuildTen6%totalData=size(DTensor_Data)
		DbuildTen6%flag=.true.
		return
	end function		
	type(DTensor) function DbuildTen8(TenDim,DTensor_data)
		type(Dimension),intent(in)::TenDim
		integer,allocatable::Dimen(:)
		real*8,intent(in)::DTensor_data(:)
		Dimen=TenDim
		DbuildTen8%rank=size(Dimen)
		DbuildTen8%TenDim=TenDim
		call DstoreTenData(DbuildTen8,DTensor_Data)
		DbuildTen8%totalData=size(DTensor_Data)
		DbuildTen8%flag=.true.
		
		if(product(Dimen).ne.size(DTensor_Data)) then
			write(*,*)"ERROR in dimension and input data in function DbuildTen"
			write(*,*)"stop"
		end if
		return
	end function	

!********************* Dset *********************
!   	Dset the data in a DTensor 	 	c
	function Dset1(dim_,Ten_data) result (T)
		type(DTensor) :: T
		real*8,intent(in) :: Ten_data(:)
		integer,intent(in)	:: dim_(:)
		integer :: total,i,rank,test_total
		type(Dimension) :: TenDim
		rank=size(dim_)
		total=size(Ten_data)
		test_total=product(dim_)
		if(test_total.ne.total) then
			write(*,*)"ERROR in Dset1,stop"
			stop
		end if
		TenDim=dim_
		call DstoreTenData(T,Ten_data)
		T%rank=rank
		T%totalData=total
		T%TenDim=TenDim
		T%flag=.true.
		return
	end function
		 
	function Dset2(TenDim,Ten_data) result (T)
		type(DTensor) :: T
		real*8,intent(in) :: Ten_data(:)
		type(Dimension),intent(in)	:: TenDim
		integer :: total,i,rank
		rank=DimSize(TenDim)
		total=size(Ten_data)
!		T%DTensor_data=Ten_data
!		call assignment_com_dim1(T%DTensor_data,Ten_data)
		call DstoreTenData(T,Ten_data)
		T%rank=rank
		T%totalData=total
		T%TenDim=TenDim
		T%flag=.true.
		return
	end function
!*********************  DgetRank	 **********************
	integer function DgetRank(T)
		type(DTensor),intent(in) :: T
		DgetRank=T%rank
	end function
!*********************  DgetTotalData	 **********************
	integer function DgetTotalData(T)
		type(DTensor),intent(in) :: T
		DgetTotalData=T%totalData
	end function	
!*********************  Dgetflag	 **********************
	logical function Dgetflag(T)
		type(DTensor),intent(in) :: T
		Dgetflag=T%flag
	end function		
!*********************  generate *********************
!		generate a DTensor with random number
	type(DTensor) function Dgenerate(Tdim) result (T)
		integer,intent(in) :: Tdim(:)
		integer :: rank,totalData,i
		real*8,allocatable :: DTensor_data(:)
		real*8 ::temp_real
		type(Dimension):: TenDim
		call init_random_seed()
		rank=size(Tdim)
		totalData=product(Tdim)
		allocate(DTensor_data(totalData))
		do i=1,totalData
			call random_number(temp_real)
			DTensor_data(i)=temp_real
		end do
		TenDim=Tdim
		call DstoreTenData(T,DTensor_data)
		T%rank=rank
		T%totalData=totalData
		T%TenDim=TenDim
		T%flag=.true.
		return
	end function
!*********************  zeroTen *********************
!		generate a Tensor with all element is 0d0
	type(DTensor) function DzeroTen(Tdim) result (T)
		integer,intent(in) :: Tdim(:)
		integer :: rank,totalData,i
		rank=size(Tdim)
		totalData=product(Tdim)
		allocate(T%DTensor_data(totalData))
		T%DTensor_data=0D0
		T%rank=rank
		T%totalData=totalData
		T%TenDim=Tdim
		T%flag=.true.
		return
	end function
!*********************  get DTensor dimenison *********************
	integer function DgetTenDim_i(T,inde)
		type(DTensor),intent(in) :: T
		integer,intent(in) :: inde
		DgetTenDim_i=Dim_i(T%TenDim,inde)
		return
	end function
	type(Dimension) function DgetTenDim(T)
		type(DTensor),intent(in) :: T
		DgetTenDim=T%TenDim
		return
	end function
	type(Dimension) function DgetTenSubDim(T,inde)
		type(DTensor),intent(in)::T
		integer,intent(in)::inde
		DgetTenSubDim=T%TenDim.sub.inde
		return
	end function
!********************* print_DTensor *********************       
	subroutine DTprint1(T)
		type(DTensor),intent(in) :: T
		if(T%flag) then
			write(*,*) "The rank of the DTensor is"
			write(*,*) T%rank
			write(*,*) "The number of  data of the DTensor is"
			write(*,*) T%totalData
			write(*,*) "The data of the DTensor is"
			write(*,*) T%DTensor_data
			call Dprint0(T%TenDim)
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
!********************* print_dim *********************            
	subroutine DTDprint(T)
		type(DTensor) :: T
		if(T%flag) then
			write(*,*) "*** START ***"
			write(*,*) "The rank of the DTensor is"
			write(*,*) T%rank
			write(*,*) "The number of  data of the DTensor is"
			write(*,*) T%totalData
			call Dprint0(T%TenDim)
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
	subroutine DTprint2(T,printtype)
		type(DTensor),intent(in) :: T
		integer,intent(in)::printtype
		if(T%flag) then
			write(*,*) "The rank of the DTensor is"
			write(*,*) T%rank
			write(*,*) "The number of  data of the DTensor is"
			write(*,*) T%totalData
			write(*,*) "The data of the DTensor is"
			select case (printtype)
			case (1)
				write(*,*) real(T%DTensor_data)
			case (2)
				write(*,'(99999F12.8)') real(T%DTensor_data)
			case (3)
				write(*,'(99999F10.2)') real(T%DTensor_data)
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
	subroutine DTMprint1(T)
		type(DTensor) T
		real*8,allocatable :: Tdata(:,:)
		real*8,allocatable :: Tdata3(:,:,:)
		real*8,allocatable :: Tdata4(:,:,:,:)
		integer,allocatable ::Tdim(:)
		integer :: i,j,k
		if(T%flag) then!if 1
			write(*,*) "*** START ***"
			select case(T%rank)
				case(1)
					write(*,*) T%DTensor_Data
					write(*,*) "*** END ***"
				case(2)
					allocate(Tdim(2))
					Tdim=T%TenDim
					allocate(Tdata(Tdim(1),Tdim(2)))
					Tdata=reshape(T%DTensor_Data,shape(Tdata))
						do i=1,Tdim(1)
							write (*,*) Tdata(i,:)
							write(*,*) ""
						end do 
						write(*,*) "*** END ***"
				case(3)
					allocate(Tdim(3))
					Tdim=T%TenDim
					allocate(Tdata3(Tdim(1),Tdim(2),Tdim(3)))
					Tdata3=reshape(T%DTensor_Data,shape(Tdata3))
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
					Tdata4=reshape(T%DTensor_Data,shape(Tdata4))
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
					write(*,*) "rank of the DTensor is large than 4"
					write(*,*) "The data of the DTensor is"
					write(*,*) T%DTensor_data
					write(*,*) "***end***"
					write(*,*) "The dimension of the DTensor is"
					call Dprint0(T%TenDim)
					write(*,*) "The rank,total data are"
					write(*,*) T%rank,T%totalData
			end select
		else!if 1
			write(*,*) "There is no data"
		end if!if 1
		return
	end subroutine
	subroutine DTMprint2(T,printtype)
		type(DTensor),intent(in) :: T
		integer,intent(in)::printtype
		real*8,allocatable :: Tdata(:,:)
		real*8,allocatable :: Tdata3(:,:,:)
		real*8,allocatable :: Tdata4(:,:,:,:)
		integer,allocatable ::Tdim(:)
		integer :: i,j,k
		if(T%flag) then!if 1
			write(*,*) "*** START ***"
			select case(T%rank)
				case(1)
					select case (printtype)
					case (1)
						write(*,*) real(T%DTensor_data)
					case (2)
						write(*,'(99999F12.8)') real(T%DTensor_data)
					case (3)
						write(*,'(99999F10.2)') real(T%DTensor_data)
					end 	select
					write(*,*) "*** END ***"
				case(2)
					allocate(Tdim(2))
					Tdim=T%TenDim
					allocate(Tdata(Tdim(1),Tdim(2)))
					Tdata=reshape(T%DTensor_Data,shape(Tdata))
					select case (printtype)
					case (1)
						do i=1,Tdim(1)
							write (*,*) real(Tdata(i,:))
							write(*,*) ""
						end do 
					case (2)
						do i=1,Tdim(1)
							write (*,'(99999F12.8)') real(Tdata(i,:))
							write(*,*) ""
						end do 
					case (3)
						do i=1,Tdim(1)
							write (*,'(99999F10.2)') real(Tdata(i,:))
							write(*,*) ""
						end do 
					end 	select
						write(*,*) "*** END ***"
				case(3)
					allocate(Tdim(3))
					Tdim=T%TenDim
					allocate(Tdata3(Tdim(1),Tdim(2),Tdim(3)))
					Tdata3=reshape(T%DTensor_Data,shape(Tdata3))
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
								write (*,'(99999F12.8)') real(Tdata3(i,:,j))
								write(*,*) ""
							end do
						end do 
					case (3)
						do j=1,Tdim(3)
							write(*,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write (*,'(99999F10.2)') real(Tdata3(i,:,j))
								write(*,*) ""
							end do
						end do 
					end 	select
					write(*,*) "*** END ***"
				case(4)
					allocate(Tdim(4))
					Tdim=T%TenDim
					allocate(Tdata4(Tdim(1),Tdim(2),Tdim(3),Tdim(4)))
					Tdata4=reshape(T%DTensor_Data,shape(Tdata4))
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
								write (*,'(99999F12.8)') real(Tdata4(i,:,k,j))
								write(*,*) ""
								end do
							end do 
						end do
					case (3)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(*,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write (*,'(99999F10.2)') real(Tdata4(i,:,k,j))
								write(*,*) ""
								end do
							end do 
						end do
					end 	select
					write(*,*) "*** END ***"
				case default
					call DTprint2(T,printtype)
			end select
		else!if 1
			write(*,*) "There is no data"
		end if!if 1
		return
	end subroutine


!cccccccccccccccc add          cccccccccccccccccc
	type(DTensor) function Dadd(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		type(Dimension)::dim1,dim2
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*)"The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		Dadd=Dset(dim1,T1%DTensor_data+T2%DTensor_data)
		return
	end function
	type(DTensor) function Dadd_num(T1,num)
		type(DTensor),intent(in) :: T1
		real*8,intent(in)::num
		type(DTensor)::T2
		integer,allocatable::dimT(:)
		dimT=T1%TenDim
		T2=Deye(dimT)
		Dadd_num=Dset(T1%TenDim,T1%DTensor_data+T2%DTensor_data)
		return
	end function
	type(DTensor) function Dadd_num_(num,T1)
		type(DTensor),intent(in) :: T1
		real*8,intent(in)::num
		type(DTensor)::T2
		integer,allocatable::dimT(:)
		dimT=T1%TenDim
		T2=Deye(dimT)
		Dadd_num_=Dset(T1%TenDim,T1%DTensor_data+T2%DTensor_data)
		return
	end function
!cccccccccccccccc    Dminus       cccccccccccccccccc		
	type(DTensor) function Dminus(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		type(Dimension)::dim1,dim2
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*)"The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		Dminus=Dset(dim1,T1%DTensor_data-T2%DTensor_data)
		return
	end function
	type(DTensor) function Dminus_num(T1,num)
		type(DTensor),intent(in) :: T1
		real*8,intent(in)::num
		type(DTensor)::T2
		integer,allocatable::dimT(:)
		dimT=T1%TenDim
		T2=Deye(dimT)
		Dminus_num=Dset(T1%TenDim,T1%DTensor_data-T2%DTensor_data)
		return
	end function
	type(DTensor) function Dminus_num_(num,T1)
		type(DTensor),intent(in) :: T1
		real*8,intent(in)::num
		type(DTensor)::T2
		integer,allocatable::dimT(:)
		dimT=T1%TenDim
		T2=Deye(dimT)
		Dminus_num_=Dset(T1%TenDim,T1%DTensor_data-T2%DTensor_data)
		return
	end function
!cccccccccccccccc Ddivide          cccccccccccccccccc	
	type(DTensor) function Ddivide(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		type(Dimension)::dim1,dim2
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*) "The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		Ddivide=Dset(dim1,T1%DTensor_data/T2%DTensor_data)
		return
	end function
!cccccccccccccccc Ddivide_number          cccccccccccccccccc	
	type(DTensor) function Ddivide_number(T1,num)
		type(DTensor),intent(in) :: T1
		real*8,intent(in) :: num
		Ddivide_number%rank=T1%rank
		Ddivide_number%TenDim=T1%TenDim
		call DstoreTenData(Ddivide_number,T1%DTensor_data/num)
		Ddivide_number%totalData=T1%totalData
		Ddivide_number%flag=.true.
		return
	end function
!ccccccccccccccccc Ddivide_number          cccccccccccccccccc	
	type(DTensor) function Ddivide_numberReal4(T1,num)
		type(DTensor),intent(in) :: T1
		real*4,intent(in) :: num
		Ddivide_numberReal4%rank=T1%rank
		Ddivide_numberReal4%TenDim=T1%TenDim
		call DstoreTenData(Ddivide_numberReal4,T1%DTensor_data/num)
		Ddivide_numberReal4%totalData=T1%totalData
		Ddivide_numberReal4%flag=.true.
		return
	end function			
!cccccccccccccccc Dmultiply_DTensor          cccccccccccccccccc	
	type(DTensor) function Dmultiply_DTensor(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		type(Dimension)::dim1,dim2
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*) "The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		Dmultiply_DTensor=Dset(dim1,T1%DTensor_data*T2%DTensor_data)
		return
	end function
!cccccccccccccccc Dmultiply_number          cccccccccccccccccc	
	type(DTensor) function Dmultiply_number(T1,num)
		type(DTensor),intent(in) :: T1
		real*8,intent(in) ::   num
		Dmultiply_number%rank=T1%rank
		Dmultiply_number%TenDim=T1%TenDim
		call DstoreTenData(Dmultiply_number,T1%DTensor_data*num)
		Dmultiply_number%totalData=T1%totalData
		Dmultiply_number%flag=.true.
		return
	end function
		type(DTensor) function Dmultiply_number_(num,T1)
		type(DTensor),intent(in) :: T1
		real*8,intent(in) ::   num
		Dmultiply_number_%rank=T1%rank
		Dmultiply_number_%TenDim=T1%TenDim
		call DstoreTenData(Dmultiply_number_,T1%DTensor_data*num)
		Dmultiply_number_%totalData=T1%totalData
		Dmultiply_number_%flag=.true.
		return
	end function
!ccccccccccccccccc Dmultiply_real4          cccccccccccccccccc	
	type(DTensor) function Dmultiply_real4(T1,num)
		type(DTensor),intent(in) :: T1
		real*4,intent(in) ::   num
		Dmultiply_real4%rank=T1%rank
		Dmultiply_real4%TenDim=T1%TenDim
		call DstoreTenData(Dmultiply_real4,T1%DTensor_data*num)
		Dmultiply_real4%totalData=T1%totalData
		Dmultiply_real4%flag=.true.
		return
	end function			
!ccccccccccccccccc Dpermute_rank3   cccccccccccccccccc		
	type(DTensor) function Dpermute_rank3(T1,index_not_Dpermute)
		type(DTensor),intent(in) :: T1
		integer,intent(in) ::   index_not_Dpermute
		integer ::i,newdimen(3),totaldata
		integer,allocatable :: dimen(:)
		type(Dimension) ::newTDim
		real*8,allocatable :: Tdata(:,:,:),newdata(:,:,:)
		dimen=T1%TenDim
		allocate(Tdata(dimen(1),dimen(2),dimen(3)))
		Tdata=T1
		if(index_not_Dpermute.eq.1) then
			newdimen(1)=dimen(1)
			newdimen(2)=dimen(3)
			newdimen(3)=dimen(2)
			newTDim=Dimpermute(T1%TenDim,(/1,3,2/))
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(1)
				newdata(i,:,:)=transpose(Tdata(i,:,:))
			end do
			call DstoreTenData(Dpermute_rank3,newdata)
			Dpermute_rank3%rank=3
			Dpermute_rank3%totalData=T1%totalData
			Dpermute_rank3%TenDim=newTDim
			Dpermute_rank3%flag=.true.
		end if
		if(index_not_Dpermute.eq.2) then
			newdimen(1)=dimen(3)
			newdimen(2)=dimen(2)
			newdimen(3)=dimen(1)
			newTDim=Dimpermute(T1%TenDim,(/3,2,1/))
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(2)
				newdata(:,i,:)=transpose(Tdata(:,i,:))
			end do
			call DstoreTenData(Dpermute_rank3,newdata)
			Dpermute_rank3%rank=3
			Dpermute_rank3%totalData=T1%totalData
			Dpermute_rank3%TenDim=newTDim
			Dpermute_rank3%flag=.true.
		end if
		if(index_not_Dpermute.eq.3) then
			newdimen(1)=dimen(2)
			newdimen(2)=dimen(1)
			newdimen(3)=dimen(3)
			newTDim=Dimpermute(T1%TenDim,(/2,1,3/))
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(3)
				newdata(:,:,i)=transpose(Tdata(:,:,i))
			end do
			call DstoreTenData(Dpermute_rank3,newdata)
			Dpermute_rank3%rank=3
			Dpermute_rank3%totalData=T1%totalData
			Dpermute_rank3%TenDim=newTDim
			Dpermute_rank3%flag=.true.
		end if
	 end function

!cccccccccccccccc Dpermute_rank2   cccccccccccccccccc			 
	type(DTensor) function Dpermute_rank2 (T)
		type(DTensor),intent(in) :: T
		real*8,allocatable :: Tdata(:,:),newdata(:,:)
		integer,allocatable ::dimen(:)
		type(Dimension) ::newTDim
		dimen=T%TenDim
		newTDim=Dimpermute(T%TenDim,(/2,1/))
		allocate(Tdata(dimen(1),dimen(2)))
		Tdata=T
		allocate(newdata(dimen(2),dimen(1)))
		newdata=transpose(Tdata)
		call DstoreTenData(Dpermute_rank2,newdata)
		Dpermute_rank2%rank=2
		Dpermute_rank2%totalData=T%totalData
		Dpermute_rank2%TenDim=newTDim
		Dpermute_rank2%flag=.true.
	end function
!cccccccccccccccc Dpermutation   cccccccccccccccccc
	type(DTensor) function Dpermutation(T,newOrder)
		type(DTensor),intent(in) :: T
		integer,intent(in)::newOrder(:)
		integer,allocatable ::inde(:),newinde(:),maxdim(:),newmaxdim(:)
		integer::i,total,j,lendim,newaddress
		logical::fastoutput!in no Dpermutation,output
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
			Dpermutation=T
			return
		end if
		total=T%TotalData
		Dpermutation%TenDim=Dimpermute(T%TenDim,newOrder)
		allocate(Dpermutation%DTensor_Data(total))
		maxdim=T%TenDim
		newmaxdim=Dpermutation%TenDim
		allocate(newinde(lendim))
		i=1
		do i=1,total
			call IndesToaddress(maxdim,inde,i)
			do j=1,lendim
				newinde(j)=inde(newOrder(j))
			end do
			newaddress=addressToIndes2(newmaxdim,newinde)
			Dpermutation%DTensor_Data(newaddress)=T%DTensor_Data(i)
		end do
		Dpermutation%rank=T%rank
		Dpermutation%totalData=total
		Dpermutation%flag=T%flag
	end function

!ccccccccccccccccc   productTen   ccccccccccccccccccc
!			the dimension of T1 and T2 should be 2 or 1
!			choose the righ form of
!						DTensorProduct
!					DTensorProduct1Dim
!					DTensorProduct1Dim2
!		this is the old version.ProductDTensor can do the same with
!	  higher perfromence.
	type(DTensor)	function DproductTen(T1,T2)
		type(DTensor),intent(in) :: T1
		type(DTensor),intent(in) :: T2
		integer :: rank1,rank2
		character*50 ::wt
		rank1=T1%rank
		rank2=T2%rank
		if((rank1.le.2).and.(rank2.le.2)) then
			if((rank1.eq.2).and.(rank2.eq.2)) then
				DproductTen=DTensorProduct(T1,T2)
			else if((rank1.eq.1).and.(rank2.eq.1)) then
				DproductTen=DTensorProduct1Dim1(T1,T2)
			else if((rank1.eq.1).or.(rank2.eq.1)) then
				DproductTen=DTensorProduct1Dim2 (T1,T2) 
			else
				write(*,*) "error in the dimension of the DTensor"
			end if
		else
		end if
		return
	end function
				 	
!cccccccccccccccc DTensorProduct   cccccccccccccccccc			 
!		both T1 and T2 should be a 2 \times 2 DTensor
	type(DTensor) function DTensorProduct(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		integer :: totalData
		integer ::T1m,T1n,T2m,T2n
		type(Dimension) :: TenDim,dim1,dim2,newdim
		if((T1%rank.eq.2).and.(T2%rank.eq.2)) then
			T1m=T1.dim.1
			T1n=T1.dim.2
			T2m=T2.dim.1
			T2n=T2.dim.2
			totalData=T1m*T2n
			allocate(DTensorProduct%DTensor_data(totalData))
			call DGEMM('N','N',T1m,T2n,T1n,1D0,T1%DTensor_Data,T1m,&
				T2%DTensor_Data,T2m,0D0,DTensorProduct%DTensor_data,T1m)
			dim1=T1%TenDim.sub.1
			dim2=T2%TenDim.sub.2
			newdim=dim1+dim2
			DTensorProduct%rank=2
			DTensorProduct%totalData=totalData
			DTensorProduct%TenDim=newdim
			DTensorProduct%flag=.true.
		else
			write(*,*)"The two DTensor input for product Should both be matrix"
		end if
	end function
!cccccccccccccccc  DTensorProduct1Dim  cccccccccccccccccc
	real*8 function DTensorProduct1Dim(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: T1data(:,:)
		real*8,allocatable :: T2data(:,:)
		real*8 :: Ndata
		integer ::T1m,T1n,T2m,T2n
		if((T1%rank.eq.1).and.(T2%rank.eq.1)) then
			T1n=T1.dim.1
			T2m=T1.dim.2
			if(T1n.ne.T2m) then
				write(*,*)"ERROR in DTensorProduct1Dim1,stop"
				stop
			end if
			Ndata=Ddot(T1n,T1%DTensor_Data,1,T2%DTensor_Data,1)
			DTensorProduct1Dim=Ndata
		else
			write(*,*)"The two DTensor input for product Should both be a vector"
		end if
	end function
!cccccccccccccccc  DTensorProduct1Dim1  cccccccccccccccccc				 
	type(DTensor) function DTensorProduct1Dim1(T1,T2)  
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: T1data(:,:)
		real*8,allocatable :: T2data(:,:)
		real*8 :: Ndata(1,1)
		integer ::T1m,T1n,T2m,T2n
		if((T1%rank.eq.1).and.(T2%rank.eq.1)) then
			T1n=T1.dim.1
			T2m=T2.dim.1
			if(T1n.ne.T2m) then
				write(*,*)"ERROR in DTensorProduct1Dim1,stop"
				stop
			end if
			Ndata=Ddot(T1n,T1%DTensor_Data,1,T2%DTensor_Data,1)
			DTensorProduct1Dim1=Dset((/1/),(/Ndata/))
		else
			write(*,*)"The two DTensor input for product Should both be a vector"
		end if
		return
	end function
!cccccccccccccccc	 DTensorProduct1Dim2   ccccccccccccccccc
	type(DTensor) function DTensorProduct1Dim2 (T1,T2) result(ResultT)
		integer :: m,n,tempDim(2),RDim(1)
		type(DTensor),intent(in) :: T1,T2
		type(DTensor) :: temp
		type(Dimension) ::diml,dimr,newdim1,newdim2
		real*8,allocatable::reData(:)
		if((T1%rank.eq.1).and.(T2%rank.eq.2)) then
			M=T2.dim.1
			N=T2.dim.2
			if((T1.dim.1) .ne. M) then
				write(*,*)"ERROR in DTensorProduct1Dim2,stop"
				stop
			end if
			allocate(reData(N))
			call DGEMV('T',M,N,1D0,T2%DTensor_Data,M,T1%DTensor_Data,1,0D0,reData,1)
			ResultT=Dset(T2%TenDim.sub.2,reData)
		else if((T1%rank.eq.2).and.(T2%rank.eq.1)) then
			M=T1.dim.1
			N=T1.dim.2
			if((T2.dim.1) .ne. N) then
				write(*,*)"ERROR in DTensorProduct1Dim2,stop"
				stop
			end if
			allocate(reData(M))
			call DGEMV('N',M,N,1D0,T1%DTensor_Data,M,T2%DTensor_Data,1,0D0,reData,1)
			ResultT=Dset(T1%TenDim.sub.1,reData)
		else
			write(*,*)"Error in the dimension of the input DTensors"
		end if
		return
	end function
!**************** ProductDTensor  ***************************
!		ProductDTensor regard the last index of T1 and the first index
!	  of T2 as the dimenion for matrix-product,other index will be see
!	  as another dimenison.T1 and T2 can be any rank,but the last dimenion
!	  of T1 and the first diemsion of T2 should be equal.
	type(DTensor) function ProductDTensor (T1,T2)
		type(DTensor),intent(in) :: T1,T2
		integer::rank1,rank2,flag,T1m,T1n,T2m,T2n,T1l,T2l
		type(Dimension)::D1,D2,newD
		real*8::Ndata
		integer::i,totaldata
		rank1=DgetRank(T1)
		rank2=DgetRank(T2)
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
			write(*,*)"ERROR in ProductDTensor",rank1,rank2
			stop
		end if
		select case (flag)
			case (1)
				T1m=T1.dim.1
				T2n=T2.dim.1
				if(T1m.ne.T2n) then
					write(*,*)"ERROR in ProductDTensor,case 1,stop"
					stop
				end if
				Ndata=Ddot(T1m,T1%DTensor_Data,1,T2%DTensor_Data,1)
				ProductDTensor=Dset((/1/),(/Ndata/))
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
					write(*,*)"ERROR in ProductDTensor,case 2,stop"
					stop
				end if
				allocate(ProductDTensor%DTensor_Data(T2n))
				call DGEMV('T',T2m,T2n,1D0,T2%DTensor_Data,T2m,T1%DTensor_Data,&
							1,0D0,ProductDTensor%DTensor_Data,1)
				ProductDTensor%totalData=T2n
				ProductDTensor%TenDim=newD
				ProductDTensor%rank=DimSize(newD)
				ProductDTensor%flag=.true.
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
					write(*,*)"ERROR in ProductDTensor,case 3,stop"
					stop
				end if
				allocate(ProductDTensor%DTensor_Data(T1m))
				call DGEMV('N',T1m,T1n,1D0,T1%DTensor_Data,T1m,T2%DTensor_Data,&
						1,0D0,ProductDTensor%DTensor_Data,1)
				ProductDTensor%totalData=T1m
				ProductDTensor%TenDim=newD
				ProductDTensor%rank=DimSize(newD)
				ProductDTensor%flag=.true.
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
					write(*,*)"error ProductDTensor,dimension"
					stop
				end if
				T1m=D1.i.1
				T1n=D1.i.2
				T2m=D2.i.1
				T2n=D2.i.2
				totaldata=T1m*T2n
				allocate(ProductDTensor%DTensor_data(totalData))
				call DGEMM('N','N',T1m,T2n,T1n,1D0,T1%DTensor_Data,T1m,&
					T2%DTensor_Data,T2m,0D0,ProductDTensor%DTensor_data,T1m)
				ProductDTensor%totalData=totaldata
				ProductDTensor%TenDim=newD
				ProductDTensor%rank=DimSize(newD)
				ProductDTensor%flag=.true.
		case default 
			write(*,*) "ERROR in ProductDTensor,no such data"
			stop
		end 	select
		return
	end function
!************************  operateDTensor  ****************************************
!	Operar is a matrix,the first DElement of every row specify
!		0: do nothing
!		1:Dcontracts
!		2:Ddecompose
!		3:DdecomposeAll
!	other DElement of every row is input parameter of the function
!		at last A*B .when flag=1,result store in A with B no change when productflag=0,
!	 the result will store in B
!		this subroutine run for big inoutA and B,do as less as possible operation on the
!	 DTensor_Data of inoutA and B.
	subroutine operateDTensor2(A,B,oA,oB,productflag)
		type(DTensor),intent(inout)::A
		type(DTensor),intent(inout)::B
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
			write(*,*)"ERROR in operateDTensor,productflag",productflag
			stop
		end if
		if(oA(1,1).ne.0) then
			call Dcompose_decompse_subroutine2(A%TenDim,A%rank,oA)
		end if
		if(oB(1,1).ne.0) then
			call Dcompose_decompse_subroutine2(B%TenDim,B%rank,oB)
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
			write(*,*)"ERROR in operateDTensor"
		end if
		return
	end subroutine
	subroutine operateDTensor1(A,B,oA,oB,productflag)
		type(DTensor),intent(inout)::A
		type(DTensor),intent(inout)::B
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
			write(*,*)"ERROR in operateDTensor,productflag",productflag
			stop
		end if
		if(oA(1).ne.0) then
			call Dcompose_decompse_subroutine1(A%TenDim,A%rank,oA)
		end if
		if(oB(1).ne.0) then
			call Dcompose_decompse_subroutine1(B%TenDim,B%rank,oB)
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
			write(*,*)"ERROR in operateDTensor"
		end if
		return
	end subroutine
		
		
		
		
		
		
		
		
		
!**************** Dcontract   ****************
!		combine two index of the DTensor,which is con_index and con_index+1
	type(DTensor) function Dcontract(T1,con_index)
		integer,intent(in) :: con_index
		type(DTensor),intent(in) :: T1
		type(Dimension) ::newdim
		newdim=DimConstract(T1%TenDim,con_index,1)
		call DstoreTenData(Dcontract,T1%DTensor_Data)
		Dcontract%rank=DimSize(newDim)
		Dcontract%totalData=T1%totalData
		Dcontract%TenDim=newdim
		Dcontract%flag=.true.
		return
	end function
!**************** Dcontracts   ****************
!		combine two index of the DTensor,which is index1 index1+1,..,index2-1,index2
!		if index2 larger than rank,the function will Dcontract the index of index1 -> rank
	type(DTensor) function Dcontracts(T1,index1,index2)
		integer,intent(in) :: index1,index2
		type(DTensor),intent(in) :: T1
		type(Dimension) ::newdim
		integer ::num
		num=index2-index1
		newdim=DimConstract(T1%TenDim,index1,num)
		call DstoreTenData(Dcontracts,T1%DTensor_Data)
		Dcontracts%rank=DimSize(newDim)
		Dcontracts%totalData=T1%totalData
		Dcontracts%TenDim=newdim
		Dcontracts%flag=.true.
		return
	end function			
	type(DTensor) function Dcontractsv(T1,vector)
		integer,intent(in) ::vector(2)
		type(DTensor),intent(in) :: T1
		type(Dimension) ::newdim
		integer ::num
		num=vector(2)-vector(1)
		newdim=DimConstract(T1%TenDim,vector(1),num)
		call DstoreTenData(Dcontractsv,T1%DTensor_Data)
		Dcontractsv%rank=DimSize(newDim)
		Dcontractsv%totalData=T1%totalData
		Dcontractsv%TenDim=newdim
		Dcontractsv%flag=.true.
		return
	end function		
!*****************  Ddecompose  *****************
! Ddecompose the de_index index of the DTensor into n(1),n(2)
!		for example the de_index index is (1,2,3,4,..inde,inde+1,...rank)
!		(1,2,3,4,..inde,inde+1,...rank)-->(1,2,3,4,..inde),(inde+1,...rank)		
!		if inde larger than rank ,the function will return no change	
	type(DTensor) function Ddecompose(T1,de_index,inde)
		type(DTensor),intent(in) :: T1
		integer,intent(in) :: de_index,inde
		type(Dimension) :: newDim
		newDim=Dimdecompose(T1%TenDim,de_index,inde)
		call DstoreTenData(Ddecompose,T1%DTensor_Data)
		Ddecompose%rank=DimSize(newDim)
		Ddecompose%totalData=T1%totalData
		Ddecompose%TenDim=newdim
		Ddecompose%flag=.true.
		return
	end function
	type(DTensor) function Ddecomposev(T1,vector)
		type(DTensor),intent(in) :: T1
		integer,intent(in) ::vector(2)
		integer:: de_index,inde
		type(Dimension) :: newDim
		de_index=vector(1)
		inde=vector(2)
		newDim=Dimdecompose(T1%TenDim,de_index,inde)
		call DstoreTenData(Ddecomposev,T1%DTensor_Data)
		Ddecomposev%rank=DimSize(newDim)
		Ddecomposev%totalData=T1%totalData
		Ddecomposev%TenDim=newdim
		Ddecomposev%flag=.true.
		return
	end function
			
	type(DTensor) function Ddecompose1(T1,de_index)
		type(DTensor),intent(in) :: T1
		integer,intent(in) :: de_index
		type(Dimension) :: newDim
		newDim=Dimdecompose(T1%TenDim,de_index,1)
		call DstoreTenData(Ddecompose1,T1%DTensor_Data)
		Ddecompose1%rank=DimSize(newDim)
		Ddecompose1%totalData=T1%totalData
		Ddecompose1%TenDim=newdim
		Ddecompose1%flag=.true.
		return
	end function
			
	type(DTensor) function DdecomposeAll(T1)
		type(DTensor),intent(in) :: T1
		integer:: i
		type(Dimension) :: newDim
		newDim=DimdecomposeAll(T1%TenDim)
		call DstoreTenData(DdecomposeAll,T1%DTensor_Data)
		DdecomposeAll%rank=DimSize(newDim)
		DdecomposeAll%totalData=T1%totalData
		DdecomposeAll%TenDim=newdim
		DdecomposeAll%flag=.true.
		return
	end function
!**********************  Dcompose_decompse  *********************
!		Dcompose_decompse do the Dcontract or Ddecompose on the DTensor
!		Operar is a matrix,the first DElement of every row specify
!			1:Dcontracts
!			2:Ddecompose
!			3:DdecomposeAll
!		other DElement of every row is input parameter of the function
	type(DTensor) function Dcompose_decompse(T1,Operar)
		integer,intent(in)::Operar(:,:)
		type(DTensor),intent(in) :: T1
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
				dimen=Dimdecompose(dimen,de_index,inde)
			else if(Operar(i,1).eq.3) then
				dimen=DimdecomposeAll(dimen)
			else
				write(*,*) "error in Dcompose_decompse"
				stop
			end if
		end do
		call DstoreTenData(Dcompose_decompse,T1%DTensor_Data)
		Dcompose_decompse%TenDim=dimen
		Dcompose_decompse%rank=DimSize(dimen)
		Dcompose_decompse%flag=.true.
		Dcompose_decompse%totalData=T1%totalData
		return
	end function
!********************  DdimOperations  ****************************
!		 DdimOperations do the Dcontract or Ddecompose on the DTensor with
!	  no change on the Tensot_data in the DTensor.
!		Operar is a matrix,the first DElement of every row specify
!		 1:Dcontracts
!		 2:Ddecompose
!		 3:DdecomposeAll
!		other DElement of every row is input parameter of the function
	subroutine DdimOperations(T1,Operar)
		integer,intent(in)::Operar(:,:)
		type(DTensor),intent(inout) :: T1
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
				T1%TenDim=Dimdecompose(T1%TenDim,de_index,inde)
				T1%rank=DimSize(T1%TenDim)
			else if(Operar(i,1).eq.3) then
				T1%TenDim=DimdecomposeAll(T1%TenDim)
				T1%rank=DimSize(T1%TenDim)
			else
				write(*,*) "error in Dcompose_decompse"
				stop
			end if
		end do
		return
	end subroutine
	subroutine DdimOperation1(T1,Operat)
		integer,intent(in)::Operat(:)
		type(DTensor),intent(inout) :: T1
		integer,allocatable::Operat2(:,:)
		allocate(Operat2(1,size(Operat)))
		Operat2(1,:)=Operat
		call DdimOperations(T1,Operat2)
		return
	end subroutine
	subroutine Dcompose_decompse_subroutine2(TenDim,rank,Operar)
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
				TenDim=Dimdecompose(TenDim,de_index,inde)
				rank=DimSize(TenDim)
			else if(Operar(i,1).eq.3) then
				TenDim=DimdecomposeAll(TenDim)
				rank=DimSize(TenDim)
			else
				write(*,*) "error in Dcompose_decompse_subroutine"
				stop
			end if
		end do
		return
	end subroutine
	subroutine Dcompose_decompse_subroutine1(TenDim,rank,Operar)
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
			TenDim=Dimdecompose(TenDim,de_index,inde)
			rank=DimSize(TenDim)
		else if(Operar(1).eq.3) then
			TenDim=DimdecomposeAll(TenDim)
			rank=DimSize(TenDim)
		else
			write(*,*) "error in Dcompose_decompse_subroutine"
			stop
		end if
		return
	end subroutine

	
!***************** equal_of_DTensor *****************
	logical function equal_of_DTensor(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		integer :: l
		equal_of_DTensor=.true.
		if(T1%rank.ne.T2%rank) then
			equal_of_DTensor=.false.
			return
		end if
		if(.not.equal_of_dim(T1%TenDim,T2%TenDim)) then
			equal_of_DTensor=.false.
			return
		end if
		l=count(abs(T1%DTensor_data-T2%DTensor_data).gt.1d-10)
		if(l.gt.0) then
			equal_of_DTensor=.false.
			return
		end if
		return
	end function
!*****************  maxDElement  *****************
	real*8 function dmaxabsElement(T)
		type(DTensor),intent(in) :: T
		dmaxabsElement=maxval(abs(T%DTensor_Data))
		return
	end function
!*****************  maxRealDElement  *****************
	real*8 function dmaxElement(T)
		type(DTensor),intent(in) :: T
		dmaxElement=maxval(T%DTensor_Data)
		return
	end function
!*****************  minDElement  *****************
	real*8 function dminabsElement(T)
		type(DTensor),intent(in) :: T
		dminabsElement=minval(abs(T%DTensor_Data))
		return
	end function
!*****************  minRealDElement  *****************
	real*8 function dminElement(T)
		type(DTensor),intent(in) :: T
		dminElement=minval(T%DTensor_Data)
		return
	end function
!ccccccccccccccccc svd cccccccccccccccccc	

!*****************  Dsvdright  *****************
! 			T should be two dimension		
!			T=U*s*V^T
!			return the (U*s) in T1
!			and V^T in T2
!			the data s(i) in the matrix s will be discard if
!			s(i)<s(1)*check_discard_svd.
!			if check_discard_svd<=0,all the data in s will be keep
!			note:I do not write the faster version of this part
	subroutine Dsvdright(T,T1,T2,check_discard)
		type(DTensor),intent(in) :: T
		type(DTensor),intent(out) :: T1,T2
		real*8,intent(in) :: check_discard
		real*8,allocatable :: T1data(:,:)
		real*8,allocatable :: T2data(:,:)
		real*8,allocatable :: newT2data(:,:)
		real*8,allocatable :: newT1data(:,:)
		real*8,allocatable :: Tdata(:,:)
		real*8,allocatable :: work(:)
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
				write(*,*)"error in Dsvdright"
			end if
			allocate(work(lw))
			Tdata=reshape(T%DTensor_Data,(/m,n/))
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
				T1=Dset(T1Dim,reshape(T1data,(/m*m/)))
				T2=Dset(T2Dim,reshape(newT2data,(/n*ms_max/)))
			else
				allocate(newT1data(m,ms_max))
				newT1data=T1data(:,1:ms_max)
				T1Dim=(/ms_max/)
				T1Dim=(T%TenDim.sub.1)+T1Dim
				T2Dim=(/ms_max/)
				T2Dim=T2Dim+(T%TenDim.sub.2)
				T1=Dset(T1Dim,reshape(newT1data,(/m*ms_max/)))
				T2=Dset(T2Dim,reshape(newT2data,(/n*ms_max/)))
			end if
		else
			write(*,*) "Input DTensor should be 2 dimension"
		end if
		return
	end subroutine												
!cccccccccccccccc Dsvdleft cccccccccccccccccc	
! 				T should be two dimension		
!				T=U*s*V^T
!				return the transpose U in T1
!			and (s*V^T) in T2
!				the data s(i) in the matrix s will be discard if
!			s(i)<s(1)*check_discard_svd.
!				if check_discard_svd<0,all the data in s will be keep
!			wornning:I do not check the part check_discard>0
	subroutine Dsvdleft(T,T1,T2,check_discard)
		type(DTensor),intent(in) :: T
		type(DTensor),intent(out) :: T1,T2
		real*8,intent(in) :: check_discard
		real*8,allocatable :: T1data(:)
		real*8,allocatable :: T2data(:)
		real*8,allocatable :: newT2data(:)
		real*8,allocatable :: Tdata(:)
		real*8,allocatable :: work(:)
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
				write(*,*)"error in Dsvdleft"
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
				T2=Dset(T2Dim,T2data)
				T1=Dset(T1Dim,T1data)
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
				T2=Dset(T2Dim,newT2data)
				T1=Dset(T1Dim,T1data)
			end if
		else
			write(*,*) "Input DTensor should be 2 dimension"
		end if
		return
	end subroutine
!Solves a general system of linear equations
!A*X=B,find X
!X is a vector or a matrix,the dimension of which is the same as B
	type(DTensor) function Dlinequ(A,B)
		type(DTensor),intent(in)::A,B
		real*8,allocatable::Xdata(:),Adata(:)
		type(dimension)::Xdim
		integer,allocatable::IPIV(:)
		integer::Na,Nb,INFO
		Na=A.dim.1
		if(Na.ne.(A.dim.2)) then
			write(*,*)"error in Dlinequ,dimension of A"
			stop
		end if
		if(DgetRank(B).eq.1)then
			Nb=1
		else
			Nb=B.dim.2
		end if
		if(Na.ne.(B.dim.1)) then
			write(*,*)"error in Dlinequ,dimension of A and B"
			stop
		end if
		allocate(IPIV(Na))
		Xdata=B
		Adata=A!the subroutine ZGESV will change A
		call DGESV(Na,Nb,Adata,Na,IPIV,Xdata,Na,INFO)
		if(INFO.ne.0)then
			write(*,*)"ZGESV is not successful"
			write(*,*)"INFO",INFO
			stop
		end if
		
		Xdim=(A%TenDim.sub.2)
		if(DgetRank(B).ne.1)then
			Xdim=Xdim+(B%TenDim.sub.2)
		end if
		Dlinequ=Dset(Xdim,Xdata)
		return
	end function
!Solves a general system of linear equations
!A*X=B,find X
!X is a vector or a matrix,the dimension of which is the same as B
!on output,A will change and X is in B	
	subroutine Dlinequ_routine(A,B)
		type(DTensor),intent(inout)::A,B
		integer,allocatable::IPIV(:)
		integer::Na,Nb,INFO
		type(dimension)::Xdim
		Na=A.dim.1
		if(Na.ne.(A.dim.2)) then
			write(*,*)"error in Dlinequ,dimension of A"
			stop
		end if
		if(DgetRank(B).eq.1)then
			Nb=1
		else
			Nb=B.dim.2
		end if
		if(Na.ne.(B.dim.1)) then
			write(*,*)"error in Dlinequ,dimension of A and B"
			stop
		end if
		allocate(IPIV(Na))
		call DGESV(Na,Nb,A%DTensor_Data,Na,IPIV,B%DTensor_Data,Na,INFO)
		if(INFO.ne.0)then
			write(*,*)"ZGESV is not successful"
			write(*,*)"INFO",INFO
			stop
		end if
		Xdim=(A%TenDim.sub.2)
		if(DgetRank(B).ne.1)then
			Xdim=Xdim+(B%TenDim.sub.2)
		end if
		call Dresetdim(B,Xdim)
		return
	end subroutine
	
!*****************  Deye  *****************
	type(DTensor) function Deye_Ten1(T)
		type(DTensor),intent(in) :: T
		real*8,allocatable :: sdata(:,:)
		integer :: i,j,lens,n
		if(T%rank.ne.1) then
			write(*,*)"ERROR in Deye_Ten,input should be a vector"
			stop
		end if
		n=T%TotalData
		allocate(sdata(n,n))
		sdata=0d0
		do i=1,n
			sdata(i,i)=T%DTensor_Data(i)
		end do
		call DstoreTenData(Deye_Ten1,sdata)
		Deye_Ten1%rank=2
		Deye_Ten1%totalData=n*n
		Deye_Ten1%TenDim=(/n,n/)
		Deye_Ten1%flag=.true.
		return
	end function
	type(DTensor) function Deye_Ten2(T,m,n)
		type(DTensor),intent(in) :: T
		integer,intent(in) :: m,n
		real*8,allocatable :: sdata(:,:)
		integer :: i,j,lens,k
		if(T%rank.ne.1) then
			write(*,*)"ERROR in Deye_Ten,input should be a vector"
			stop
		end if
		lens=T%TotalData
		allocate(sdata(m,n))
		sdata=0d0
		do i=1,min(m,n)
			if(lens.gt.i) then
				sdata(i,i)=T%DTensor_Data(i)
			end if
		end do
		call DstoreTenData(Deye_Ten2,sdata)
		Deye_Ten2%rank=2
		Deye_Ten2%totalData=m*n
		Deye_Ten2%TenDim=(/m,n/)
		Deye_Ten2%flag=.true.
		return
	end function
!*****************  Deye_real  *****************
	type(DTensor) function Deye_real(s,m,n)
		real*8,intent(in) :: s(:)
		integer,intent(in) :: m,n
		real*8,allocatable :: sdata(:,:)
		integer :: i,j,lens
		lens=size(s)
		allocate(sdata(m,n))
		sdata=0d0
		do i=1,min(m,n)
			if(lens.gt.i) then
				sdata(i,i)=s(i)
			end if
		end do
		call DstoreTenData(Deye_real,sdata)
		Deye_real%rank=2
		Deye_real%totalData=m*n
		Deye_real%TenDim=(/m,n/)
		Deye_real%flag=.true.
		return
	end function
!*****************  one_Ten  *****************
	type(DTensor) function Done(m,n)
		integer,intent(in) :: m,n
		real*8,allocatable :: sdata(:,:)
		integer :: i,j
		allocate(sdata(m,n))
		sdata=0d0
		do i=1,min(m,n)
				sdata(i,i)=1d0
		end do
		call DstoreTenData(Done,sdata)
		Done%rank=2
		Done%totalData=m*n
		Done%TenDim=(/m,n/)
		Done%flag=.true.
		return
	end function	
!****************  Done_Ten  *************
	type(DTensor) function Done_Ten(dimen)
		integer,intent(in) :: dimen(:)
		integer :: i,total,lenDim,min_dim,address,dim_i
		real*8,allocatable::TenData(:)
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
			TenData(address)=1d0
		end do
		Done_Ten=DbuildTen(dimen,TenData)
		return
	end function	
!*****************  Dtrace  *********************
! return Dtrace(A),A is a matrix

	real*8 function Dtrace(T)
		type(DTensor),intent(in) :: T
		integer::rank,i
		rank=DgetRank(T)
		if(rank.ne.2) then
			write(*,*)"error in Dtrace"
			write(*,*)"input DTensor should be a matrix"
			stop
		end if
		if((T.dim.1).ne.(T.dim.2)) then
			write(*,*)"error in Dtrace"
			write(*,*)"input DTensor should be a matrix"
			write(*,*)(T.dim.1),(T.dim.2)
			stop
		end if
		Dtrace=dcmplx(0d0,0d0)
		do i=1,(T.dim.1)
			Dtrace=Dtrace+(T.i.(/i,i/))
		end do
		return
	end function	
		
		
!*****************  Dvalue  *****************
!			return the Dvalue of <phi|F|phi>,where F is an operator
!			F=<I',s1',s2',..|F|I,s1,s2,..>
!			update at 2013.10.29
	real*8 function Dvalue(F_,Tr_)
		type(DTensor),intent(in) :: F_,Tr_
		type(DTensor) :: Tr,Tl,F,midResult
		integer :: i,rank,oper(2,3)
		rank=Tr_%rank
		Tr=Tr_.c.(/1,rank/)
		oper(1,:)=(/1,1,rank/)
		oper(2,:)=(/1,2,rank+1/)
		F=F_.cd.oper
		Dvalue=Tl*F*Tr
		return
	end function		
!***************   dot   ******************
!		return <phi1|phi2>
!	
	real*8 function Doudot(phi1,phi2)
		Type(DTensor),intent(in)::phi1,phi2
		integer::N1,N2
		N1=DgetTotalData(phi1)
		N2=DgetTotalData(phi2)
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			stop
		end if
		Doudot=Ddot(N1,phi1%dTensor_Data,1,phi2%dTensor_Data,1)
		RETURN
	end function	

!*****************  Dnorm   *****************
!			return  <phi|phi>
!			update at 2013.10.29
	real*8 function Dnorm2(Tr)
		type(DTensor),intent(in) :: Tr
		Dnorm2=dnrm2(Tr%TotalData,Tr%DTensor_Data,1)
		Dnorm2=Dnorm2*Dnorm2
		return
	end function	
!		return  sqrt(<phi|phi>	)
	real*8 function Dnorm(Tr)
		type(DTensor),intent(in) :: Tr
		Dnorm=dnrm2(Tr%TotalData,Tr%DTensor_Data,1)
		return
	end function		
!*****************  addressToIndes   *****************
	integer function addressToIndes(T,Adim)
		type(DTensor),intent(in) :: T
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
		
! 		Dmodify the code below in 2013.5.18		
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
! 		Dmodify the code below in 2013.5.18
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
!*****************  DElement   *****************
	real*8 function DElement(T,Tdim)
		type(DTensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		integer::Inde
		inde=addressToIndes(T,Tdim)
		DElement=T%DTensor_Data(inde)
		return
	end function		  
	real*8 function DElement2(T,inde)
		type(DTensor),intent(in) ::T
		integer,intent(in)::inde
		DElement2=T%DTensor_Data(inde)
		return
	end function
!		inde=[-1,inde_min,inde_max] output data(:,inde_min:inde_max)
!	or[-2,inde_min,inde_max],data(inde_min:inde_max,:)
!	or [-1,inde_row] [-2,inde_col],output row or col

	type(DTensor) function DsubTen(T,inde)	
		type(DTensor),intent(in) ::T
		integer,intent(in)::inde(:)
		real*8,allocatable::Tdata(:,:),DsubTen1(:),DsubTen2(:,:)
		integer::num
		if(DgetRank(T).ne.2)then
			write(*,*)"error in DsubTen,only matrix is allowed"
			stop
		end if
		Tdata=T
		if(size(inde).eq.2) then
			select case (inde(1))
				case (-1)!output row
					allocate(DsubTen1(T.dim.2))
					DsubTen1=Tdata(inde(2),:)
				case (-2)!
					allocate(DsubTen1(T.dim.1))
					DsubTen1=Tdata(:,inde(2))
				case default 
					write(*,*) "no such case in DsubTen"
					write(*,*)inde
					stop
			end 	select
			DsubTen=DsubTen1
			return
		end if
		if(size(inde).ne.3) then
			write(*,*) "no such case in DsubTen"
			write(*,*) "length of inde is ",size(inde)
			write(*,*)inde
			stop
		end if
		num=inde(3)-inde(2)+1
		if(num.le.0) then
			write(*,*)"error DsubTen"
			stop
		end if
		select case (inde(1))
			case (-1)!output row
				allocate(DsubTen2(num,T.dim.2))
				DsubTen2=Tdata(inde(2):inde(3),:)
			case (-2)!
				allocate(DsubTen2(T.dim.1,num))
				DsubTen2=Tdata(:,inde(2):inde(3))
			case default 
				write(*,*) "no such case in DsubTen"
				write(*,*)inde
				stop
		end 	select
		DsubTen=DsubTen2
		return
	end function
!***************  DdirectProduct  *********************
	type(DTensor) function DdirectProduct(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: T1data(:,:),T2data(:,:)
		real*8,allocatable :: newdata(:,:,:,:)
		integer :: m1,n1,m2,n2,i,j,k,l,rank(2)
		type(Dimension):: D1,Dtemp
		if((DgetRank(T1).eq.1).and.(DgetRank(T2).eq.1)) then
			DdirectProduct=DdirectProduct1(T1,T2)
			return
		end if
		if((DgetRank(T1).eq.2).and.(DgetRank(T2).eq.2)) then
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
			T1data=reshape(T1%DTensor_Data,(/m1,n1/))
			T2data=reshape(T2%DTensor_Data,(/m2,n2/))
			do l=1,n2
				do k=1,n1
					 do j=1,m2
						 do i=1,m1
							newdata(i,j,k,l)=T1data(i,k)*T2data(j,l)
						end do
					end do
				end do
			end do
			DdirectProduct=Dset(D1,reshape(newdata,(/m1*m2*n1*n2/)))
		else
			write(*,*)"The DdirectProduct is just for matrix"
			stop
		end if
	return
	end function
!********************  DdirectProductM  ************************
	type(DTensor) function DdirectProductM(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: T1data(:,:),T2data(:,:)
		real*8,allocatable :: newdata(:,:,:,:)
		integer :: m1,n1,m2,n2,i,j,k,l
		type(Dimension):: D1,Dtemp
		if((DgetRank(T1).eq.1).and.(DgetRank(T2).eq.1)) then
			DdirectProductM=DdirectProduct1V(T1,T2)
			return
		end if
		if((DgetRank(T1).eq.2).and.(DgetRank(T2).eq.2)) then
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
			T1data=reshape(T1%DTensor_Data,(/m1,n1/))
			T2data=reshape(T2%DTensor_Data,(/m2,n2/))
			do l=1,n2
				do k=1,n1
					do j=1,m2
						do i=1,m1
							newdata(i,j,k,l)=T1data(i,k)*T2data(j,l)
						end do
					end do
				end do
			end do
			DdirectProductM=Dset(D1,reshape(newdata,(/m1*m2*n1*n2/)))
		else
			write(*,*)"The DdirectProductM is just for matrix"
			stop
		end if
		return
	end function
!************* DdirectProduct1 **************************
!
!       for two DTensor whose ranks are 1		
!
	type(DTensor) function DdirectProduct1V(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: newdata(:,:)
		integer :: m1,n1,m2,n2,i,j,k,l
		type(Dimension):: D1,Dtemp
		if((DgetRank(T1).eq.1).and.(DgetRank(T2).eq.1)) then
			D1=T1%TenDim
			D1=D1+T2%TenDim
			D1=DimConstract(D1,1,1)
			n1=T1.dim.1
			n2=T2.dim.1
			allocate(newdata(n1,n2))
			do l=1,n2
				do k=1,n1
					newdata(k,l)=T1%DTensor_Data(k)*T2%DTensor_Data(l)
				end do
			end do
			DdirectProduct1V=Dset(D1,reshape(newdata,(/n1*n2/)))
		else
			write(*,*)"The DdirectProduct1 is just for vector"
			stop
		end if
		return
	end function	
	type(DTensor) function DdirectProduct1(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: newdata(:,:)
		integer :: m1,n1,m2,n2,i,j,k,l
		type(Dimension):: D1,Dtemp
		if((DgetRank(T1).eq.1).and.(DgetRank(T2).eq.1)) then
			D1=T1%TenDim
			D1=D1+T2%TenDim
			n1=T1.dim.1
			n2=T2.dim.1
			allocate(newdata(n1,n2))
			do l=1,n2
				do k=1,n1
					newdata(k,l)=T1%DTensor_Data(k)*T2%DTensor_Data(l)
				end do
			end do
			DdirectProduct1=Dset(D1,reshape(newdata,(/n1*n2/)))
		else
			write(*,*)"The DdirectProduct1 is just for vector"
			stop
		end if
		return
	end function
!*******************  direct sum  **********************
	type(DTensor) function  Ddirect_sum(T1,T2) result(T)
		type(DTensor),intent(in) :: T1,T2
		type(DTensor) ::temp1,temp2
		type(Dimension) ::dimT1,dimT2
		integer :: m,n
		dimT1=T1%TenDim
		dimT2=T2%TenDim
		if(dimT1.equ.dimT2) then
			temp1=T1.mxx.deye(m,n)
			temp2=deye(m,n).mxx.T2
			T=temp1+temp2
			T=.dc.T
		else
			write(*,*) "error in Ddirect_sumM"
		end if
		return 
	end function
!****************  direct sum returning matrix *******************
	type(DTensor) function  Ddirect_sumM(T1,T2) result(T)
		type(DTensor),intent(in) :: T1,T2
		type(DTensor) ::temp1,temp2
		type(Dimension) ::dimT1,dimT2
		integer :: m,n
		dimT1=T1%TenDim
		dimT2=T2%TenDim
		if(dimT1.equ.dimT2) then
			m=T1.dim.1
			n=T1.dim.2
			temp1=T1.mxx.deye(m,n)
			temp2=deye(m,n).mxx.T2
			T=temp1+temp2
		else
			write(*,*) "error in Ddirect_sumM"
		end if
		return 
	end function
!**************   Dconbine   ********************
!
!                    / T1 \
!	Dconbine(T1,T2)	=  |----|
!		               \ T2 /
!	T1 :a [...,m,n,l] matrix
!	T2 :a [...,m,n,l] matrix
!	Dconbine(T1,T2):a [...,m,n,l,2] matrix
!or 
!	T1 :a [...,m,n,l] matrix
!	T2 :a [...,m,n] matrix
!	Dconbine(T1,T2):a [...,m,n,l+1] matrix
	type(DTensor) function Dconbine(T1,T2)
		type(DTensor),intent(in)::T1,T2
		integer,allocatable::dim1(:),dim2(:)
		integer::total,total1,i,dim_n
		real*8,allocatable::newdata(:),olddata(:,:)
		type(Dimension)::newDim,dimen1,dimen2
		if(T1%rank.eq.T2%rank) then
			dim1=.subDim.T1
			dim2=.subDim.T2
			if(.not.(dim1.equ.dim2)) then
				write(*,*)"can not conbie two DTensor"
				write(*,*)dim1
				write(*,*)dim2
				write(*,*)"stop"
				stop
				return
			end if
			total1=DgetTotalData(T1)
			total=total1+total1
			allocate(newdata(total1*2))
			newdata(:total1)=T1%DTensor_Data
			newdata(total1+1:)=T2%DTensor_Data
			call DstoreTenData(Dconbine,newdata)
			Dconbine%Rank=DgetRank(T1)+1
			newDim=.subDim.T1
			newDim=newDim+(/2/)
			Dconbine%TenDim=newDim
			Dconbine%totalData=total
			Dconbine%flag=.true.
			return
		end if
		dimen1=.subDim.T1
		dimen1=DimConstract(dimen1,1,T1%rank-2)
		dimen2=.subDim.T2
		dimen2=DimConstract(dimen2,1,T2%rank)
		if(.not.((dimen1.sub.1).equ.dimen2)) then
			write(*,*)"can not conbie two DTensor"
			write(*,*)"stop2"
			stop
			return
		end if
		dim_n=dimen1.i.2
		total1=DgetTotalData(T1)
		total=DgetTotalData(T2)
		allocate(newdata(total+total1))
		newdata(1:total1)=T1%DTensor_Data!olddata
		newdata(total1+1:)=T2%DTensor_Data
		call DstoreTenData(Dconbine,newdata)
		Dconbine%Rank=DgetRank(T1)
		newDim=(.subDim.T2)
		newDim=newDim+(/dim_n+1/)
		Dconbine%TenDim=newDim
		Dconbine%totalData=total+total1
		Dconbine%flag=.true.
	end function
	type(DTensor) function Dconbinerow(T1,T2)
		type(DTensor),intent(in)::T1,T2
		integer,allocatable::dim1(:),dim2(:)
		integer::total,total1,i,dim_n
		real*8,allocatable::newdata(:,:),olddata(:,:)
		type(Dimension)::newDim,dimen1,dimen2
		if(T1%rank.eq.T2%rank) then
			dim1=.subDim.T1
			dim2=.subDim.T2
			if(.not.(dim1.equ.dim2)) then
				write(*,*)"can not Dconbinerow two DTensor"
				write(*,*)dim1
				write(*,*)dim2
				write(*,*)"stop"
				stop
				return
			end if
			total1=DgetTotalData(T1)
			total=total1+total1
			allocate(newdata(2,total1))
			newdata(1,:)=T1%DTensor_Data
			newdata(2,:)=T2%DTensor_Data
			call DstoreTenData(Dconbinerow,newdata)
			Dconbinerow%Rank=DgetRank(T1)+1
			newDim=.subDim.T1
			newDim=(/2/)+newDim
			Dconbinerow%TenDim=newDim
			Dconbinerow%totalData=total
			Dconbinerow%flag=.true.
			return
		end if
		dimen1=.subDim.T1
		dimen1=DimConstract(dimen1,2,T1%rank)
		dimen2=.subDim.T2
		dimen2=DimConstract(dimen2,1,T2%rank)
		if(.not.((dimen1.sub.2).equ.dimen2)) then
			write(*,*)"can not Dconbinerow two DTensor"
		!	write(*,*)T1%rank,T2%rank
			write(*,*)"stop2"
			stop
			return
		end if
		dim_n=dimen1.i.1
		total1=DgetTotalData(T1)
		total=DgetTotalData(T2)
		allocate(newdata(dim_n+1,total))
		allocate(olddata(dimen1.i.1,dimen1.i.2))
		call dcopy(DgetTotalData(T1),T1%DTensor_Data,1,olddata,1)
		newdata(1:dim_n,:)=olddata
		newdata(dim_n+1,:)=T2%DTensor_Data
		call DstoreTenData(Dconbinerow,newdata)
		Dconbinerow%Rank=DgetRank(T1)
		newDim=(.subDim.T2)
		newDim=(/dim_n+1/)+newDim
		Dconbinerow%TenDim=newDim
		Dconbinerow%totalData=total+total1
		Dconbinerow%flag=.true.
	end function			
				
!****************  DTensor1  ************************
!	dimen is the dimension of the DTensor
!	return the DTensor with all the data are num
!	num is a real*8 type
	type(DTensor)function DTensor1(dimen,num)
		type(Dimension),intent(in) :: dimen
		real*8,intent(in) :: num
		integer::total,i
		total=1
		do i=1,DimSize(dimen)
			total=total*(dimen.i.i)
		end do
		allocate(DTensor1%DTensor_data(total))
		do i=1,total
			DTensor1%DTensor_data(i)=num
		end do
		DTensor1%rank=DimSize(dimen)
		DTensor1%totalData=total
		DTensor1%TenDim=dimen
		DTensor1%flag=.true.
		return
	end function
!******************  Dmodify the data in the DTensor  **************
! Dmodify the data in dimen of the DTensor
	subroutine DmodifyTen_val(Ten,dimen,val)
		type(DTensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		real*8,intent(in)::val
		integer::addre
		addre=addressToIndes(Ten,dimen)
		Ten%DTensor_Data(addre)=val
		return
	end subroutine
	subroutine DmodifyTen_val2(Ten,inde,val)
		type(DTensor),intent(inout)::Ten
		integer,intent(in)::inde
		real*8,intent(in)::val
		integer::totalData
		real*8,allocatable::Tdata(:)
		totalData=Ten%totalData
		if(inde.le.totalData) then
			Ten%DTensor_Data(inde)=val
		else
			if(Ten%rank.gt.1)then
				write(*,*)"the index is larger than the length of the Tensor"
				write(*,*)"Only rank=1 is allow for re allocate"
				write(*,*)"stop"
				stop
			end if
			Tdata=Ten
			if(Ten%flag)then
				deallocate (Ten%DTensor_Data)
				allocate(Ten%DTensor_Data(inde))
				Ten%DTensor_Data(1:totalData)=Tdata
				Ten%DTensor_Data(inde)=val
				Ten%totalData=inde
				Ten%TenDim=(/inde/)
			else
				allocate(Ten%DTensor_Data(inde))
				Ten%DTensor_Data(inde)=val
				Ten%totalData=inde
				Ten%TenDim=(/inde/)
				Ten%flag=.true.
			end if
		end if
		return
	end subroutine
! Dmodify all the data in the DTensor
	subroutine DmodifyTen_dat(Ten,val)
		type(DTensor),intent(inout)::Ten
		real*8,intent(in)::val(:)
		if(size(val).eq.Ten%TotalData) then
			write(*,*)"ERROR in Dmodify DTensor data"
			write(*,*)"stop"
			stop
		end if
		call dcopy(size(val),val,1,Ten%DTensor_Data,1)
		return
	end subroutine
! Dmodify all the data in the DTensor
	subroutine DmodifyTen_dat2(Ten,val)
		type(DTensor),intent(inout)::Ten
		real*8,intent(in)::val(:,:)
		if(Ten%rank.ne.2) then
			write(*,*)"ERROR in Dmodify DTensor data"
			write(*,*)"The rank is 2,stop"
			stop
		end if
		if(product(val).eq.Ten%TotalData) then
			write(*,*)"ERROR in Dmodify DTensor data"
			write(*,*)"stop"
			stop
		end if
		call DstoreTenData(Ten,val)
		return
	end subroutine	
! Dmodify all the data in the DTensor
	subroutine DmodifyTen_dat3(Ten,val)
		type(DTensor),intent(inout)::Ten
		real*8,intent(in)::val(:,:,:)
		if(Ten%rank.ne.3) then
			write(*,*)"ERROR in Dmodify DTensor data"
			write(*,*)"The rank is 3,stop"
			stop
		end if
		if(product(val).eq.Ten%TotalData) then
			write(*,*)"ERROR in Dmodify DTensor data"
			write(*,*)"stop"
			stop
		end if
		call DstoreTenData(Ten,val)
		return
	end subroutine	
! Dmodify all the dimension in the DTensor
	subroutine Dresetdim1(Ten,dimen)
		type(DTensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		type(Dimension)::dimensi
		if(product(dimen).ne.Ten%TotalData) then
			write(*,*)"ERROR in Dresetdim"
			write(*,*)product(dimen),Ten%TotalData
			write(*,*)"stop"
			stop
		end if
		dimensi=dimen
		Ten%TenDim=dimensi
		Ten%rank=size(dimen)
		return
	end subroutine	
	subroutine Dresetdim2(Ten,dimen)
		type(DTensor),intent(inout)::Ten
		type(dimension),intent(in)::dimen
		integer,allocatable::dimensi(:)
		dimensi=dimen
		if(product(dimensi).ne.Ten%TotalData) then
			write(*,*)"ERROR in Dresetdim"
			write(*,*)product(dimensi),Ten%TotalData
			write(*,*)"stop"
			stop
		end if
		Ten%TenDim=dimen
		Ten%rank=size(dimensi)
		return
	end subroutine	
			
!******************  DTen_Ten  *********************
!		Dcontract The i1 index of T1 and the i2 index of t2
!		T1: (D1,D2,..D_i1,...),T2 :(E1,E2,...,E_i2,....),then the result
!	  will be T:(D1,D2,..D_{i1-1},D_{i1+1},...,D_{rank1},E1,E2,...,E_{i2-1},E_{i2+1},...,E_{rank2})	
!		if if Ddecompose=1, T will be Ddecompose
	type(DTensor) function DTen_Ten(T1_,T2_,i1,i2,ifDdecompose) result(T)
		type(DTensor),intent(in) :: T1_,T2_
		integer,intent(in) :: i1,i2,ifDdecompose
		integer :: i,j,k,D2,decompoD1,decompoD2!Di use for decompese
		type(DTensor) :: T1,T2
		integer :: decompoD1keep,decompoD2keep,rank1,&
			Tdim1(T1_%rank),Tdim2(T2_%rank),rank2,oper(2,3)
		rank1=T1_%rank
		rank2=T2_%rank
		T1=T1_
		T2=T2_
		!T1
		if((rank1.ne.1).and.(rank1.ne.2)) then
			if(i1.eq.1) then
				T1=Dcontracts(T1,2,rank1)
				T1=.p.T1
			else if (i1.eq.rank1) then
				T1=Dcontracts(T1,1,rank1-1)
			else
				oper(1,:)=(/1,1,i1-1/)
				oper(2,:)=(/1,3,rank1/)
				!T1=Dcontracts(T1,1,i1-1)
				!T1=Dcontracts(T1,3,rank1)
				T1=T1.cd.oper
				T1=T1.p.1
				call DdimOperation(T1,(/1,1,2/))
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
				T2=Dcontracts(T2,2,rank2)
			else if (i2.eq.rank2) then
				T2=Dcontracts(T2,1,rank2-1)
				T2=.p.T2
			else
				oper(1,:)=(/1,1,i2-1/)
				oper(2,:)=(/1,3,rank2/)
				T2=T2.cd.oper
				!T2=Dcontracts(T2,1,i2-1)
				!T2=Dcontracts(T2,3,rank2)
				T2=T2.p.3
				call DdimOperation(T2,(/1,2,3/))
				!T2=T2.c.2
			end if
		end if
		if(rank2.eq.2) then
			if(i2.eq.2) then
				T2=.p.T2
			end if
		end if
		T=T1 * T2
		if(ifDdecompose.eq.1) then
			call DdimOperation(T,(/3/))
			!T=.dc.T
		end if
		return
	end function
!*****************Dpermutefo******************************
!		T_{1,2,3,..,i,..,n},Dpermutefo(T,i)=_{i,1,2,3,..,i-1,i+1,..,n}
!
	type(DTensor) function Dpermutefo(T,inde)
		type(DTensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		integer::oper(2,3)
		rank=DgetRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function Dpermutefo"
			write(*,*)"stop"
			stop
		end if
		if(inde.eq.1) then
			Dpermutefo=T
			return
		end if
		if(inde.eq.rank) then
			Dpermutefo=Dcontracts(T,1,rank-1)
			Dpermutefo=.p.Dpermutefo
			Dpermutefo%TenDim=DimdecomposeAll(Dpermutefo%TenDim)
			Dpermutefo%rank=DimSize(Dpermutefo%TenDim)
			return
		end if
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,3,rank/)
		Dpermutefo=T.cd.oper
		Dpermutefo=Dpermutefo.p.3
		Dpermutefo%TenDim=DimdecomposeAll(Dpermutefo%TenDim)
		Dpermutefo%rank=DimSize(Dpermutefo%TenDim)
		return
	end function
!****************** DpermuteInde**************************
!		T_{1,2,3,..,i,..,n},DpermuteInde(T,i)=_{2,3,..,i,1,i+1,..,n}
!
	type(DTensor) function DpermuteInde(T,inde)
		type(DTensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		integer::oper(2,3)
		rank=DgetRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function DpermuteInde",inde,rank
			write(*,*)"stop"
			stop
		end if
		if(inde.eq.1) then
			DpermuteInde=T
			return
		end if
		if(inde.eq.rank) then
			DpermuteInde=Dcontracts(T,2,rank)
			DpermuteInde=.p.DpermuteInde
			DpermuteInde%TenDim=DimdecomposeAll(DpermuteInde%TenDim)
			DpermuteInde%rank=DimSize(DpermuteInde%TenDim)
			return
		end if
		oper(1,:)=(/1,2,inde/)
		oper(2,:)=(/1,3,rank/)
		DpermuteInde=T.cd.oper
		DpermuteInde=DpermuteInde.p.3
		DpermuteInde%TenDim=DimdecomposeAll(DpermuteInde%TenDim)
		DpermuteInde%rank=DimSize(DpermuteInde%TenDim)
		return
	end function



!****************************************************************************************
	 
! 	DTensorlink and DTensornode 	  





!	add a DTensor to the end of the link
!	there is no vector index in the  input DTensor
	subroutine Dpush_backTen(h,T)
		type(DTensorlink),intent(inout) ::h
		type(DTensor),intent(in):: T
		type(DTensornode),pointer ::node
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
!	add a DTensor to the begin of the link	
!	there is no vector index in the  input DTensor
	subroutine Dpush_fo(h,T)
		type(DTensorlink),intent(inout) ::h
		type(DTensor),intent(in):: T
		type(DTensornode),pointer ::node,p
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
!	add a DTensor to the end of the link	
!	inde is the index of the DTensor
	subroutine Dpush_backI(h,T,inde)
		type(DTensorlink),intent(inout) ::h
		type(DTensor),intent(in):: T
		integer,intent(in) :: inde(:)
		type(DTensornode),pointer ::node
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
!	add a DTensor to the begin of the link	
!	inde is the index of the DTensor	
	subroutine Dpush_foI(h,T,inde)
		type(DTensorlink),intent(inout) ::h
		type(DTensor),intent(in):: T
		integer,intent(in) :: inde(:)
		type(DTensornode),pointer ::node,p
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
!	add a DTensornode to the end of the link
!	there is no vector index in the  input DTensor
	subroutine Dpush_backnode(h,node)
		type(DTensorlink),intent(inout) ::h
		type(DTensornode),pointer,intent(inout) ::node
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
!	add a DTensornode to the begin of the link	
!	there is no vector index in the  input DTensor
	subroutine Dpush_fonode(h,node)
		type(DTensorlink),intent(inout) ::h
		type(DTensornode),pointer,intent(inout)  ::node
		type(DTensornode),pointer ::p
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
!	there is no vector index in the  input DTensor
	subroutine Dpush_backnum(h,num)
		type(DTensorlink),intent(inout) ::h
		real*8,intent(in) ::num
		type(DTensor)::Ten
		Ten=(/num/)
		call Dpush_backTen(h,Ten)
		return
	end subroutine
!	add a compelx munber to the begin of the link	
!	there is no vector index in the  input DTensor
	subroutine Dpush_fonum(h,num)
		type(DTensorlink),intent(inout) ::h
		real*8,intent(in) ::num
		type(DTensor)::Ten
		Ten=(/num/)
		call Dpush_fo(h,Ten)
		return
	end subroutine
!	add a compelx munber to the end of the link
!	there is vector index in the  input DTensor
	subroutine Dpush_backnumI(h,num,inde)
		type(DTensorlink),intent(inout) ::h
		real*8,intent(in) ::num
		integer,intent(in) :: inde(:)
		type(DTensor)::Ten
		Ten=(/num/)
		call Dpush_backI(h,Ten,inde)
		return
	end subroutine
!	add a compelx munber to the begin of the link	
!	there is vector index in the  input DTensor
	subroutine Dpush_fonumI(h,num,inde)
		type(DTensorlink),intent(inout) ::h
		real*8,intent(in) ::num
		integer,intent(in) :: inde(:)
		type(DTensor)::Ten
		Ten=(/num/)
		call Dpush_foI(h,Ten,inde)
		return
	end subroutine
!	return the inde DTensor in the link
	type(DTensor) function  DTi(h,inde)		
		type(DTensorlink),intent(in) :: h
		integer,intent(in) :: inde
		integer ::i
		type(DTensornode),pointer ::p
		if(h%length.lt.inde) then
			write(*,*)"ERROR in Ti"
		else
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			DTi=p%Ten
		end if
		return
	end  function
!	return the inde DTensor in the link
!	inde is vector	
	type(DTensor) function  DTiI(h,inde)		
		type(DTensorlink),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i
		type(DTensornode),pointer ::p
		logical::continu
		continu=.true.
		if(DlenOfIndice(h).ne.size(inde)) then
			write(*,*)"ERROR in TiI"
		else
			p=>h%head
			i=1
			do while (continu)
				if(inde.equ.p%indice) then
					DTiI=p%Ten
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
!	return the inde Tensor in the link
!	inde is vector	
! if no such Tensor,return .false.
! else retrun .true.
	logical function  DifTiI(TiI,h,inde)		
		type(DTensorlink),intent(in) :: h
		type(DTensor),intent(out)::TiI
		integer,intent(in) :: inde(:)
		integer ::i
		type(DTensornode),pointer ::p
		logical::continu
		continu=.true.
		DifTiI=.false.
		if(DlenOfIndice(h).ne.size(inde)) then
			write(*,*)"ERROR in TiI"
		else
			p=>h%head
			i=1
			do while (continu)
				if(inde.equ.p%indice) then
					TiI=p%Ten
					DifTiI=.true.
					return
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
!	return the inde DTensornode's address in the link
!	on return,output is a pointer
	 subroutine  Dnode_i_int(h,p,inde)
	 	type(DTensornode),pointer,intent(inout)::p
		type(DTensorlink),intent(in) :: h
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
!	return the inde DTensornode's address in the link
!	inde is vector	
!	on return,output is a pointer
	subroutine  Dnode_i_vec(h,p,inde) 
		type(DTensorlink),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i
		type(DTensornode),pointer,intent(inout) ::p
		logical::continu
		continu=.true.
		if(DlenOfIndice(h).ne.size(inde)) then
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
					write(*,*)"ERROR in Dnode_i_vec"
					stop
				end if
			end do
			
		end if
		return
	end  subroutine	
! return the length of the index of the first Tesnor
	integer function DlenOfIndice(h)
		type(DTensorlink),intent(in) :: h
		type(DTensornode),pointer ::p
		p=>h%head
		DlenOfIndice=size(p%indice)
		return
	end function
!	return the index of the inde DTensor
!	T_{1,1}->T_{1,2}->T_{1,3}->T_{2,1}->T_{2,2}
!	inde=2 => indice=[1,2]
! on entry the size of indice should be equal to the one in h
	subroutine DTenIndice(h,inde,indice)
		type(DTensorlink),intent(in) :: h
		integer,intent(in) :: inde
		integer,intent(out) :: indice(:)
		integer :: i,lenOfind
		type(DTensornode),pointer ::p
		lenOfind=size(indice)
		if(lenOfind.ne.DlenOfIndice(h)) then
			write(*,*) "Error in DTenIndice"
			write(*,*) lenOfind,DlenOfIndice(h)
			p=>h%head
			p=>p%Next
			write(*,*)p%indice
			write(*,*)"stop"
			stop
		end if
		if(h%length.lt.inde) then
			write(*,*)"ERROR in DTenIndice"
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
	logical function DifDnextnode(p)
		type(DTensornode),pointer,intent(inout) ::p
		if(associated(p%next)) then
			DifDnextnode=.true.
			p=>p%next
		else
			DifDnextnode=.false.
		end if
	end function
!	p is a pointer of DTensornode,p points to next
	subroutine Dnextnode(p)
		type(DTensornode),pointer,intent(inout) ::p
		if(associated(p%next)) then
			p=>p%next
		else
			write(*,*)"Error in Dnextnode,p is pointing to the end of the link"
			stop
		end if
	end subroutine
!	p is a pointer of DTensornode,p points to head of the h	
	subroutine Dheadnode(h,p)
		type(DTensorlink),intent(in)::h
		type(DTensornode),pointer,intent(inout) ::p
		p=>h%head
	end subroutine
!	p is a pointer of DTensornode,p points to end of the h	
	subroutine Dendnode(h,p)
		type(DTensorlink),intent(in)::h
		type(DTensornode),pointer,intent(inout) ::p
		p=>h%Tend
	end subroutine
!	add a empty node to the link
	subroutine Daddnode(h)
		type(DTensorlink),intent(inout)::h
		type(DTensornode),pointer ::p
		allocate(p)
		call Dpush_back(h,p)
		return
	end subroutine
		
		
!	get the index of the inde DTensor,and return .true.
!	T_{1,1}->T_{1,2}->T_{1,3}->T_{2,1}->T_{2,2}
!	inde=2 => indice=[1,2]
!	if there no index in h,return .false.
	logical function DTenIndiceLog(h,inde,indice)
		type(DTensorlink),intent(in) :: h
		integer,intent(in) :: inde
		integer,allocatable,intent(out) :: indice(:)
		integer :: i
		type(DTensornode),pointer ::p
		if(h%length.lt.inde) then!if[2]
			write(*,*)"ERROR in DTenIndiceLog"
		else!if[2]
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			if(allocated(p%indice)) then!if[3]
				DTenIndiceLog=.true.
				allocate(indice(size(p%indice)))
				indice=p%indice
			else!if[3]
				DTenIndiceLog=.false.
			end if!if[3]
		end if!if[2]
		return
	end function
!check and output length of the link		
!check if the length of the link is equal to the Dvalue in head of the link
	integer function DChecklength(link)
		type(DTensorlink),intent(in) :: link
		type(DTensornode),pointer ::p
		integer ::num
		DChecklength=link%length
		p=>link%head
		num=0
		do while(associated(p))
			p=>p%next
			num=num+1
		end do
		if(num.ne.DChecklength)then
			write(*,*)"The length of link is",num
			write(*,*)"The Dvalue of length in link is" ,DChecklength
			write(*,*)"Error in length of link"
			write(*,*)"stop !"
			stop
		end if
		return
	end function
!	output length of the link		
	integer function dlinklength(link)
		type(DTensorlink),intent(in) :: link
		dlinklength=link%length
		return
	end function	
	
!	Dmodify the inde DElement of the link
!	The DTensor will be Dmodify
	subroutine Dmodifylink(h,T,inde)
		type(DTensorlink),intent(inout):: h
		integer,intent(in) :: inde
		type(DTensor),intent(in):: T
		integer ::i
		type(DTensornode),pointer ::p
		if(h%length.lt.inde) then
			write(*,*)"ERROR in Dmodify,length"
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
	subroutine Dcleanlink(h)
		type(DTensorlink),intent(inout) :: h
		type(DTensornode),pointer ::p,p1
		if(h%length.eq.0) then!if[2]
			nullify(h%head)
			nullify(h%Tend)
			return
		end if!if[2]
		p=>h%head
		do while(associated(p))
			p1=>p%next
			call cleanDTensor(p%Ten)
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
	subroutine DcopyLink(hout,hin)
		type(DTensorlink),intent(inout) :: hout
		type(DTensorlink),intent(in) :: hin
		integer :: i
		integer,allocatable::indes(:)
		logical :: logi
		call Dcleanlink(hout)
		do i=1,hin%length
			logi=DTenIndiceLog(hin,i,indes)
			if(logi) then!if[2]
				call Dpush_backI(hout,hin.i.i,indes)
			else
				call Dpush_back(hout,hin.i.i)
			end if!if[2]
		end do
		return
	end subroutine
	subroutine DcopyLinkhead(hout,hin)
		type(DTensorlink),intent(inout) :: hout
		type(DTensorlink),intent(in) :: hin
		integer :: i
		integer,allocatable::indes(:)
		logical :: logi
		call Dcleanlink(hout)
		hout%head=>hin%head
		hout%Tend=>hin%Tend
		hout%length=hin%length
		return
	end subroutine
!	connect two links:
!			link1 :T1->T2->T4->T5...->TN
!			link2: TT1->TT2->..->TTM
!			result link will b2 T1->T2->..->TN->TT1->TT2->...->TTM
	type(DTensorlink) function dconnectlink(link1,link2)
		type(DTensorlink),intent(in) :: link1,link2
		type(DTensorlink) :: Templink
		type(DTensornode),pointer ::p
		call DcopyLink(dconnectlink,link1)
		call DcopyLink(Templink,link2)
		p=>dconnectlink%Tend
		p%next=>Templink%head
		dconnectlink%Tend=>Templink%Tend
		dconnectlink%length=dconnectlink%length+Templink%length
		return
	end function
!***************  Ddeletelink  **********************
!	link   :	T1->T2->...->T_{inde1-1}->T_inde1->T_{inde1+1}->....->T_{inde2-1}->T_inde2->T_{inde2+1}...->TN 
!  result:		T1->T2->...->T_{inde1-1}->T_{inde1+1}....->T_{inde2-1}->->T_{inde2+1}...->TN 
!  delete the inde1 to inde2 DTensor in the link,note that the deleted DTensor are include the DTensors of inde1 and inde2
	type(DTensorlink) function Ddeletelink(Dlink_in,inde1,inde2) result(link)			
		type(DTensorlink),intent(in) :: Dlink_in
		integer,intent(in)::inde1,inde2
		type(DTensorlink) :: Templink
		type(DTensornode),pointer ::p1,p2,p
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
		call DcopyLink(link,Dlink_in)
		if(inde1.eq.1) then!if[4]
			if(inde2.eq.link%length) then!if[5]
				call Dcleanlink(link)
				return
			end if!if[5]
			p2=>link%head
			do i=1,inde2
				p=>p2
				p2=>p2%next
				call cleanDTensor(p%Ten)
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
				call cleanDTensor(p%Ten)
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
			call cleanDTensor(p%Ten)
			deallocate(p)
		end do
		p1%next=>p2
		link%length=link%length-(inde2-inde1)-1
		return
	end function		

!DTensor=DTensorlink,if every DElement(DTensor) in the link is a one-DElement DTensor
	subroutine DTenLink(Ten,link)
		type(DTensor),intent(inout)::Ten
		type(DTensorlink),intent(in)::link
		real*8,allocatable::DTensordata(:)
		integer::Tenlen,i
		logical::goon
		type(DTensornode),pointer::p
		Tenlen=dlinklength(link)
		allocate(DTensordata(Tenlen))
		call Dheadnode(link,p)
		do i=1,Tenlen
			if(DgetTotalData(p%Ten).eq.1)then
				DTensordata(i)=p%Ten.i.1
			else
				write(*,*)"error in assignment of a link to DTensor"
			end if
			goon=DifDnextnode(p)
		end do
		Ten=DTensordata
		nullify(p)
		return
	end subroutine
	
	subroutine DLprint(link)
		type(DTensorlink),intent(in) :: link
		integer::lenlink,i
		logical::goon,goon2
		integer,allocatable::indice(:)
		lenlink=DChecklength(link)
		if(lenlink.eq.0) then
			write(*,*)"There is no DTensor in the link"
			return
		end if
		write(*,*) "length of the link is",lenlink
		goon=DTenIndiceLog(link,1,indice)
		if(goon) then
			write(*,*)"indice of every DTensor in the link are"
			do i=1,lenlink
				goon2=DTenIndiceLog(link,i,indice)
				if(goon2) then
					write(*,*)indice
				else
					write(*,*)"no indece"
				end if
			end do
		else
			goon=.true.
			do i=1,lenlink
				goon2=DTenIndiceLog(link,i,indice)
				if(goon2) then
					write(*,*)"error in the link,there are some DTensor with indice,some are not"
					goon=.false.
				end if
			end do
		end if
		return
	end subroutine
	
	subroutine DLDprint(link)
		type(DTensorlink),intent(in) :: link
		type(DTensornode),pointer::p
		integer::lenlink,i
		logical::goon
		integer,allocatable::indice(:)
		lenlink=DChecklength(link)
		if(lenlink.eq.0) then
			write(*,*)"There is no DTensor in the link"
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
	 
! 	DTensorlist and Dlistnode 	  





!	add a DTensorlink to the end of the list
!	there is no vector index in the  input DTensorlink
! if pointtoL=true,than the new node will point to L
!  otherwhile ,create a new link
	subroutine Dpush_backlink(h,L,pointtoL)
		type(DTensorlist),intent(inout) ::h
		type(DTensorlink),intent(in):: L
		logical,intent(in)::pointtoL
		type(Dlistnode),pointer ::node
		allocate(node)
		if(pointtoL) then
		   call DcopyLinkhead(node%link,L)
		else
   		call DcopyLink(node%link,L)
   	end if
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
!	add a DTensorlink to the begin of the list	
!	there is no vector index in the  input DTensorlink
	subroutine Dpush_folink(h,L,pointtoL)
		type(DTensorlist),intent(inout) ::h
		type(DTensorlink),intent(in):: L
		logical,intent(in)::pointtoL
		type(Dlistnode),pointer ::node,p
		allocate(node)
		if(pointtoL)then
   		call DcopyLinkhead(node%link,L)
   	else
   	   call DcopyLink(node%link,L)
   	end if
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
!	add a DTensorlink to the end of the list	
!	inde is the index of the DTensorlink
	subroutine Dpush_backlinkI(h,L,inde,pointtoL)
		type(DTensorlist),intent(inout) ::h
		type(DTensorlink),intent(in):: L
		integer,intent(in) :: inde(:)
		logical,intent(in)::pointtoL
		type(Dlistnode),pointer ::node
		allocate(node)
		if(pointtoL)then
   		call DcopyLinkhead(node%link,L)
   	else
   	   call DcopyLink(node%link,L)
   	end if
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
!	add a DTensorlink to the begin of the list	
!	inde is the index of the DTensorlink	
	subroutine Dpush_folinkI(h,L,inde,pointtoL)
		type(DTensorlist),intent(inout) ::h
		type(DTensorlink),intent(in):: L
		integer,intent(in) :: inde(:)
		logical,intent(in)::pointtoL
		type(Dlistnode),pointer ::node,p
		allocate(node)
		if(pointtoL)then
   		call DcopyLinkhead(node%link,L)
   	else
   	   call DcopyLink(node%link,L)
   	end if
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
!	add a Dlistnode to the end of the list
!	there is no vector index in the  input DTensorlink
	subroutine Dpush_backlinknode(h,node)
		type(DTensorlist),intent(inout) ::h
		type(Dlistnode),pointer,intent(inout) ::node
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
!	add a Dlistnode to the end of the list
!	there is no vector index in the  input DTensor
	subroutine Dpush_folinknode(h,node)
		type(DTensorlist),intent(inout) ::h
		type(Dlistnode),pointer,intent(inout)  ::node
		type(Dlistnode),pointer ::p
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
!	return the inde DTensorlink in the list
	function  DLi(h,inde)		result(link)
		type(DTensorlink):: link	
		type(DTensorlist),intent(in) :: h
		integer,intent(in) :: inde
		integer ::i
		type(Dlistnode),pointer ::p
		if(h%length.lt.inde) then
			write(*,*)"ERROR in Ti"
		else
			call Dcleanlink(link)
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			link=p%link
		end if 
		return
	end  function
!	return the inde DTensorlink in the list
!	inde is vector	
	function  DLiI(h,inde)		result(link)
		type(DTensorlink):: link
		type(DTensorlist),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i
		type(Dlistnode),pointer ::p
		logical::continu
		continu=.true.
		call Dcleanlink(link)
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
	end  function
!	return the inde DTensorlink in the list
	subroutine  DcopyLi(link,h,inde)	
		type(DTensorlink),intent(inout) :: link	
		type(DTensorlist),intent(in) :: h
		integer,intent(in) :: inde
		integer ::i
		type(Dlistnode),pointer ::p
		if(h%length.lt.inde) then
			write(*,*)"ERROR in Ti"
		else
			call Dcleanlink(link)
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			call DcopyLink(link,p%link)
		end if 
		return
	end  subroutine
!	return the inde DTensorlink in the list
!	inde is vector	
	subroutine  DcopyLiI(link,h,inde)		
		type(DTensorlink),intent(inout) :: link
		type(DTensorlist),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i
		type(Dlistnode),pointer ::p
		logical::continu
		continu=.true.
		call Dcleanlink(link)
		p=>h%head
		i=1
		do while (continu)
			if(inde.equ.p%indice) then
				call DcopyLink(link,p%link)
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
!	return the inde DTensorlink in the list
	type(DTensor) function  DLiv(h,inde)	
		type(DTensorlist),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i,j
		type(Dlistnode),pointer ::p
		type(DTensornode),pointer ::Tp
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
			DLiv=Tp%Ten
		end if 
		return
	end  function
!	return the inde DTensorlink in the list
!	inde is vector	
	type(DTensor) function  DLiIv(h,inde)		
		type(DTensorlist),intent(in) :: h
		integer,intent(in) :: inde(:,:)
		integer ::i,j,size1,size2
		type(Dlistnode),pointer ::p
		type(DTensornode),pointer ::Tp
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
!search DTensor		
		if(inde(2,1).eq.0)then!if[2]
				do j=1,inde(2,2)-1
					Tp=>Tp%next
				end do
				DLiIv=Tp%Ten
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
						DLiIv=Tp%Ten
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
!	return the inde Dlistnode's address in the list
!	on return,output is a pointer
	 subroutine  DlistDnode_i_int(h,p,inde)
	 	type(Dlistnode),pointer,intent(inout)::p
		type(DTensorlist),intent(in) :: h
		integer,intent(in) :: inde
		integer ::i
		if(h%length.lt.inde) then
			write(*,*)"ERROR in DlistDnode_i_int"
		else
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
		end if
		return
	end  subroutine
!	return the inde Dlistnode's address in the link
!	inde is vector	
!	on return,output is a pointer
	subroutine  DlistDnode_i_vec(h,p,inde) 
		type(DTensorlist),intent(in) :: h
		integer,intent(in) :: inde(:)
		integer ::i
		type(Dlistnode),pointer,intent(inout) ::p
		logical::continu
		continu=.true.
		if(DlenOfIndicelist(h).ne.size(inde)) then
			write(*,*)"ERROR in DlistDnode_i_vec"
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
					write(*,*)"ERROR in DlistDnode_i_vec"
					stop
				end if
			end do
			
		end if
		return
	end  subroutine	
! return the length of the index of the first Tesnorlink
	integer function DlenOfIndicelist(h)
		type(DTensorlist),intent(in) :: h
		type(Dlistnode),pointer ::p
		p=>h%head
		DlenOfIndicelist=size(p%indice)
		return
	end function
!	return the index of the inde DTensorlink
!	T_{1,1}->T_{1,2}->T_{1,3}->T_{2,1}->T_{2,2}
!	inde=2 => indice=[1,2]
! on entry the size of indice should be equal to the one in h
	subroutine DlinkIndice(h,inde,indice)
		type(DTensorlist),intent(in) :: h
		integer,intent(in) :: inde
		integer,intent(out) :: indice(:)
		integer :: i,lenOfind
		type(Dlistnode),pointer ::p
		lenOfind=size(indice)
		if(lenOfind.ne.DlenOfIndicelist(h)) then
			write(*,*) "Error in DlinkIndice"
			write(*,*) lenOfind,DlenOfIndicelist(h)
			p=>h%head
			p=>p%Next
			write(*,*)p%indice
			write(*,*)"stop"
			stop
		end if
		if(h%length.lt.inde) then
			write(*,*)"ERROR in DlinkIndice"
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
	logical function Difnextlistnode(p)
		type(Dlistnode),pointer,intent(inout) ::p
		if(associated(p%next)) then
			Difnextlistnode=.true.
			p=>p%next
		else
			Difnextlistnode=.false.
		end if
	end function
!	p is a pointer of Dlistnode,p points to next
	subroutine DnextDlistnode(p)
		type(Dlistnode),pointer,intent(inout) ::p
		if(associated(p%next)) then
			p=>p%next
		else
			write(*,*)"Error in next;istnode,p is pointing to the end of the list"
			stop
		end if
	end subroutine
!	p is a pointer of Dlistnode,p points to head of the h	
	subroutine Dheadlist(h,p)
		type(DTensorlist),intent(in)::h
		type(Dlistnode),pointer,intent(inout) ::p
		p=>h%head
		return
	end subroutine
!	p is a pointer of Dlistnode,p points to end of the h	
	subroutine Dendlist(h,p)
		type(DTensorlist),intent(in)::h
		type(Dlistnode),pointer,intent(inout) ::p
		p=>h%Lend
	end subroutine
!	add a empty Dlistnode to the list
	subroutine Daddlist(h)
		type(DTensorlist),intent(inout)::h
		type(Dlistnode),pointer ::p
		allocate(p)
		call Dpush_backlinknode(h,p)
		return
	end subroutine
		
		
!	get the index of the inde DTensorlink,and return .true.
!	T_{1,1}->T_{1,2}->T_{1,3}->T_{2,1}->T_{2,2}
!	inde=2 => indice=[1,2]
!	if there no index in h,return .false.
	logical function DlinkIndiceLog(h,inde,indice)
		type(DTensorlist),intent(in) :: h
		integer,intent(in) :: inde
		integer,allocatable,intent(out) :: indice(:)
		integer :: i
		type(Dlistnode),pointer ::p
		if(h%length.lt.inde) then!if[2]
			write(*,*)"ERROR in DTenIndiceLog"
		else!if[2]
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			if(allocated(p%indice)) then!if[3]
				DlinkIndiceLog=.true.
				allocate(indice(size(p%indice)))
				indice=p%indice
			else!if[3]
				DlinkIndiceLog=.false.
			end if!if[3]
		end if!if[2]
		return
	end function
!check and output length of the link		
!check if the length of the link is equal to the Dvalue in head of the link
	integer function DChecklistlength(list)
		type(DTensorlist),intent(in) :: list
		type(Dlistnode),pointer ::p
		integer ::num
		DChecklistlength=list%length
		p=>list%head
		num=0
		do while(associated(p))
			p=>p%next
			num=num+1
		end do
		if(num.ne.DChecklistlength)then
			write(*,*)"The length of link is",num
			write(*,*)"The Dvalue of length in link is" ,DChecklistlength
			write(*,*)"Error in length of link"
			write(*,*)"stop !"
			stop
		end if
		return
	end function
!	output length of the link		
	integer function Dlistlength(list)
		type(DTensorlist),intent(in) :: list
		Dlistlength=list%length
		return
	end function	
	
!	Dmodify the inde DElement of the list
!	The DTensorlink will be Dmodify
	subroutine Dmodifylist(h,L,inde)
		type(DTensorlist),intent(inout):: h
		integer,intent(in) :: inde
		type(DTensorlink),intent(in):: L
		integer ::i
		type(Dlistnode),pointer ::p
		if(h%length.lt.inde) then
			write(*,*)"ERROR in Dmodifylink,length"
		else
			p=>h%head
			do i=1,inde-1
				p=>p%next
			end do
			call DcopyLink(p%link,L)
		end if
		return
	end subroutine
!	clean the list
	subroutine Dcleanlist(h)
		type(DTensorlist),intent(inout) :: h
		type(Dlistnode),pointer ::p,p1
		if(h%length.eq.0) then!if[2]
			nullify(h%head)
			nullify(h%Lend)
			return
		end if!if[2]
		p=>h%head
		do while(associated(p))
			p1=>p%next
			call Dcleanlink(p%link)
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
	subroutine DcopyList(hout,hin)
		type(DTensorlist),intent(inout) :: hout
		type(DTensorlist),intent(in) :: hin
		type(DTensorlink)::tempL
		integer :: i
		integer,allocatable::indes(:)
		logical :: logi
		call Dcleanlist(hout)
		do i=1,hin%length
			logi=DlinkIndiceLog(hin,i,indes)
			call Dlink_i(tempL,hin,i)
			if(logi) then!if[2]
				call Dpush_backlinkI(hout,tempL,indes,.false.)
			else
				
				call Dpush_backlink(hout,tempL,.false.)
			end if!if[2]
		end do
		call Dcleanlink(tempL)
		return
	end subroutine
	subroutine DcopyListhead(hout,hin)
	   type(DTensorlist),intent(inout) :: hout
		type(DTensorlist),intent(in) :: hin
		hout%head=>hin%head
		hout%Lend=>hin%Lend
		hout%length=hin%length
		return
	end subroutine
!	connect two lists:
!			list1 :T1->T2->T4->T5...->TN
!			list2: TT1->TT2->..->TTM
!			result list will b2 T1->T2->..->TN->TT1->TT2->...->TTM
	type(DTensorlist) function Dconnectlist(list1,list2)
		type(DTensorlist),intent(in) :: list1,list2
		type(DTensorlist) :: Templist
		type(Dlistnode),pointer ::p
		call DcopyList(Dconnectlist,list1)
		call DcopyList(Templist,list2)
		p=>Dconnectlist%Lend
		p%next=>Templist%head
		Dconnectlist%Lend=>Templist%Lend
		Dconnectlist%length=Dconnectlist%length+Templist%length
		return
	end function
!***************  Ddeletelink  **********************
!	list   :	T1->T2->...->T_{inde1-1}->T_inde1->T_{inde1+1}->....->T_{inde2-1}->T_inde2->T_{inde2+1}...->TN 
!  result:		T1->T2->...->T_{inde1-1}->T_{inde1+1}....->T_{inde2-1}->->T_{inde2+1}...->TN 
!  delete the inde1 to inde2 DTensor in the link,note that the deleted DTensor are include the DTensors of inde1 and inde2
	type(DTensorlist) function Ddeletelist(list_in,inde1,inde2) result(list)			
		type(DTensorlist),intent(in) :: list_in
		integer,intent(in)::inde1,inde2
		type(DTensorlist) :: Templink
		type(Dlistnode),pointer ::p1,p2,p
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
		call DcopyList(list,list_in)
		if(inde1.eq.1) then!if[4]
			if(inde2.eq.list%length) then!if[5]
				call Dcleanlist(list)
				return
			end if!if[5]
			p2=>list%head
			do i=1,inde2
				p=>p2
				p2=>p2%next
				call Dcleanlink(p%link)
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
				call Dcleanlink(p%link)
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
			call Dcleanlink(p%link)
			deallocate(p)
		end do
		p1%next=>p2
		list%length=list%length-(inde2-inde1)-1
		return
	end function		
	
	subroutine DListprint(list)
		type(DTensorlist),intent(in) :: list
		integer::lenlist,i
		logical::goon,goon2
		integer,allocatable::indice(:)
		lenlist=DChecklistlength(list)
		if(lenlist.eq.0) then
			write(*,*)"There is no DTensorlink in the list"
			return
		end if
		write(*,*) "length of the list is",lenlist
		goon=DlinkIndiceLog(list,1,indice)
		if(goon) then
			write(*,*)"indice of every DTensorlink in the list are"
			do i=1,lenlist
				goon2=DlinkIndiceLog(list,i,indice)
				if(goon2) then
					write(*,*)indice
				else
					write(*,*)"no indece"
				end if
			end do
		else
			goon=.true.
			do i=1,lenlist
				goon2=DlinkIndiceLog(list,i,indice)
				if(goon2) then
					write(*,*)"error in the list,there are some DTensorlink with indice,some are not"
					goon=.false.
				end if
			end do
		end if
		return
	end subroutine
	
	subroutine DListlinkprint(list)
		type(DTensorlist),intent(in) :: list
		type(Dlistnode),pointer::p
		integer::lenlist,i
		logical::goon
		integer,allocatable::indice(:)
		lenlist=DChecklistlength(list)
		if(lenlist.eq.0) then
			write(*,*)"There is no DTensor in the link"
			return
		end if
		p=>list%head
		call DLDprint(p%link)
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
			call DLDprint(p%link)
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
	subroutine Dwarning(T,chara)
		character*50,intent(in) :: chara
		type(DTensor),intent(in) :: T
		real*8	::maxV
		maxV=dmaxabsElement(T)
		if(maxV.ge.large_number_warning) then
			write (*,*) "***   warning   ***"
			write (*,*) chara
			write (*,*) "data in the DTensor is large,the abs Dvalue largest DElement is"
			write(*,*)		maxV
		end if
	end subroutine
	subroutine Dwarning2(T)
		type(DTensor),intent(in) :: T
		real*8	::maxV
		maxV=dmaxabsElement(T)
		if(maxV.ge.large_number_warning) then
			write (*,*) "***   warning   ***"
			write (*,*) "data in the DTensor is large,the abs Dvalue largest DElement is"
			write(*,*)		maxV
		end if
	end subroutine	
!*****************  NANwarning   *****************
!			warning		NAN
	subroutine DNANwarning(T,chara)
		character*50,intent(in) :: chara
		type(DTensor),intent(in) :: T
		real*8::val
		integer ::i,total
		total=T%TotalData
		do i=1,total
			val=T.i. i
			if(isnan(val)) then
				write (*,*) "***   NANwarning   ***"
				write (*,*) chara
				return
			end if
		end do
		return
	end subroutine
	subroutine DNANwarning2(T)
		type(DTensor),intent(in) :: T
		real*8::val
		integer ::i,total
		total=T%TotalData
		do i=1,total
			val=T.i. i
			if(isnan(val)) then
				write (*,*) "***   NANwarning   ***"
				write (*,*) "stop"
				stop
			end if
		end do
		return
	end subroutine	
!*****************  NANjudge   *****************
!			warning	NAN
	logical function DNANjudge(T)
		type(DTensor),intent(in) :: T
		integer ::i,total
		real*8::val
		DNANjudge=.false.
		total=T%TotalData
		do i=1,total
			val=T.i. i
			if(isnan(val))  then
				DNANjudge=.true.
			return
			end if
		end do
		return
	end function
	logical function DNANjudgev(val)
		real*8,intent(in) :: val(:)
		integer ::i,total
		DNANjudgev=.false.
		total=size(val)
		do i=1,total
			if(isnan(val(i))) then
				DNANjudgev=.true.
			return
			end if
		end do
		return
	end function			

   subroutine Dcheck_zero_Ten(T)!in case of the DElement of rho is -1*10^-70,this is zero
	   type(DTensor),intent(inout)::T
	   integer::i
	   do i=1,DgetTotalData(T)
         if(abs(T%DTensor_Data(i)).lt.zero_num) then
            T%DTensor_Data(i)=0d0
	      end if
	   end do
	end subroutine



!**********************************************************************
!**********************************************************************
!	the code below is for MPI
!**********************************************************************
	subroutine sent_DTensor(Ten1,Ten2,ID1,ID2,ierr)
		type(DTensor),intent(in)::Ten1
		type(DTensor),intent(inout)::Ten2
		integer,intent(in)::ID1,ID2,ierr
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE)
		tag=1
		call mpi_comm_rank(mpi_comm_world,proID,ierr)
		call mpi_comm_size(mpi_comm_world,proNum,ierr )
		if(proID.eq.ID1) then
			call mpi_send(Ten1%flag,1,MPI_logical,ID2,tag,MPI_Comm_world,ierr)
			if(.not.Ten1%flag)then
				return
			end if
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Ten2%flag,1,MPI_logical,ID1,tag,MPI_Comm_world,istatus,ierr)
			if(.not.Ten2%flag)then
				return
			end if
		end if
!**********************************************************************************			
		if(proID.eq.ID1) then
			call mpi_send(Ten1%rank,1,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Ten2%rank,1,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
		end if
!**********************************************************************************		
		call sent_Dimension(Ten1%TenDim,Ten2%TenDim,ID1,ID2,ierr)
!**********************************************************************************		
		if(proID.eq.ID1) then
			call mpi_send(Ten1%totalData,1,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Ten2%totalData,1,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
		end if
!**********************************************************************************	
		if(proID.eq.ID1) then
			call mpi_send(Ten1%flag,1,MPI_logical,ID2,tag,MPI_Comm_world,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Ten2%flag,1,MPI_logical,ID1,tag,MPI_Comm_world,istatus,ierr)
		end if
!**********************************************************************************	
		if(proID.eq.ID1) then
			call mpi_send(Ten1%DTensor_data,Ten1%totalData,MPI_double_precision,ID2,tag,MPI_Comm_world,ierr)
		end if
		if(proID.eq.ID2) then
			if(allocated(Ten2%DTensor_data)) then
				deallocate(Ten2%DTensor_data)
			end if
			allocate(Ten2%DTensor_data(Ten2%totalData))
			call mpi_recv(Ten2%DTensor_data,Ten2%totalData,MPI_double_precision,ID1,tag,MPI_Comm_world,istatus,ierr)
		end if
!**********************************************************************************	
		return
	end subroutine	
	
	subroutine BCAST_DTensor(Ten1,ID,ierr)
		type(DTensor),intent(inout)::Ten1
		integer,intent(in)::ID,ierr
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE)
		tag=1
		call mpi_comm_rank(mpi_comm_world,proID,ierr)
		call mpi_comm_size(mpi_comm_world,proNum,ierr )
		
		call MPI_BCAST(Ten1%flag,1,MPI_logical,ID,MPI_COMM_WORLD,ierr)
		call MPI_Barrier(MPI_COMM_WORLD,ierr)
		if(.not.Ten1%flag) then
			return
		end if
!**********************************************************************************		
		call MPI_BCAST(Ten1%rank,1,MPI_integer,ID,MPI_COMM_WORLD,ierr)	
!**********************************************************************************		
		call BCAST_Dimension(Ten1%TenDim,ID,ierr)
!**********************************************************************************		
		call MPI_BCAST(Ten1%totalData,1,MPI_integer,ID,MPI_COMM_WORLD,ierr)
!**********************************************************************************	
		if(proID.ne.ID) then
			if(allocated(Ten1%DTensor_data)) then
				deallocate(Ten1%DTensor_data)
			end if
			allocate(Ten1%DTensor_data(Ten1%totalData))
		end if
		call MPI_BCAST(Ten1%DTensor_data,Ten1%totalData,MPI_double_precision,ID,MPI_COMM_WORLD,ierr)
!**********************************************************************************	
		return
	end subroutine
	
	subroutine sent_DTensorlink(link1,link2,ID1,ID2,ierr)
		type(DTensorlink),intent(in)::link1
		type(DTensorlink),intent(inout)::link2
		integer,intent(in)::ID1,ID2,ierr
		integer::proID,proNum,tag,linklen,indeLen,i,istatus(MPI_STATUS_SIZE)
		type(DTensor)::Ten
		integer,allocatable:: indice(:)
		logical::is_index
		tag=1
		linklen=0
		call mpi_comm_rank(mpi_comm_world,proID,ierr)
		call mpi_comm_size(mpi_comm_world,proNum,ierr )
		if(proID.eq.ID1) then
			linklen=Dlinklength(link1)
			call mpi_send(linklen,1,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
			if(linklen.eq.0) then
				return
			end if
		end if
		if(proID.eq.ID2) then
			call mpi_recv(linklen,1,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
			if(linklen.eq.0) then
				return
			end if
			call Dcleanlink(link2)
		end if
!**********************************************************************************
		do i=1,linklen
			if(proID.eq.ID1) then
				is_index=DTenIndiceLog(link1,i,indice)
				call mpi_send(is_index,1,MPI_logical,ID2,tag,MPI_Comm_world,ierr)
			end if	
			if(proID.eq.ID2) then
				call mpi_recv(is_index,1,MPI_logical,ID1,tag,MPI_Comm_world,istatus,ierr)
			end if
			if(is_index) then
				if(proID.eq.ID1) then
					indeLen=size(indice)
					call mpi_send(indeLen,1,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
				end if
				if(proID.eq.ID2) then
					call mpi_recv(indeLen,1,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
				end if
				if(proID.eq.ID1) then
					call mpi_send(indice,indeLen,MPI_integer,ID2,tag,MPI_Comm_world,ierr)
				end if
				if(proID.eq.ID2) then
					if(allocated(indice)) then
						deallocate(indice)
					end if
					allocate(indice(indeLen))
					call mpi_recv(indice,indeLen,MPI_integer,ID1,tag,MPI_Comm_world,istatus,ierr)
				end if
				if(proID.eq.ID1) then
					Ten=link1.i.i
				end if
				call sent_DTensor(Ten,Ten,ID1,ID2,ierr)
				if(proID.eq.ID2) then
					call Dpush_backI(link2,Ten,indice)
				end if
			else
				if(proID.eq.ID1) then
					Ten=link1.i.i
				end if
				call sent_DTensor(Ten,Ten,ID1,ID2,ierr)
				if(proID.eq.ID2) then
					call Dpush_back(link2,Ten)
				end if
			end if
		end do
		return
	end subroutine
	
	
	subroutine BCAST_DTensorlink(link1,ID,ierr)
		type(DTensorlink),intent(inout)::link1
		integer,intent(in)::ID,ierr
		integer::proID,proNum,tag,linklen,indeLen,i,istatus(MPI_STATUS_SIZE)
		type(DTensor)::Ten
		integer,allocatable:: indice(:)
		logical::is_index
		tag=1
		call mpi_comm_rank(mpi_comm_world,proID,ierr)
		call mpi_comm_size(mpi_comm_world,proNum,ierr )
		if(proID.eq.ID) then
			linklen=Dlinklength(link1)
		else
			call Dcleanlink(link1)
		end if
		call MPI_BCAST(linklen,1,MPI_integer,ID,MPI_COMM_WORLD,ierr)
		if(linklen.eq.0) then
			call Dcleanlink(link1)
			return
		end if
!**********************************************************************************
		do i=1,linklen
			if(proID.eq.ID) then
				is_index=DTenIndiceLog(link1,i,indice)
			end if
			call MPI_BCAST(is_index,1,MPI_logical,ID,MPI_COMM_WORLD,ierr)
			
			if(is_index) then
			
				if(proID.eq.ID) then
					indeLen=size(indice)
				end if
				call MPI_BCAST(indeLen,1,MPI_integer,ID,MPI_COMM_WORLD,ierr)
				if(proID.ne.ID) then
					if(allocated(indice)) then
						deallocate(indice)
					end if
					allocate(indice(indeLen))
				end if
				call MPI_BCAST(indice,indeLen,MPI_integer,ID,MPI_COMM_WORLD,ierr)
				if(proID.eq.ID) then
					Ten=link1.i.i
				end if
				call BCAST_DTensor(Ten,ID,ierr)
				if(proID.ne.ID) then
					call Dpush_backI(link1,Ten,indice)
				end if
				
			else
				if(proID.eq.ID) then
					Ten=link1.i.i
				end if
				call BCAST_DTensor(Ten,ID,ierr)
				if(proID.ne.ID) then
					call Dpush_back(link1,Ten)
				end if
			end if
			
		end do
		return
	end subroutine
	
	
	
end module
!***************************************************
!***************** END OF TNESOR *************
!***************************************************


			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
