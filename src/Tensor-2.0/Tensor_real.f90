
!		 This file is modify from Tensor.f90 which is used for real*8 type
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
!			5.These files call arpack lapack and blas . only eigen use arpack .
!		gfortran name.f90 -o name -larpack -llapack -lrefblas 
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

module Tensor_real
	use Dimension_typede
	implicit none
!****************************************************	
!*********		define data type	***************
	type  DTensor
		integer,private :: rank!length of Dimension.TenDim%Dimsize
		type(Dimension),private	:: TenDim!dimenison
		real*8,allocatable :: DTensor_data(:)!The data store in one dimenson array
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
		type(DTensorlink),private ::link
		type(Dlistnode),private,pointer :: next
	end type Dlistnode
	
!*********		End define		*********************


!****************************************************
!*********		define private data	***************	
	real*8,private :: cputime(2)
	real*8,parameter,private :: II=(0,1),EE=(1,0)! II**2=-1
	real*8,parameter,private :: large_number_warning=1d80!use for checking
	real*8,parameter,private :: zero_num=1d-15!
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
	private::Dpermutation_data3,Dpermutation_data2,Dpermutefo_data,Dpermuteback_data! Use for permuation and permutefo_vec,permuteback_vec
!*********		End define	***************		
		
	interface DTMprint
		module procedure DTMprint1
		module procedure DTMprint2
		module procedure DTMprint2_file
	end interface
	interface DTprint
		module procedure DTprint1
		module procedure DTprint2
		module procedure DTprint2_file
	end interface

!**************************************************************************
!**************************************************************************
!! make LQ decompostion: A = LQ, Q is an orthomomal matrix, and L is a lower triangle matrix
!LQ is a function that return a DTensorlink
!the interface below is a subroutine, is the same function as DQRlink or DLQdecompose
!call DLQdecomposition(T,Q,R)  or DLQdecomposition(T,Q),L store in T
	interface DLQdecomposition
		module procedure DLQdecomposition1
		module procedure DLQdecomposition2
	end interface
!********************************************************************
!QR decompostion: A = QR, Q is an orthomomal matrix, and R is an upper triangle matrix	
!QR is a function that return a DTensorlink 
!the interface below is a subroutine , is the same function as QRlink or QRdecompose
!call DQRdecomposition(T,Q,R)  or  DQRdecomposition(T,R),Q store in T
	interface DQRdecomposition
		module procedure DQRdecomposition1
		module procedure DQRdecomposition2
	end interface
!Other function DQRdecompose, DLQdecompose:output a array of Tensor
!Other function DQRlink, DLQlink:output a Tensorlink
!**************************************************************************


!*********** SVD **************:

!			link=DSVD(T) or link=DSVD(T,Ncut)
!				T should be two dimension
!				T=U*s*V^T,on output,the order in the link is:[U]->[s]->[V^T]
! 				T=(link.i.1)*(eye(link.i.2))*(link.i.3)
!
	interface DSVDlink
		module procedure Dsvd_no_cut
		module procedure dsvd_cut
	end interface
	
!			Dsvddecomposition(T,U,s,V),or svddecomposition(T,U,s,V,Ncut):
!			T should be two dimension	
! 				T=U*s*V^T,on output T=U*eye(s)*V	
!
	interface DSVDdecomposition
		module procedure DSVDdecomposition1
		module procedure DSVDcutoff
	end interface
!Other function svddecompose	output a array
!**************************************************************************

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
!one can use subroutine setTensorName(T,Tensorname,nameID) to set a name to the Tensor
		module procedure DbuildTen1
		module procedure DbuildTen2
		module procedure DbuildTen3
		module procedure DbuildTen4
		module procedure DbuildTen5
		module procedure DbuildTen6
		module procedure DbuildTen8
		module procedure DbuildTen8_2
	end interface

	interface Dgenerate
		module procedure Dgenerate_NoName!generate a random Tensor (0~1)
		module procedure Dgenerate_Namechar
	end interface
	interface setDTensorName
		module procedure setDTensorName2!set all the name in Tensor call setTensorName(T,'A') or setTensorName(T,'A12')
		module procedure setDTensorName3!modify the name of ith index
		module procedure setDTensorName4!modify the index whose name is old_w to new_w
		
	end interface
	interface Dmodify
		module procedure Dmodifylink
		module procedure Dmodifylist
		module procedure DmodifyTen_val
		module procedure DmodifyTen_val2
		module procedure DmodifyTen_dat
		module procedure DmodifyTen_dat2
		module procedure DmodifyTen_dat3
		module procedure DmodifyTen_some_data1_real
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
		module procedure Deye_real!input s(:) is real*8,a m*n matrix, m,n are input parameter
		module procedure Deye_real2!input s(:) is real*8,output a n*n matrix,n=size(s)
		module procedure Done!identity matrix,return a m*n identity matrix, m,n are input parameter
		module procedure Done_Ten!return a identity DTensor ,intput a vec as the dimension
!										!example,T=one_Ten((/2,3,4/)), then T(i,i,i)=1d0,other T=0d0
		module procedure Deye_Ten1!input DTensor,output diag[T.i.1 , T.i.2 , ... , , T.i.n],rank of T should be 1
		module procedure Deye_Ten2!input DTensor,output a matrix of m*n: diag[T.i.1 , T.i.2 , ... , , T.i.n],rank of T should be 1
	end interface
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
	interface DTenproduct
		module procedure DTenproduct_noName
		module procedure DTenproduct_Name!In put dimension as character
		module procedure DTenproduct_name_int
		module procedure DTenproduct_int_name
		module procedure DTenproduct_old
		module procedure DTenproduct_name_rename1
		module procedure DTenproduct_name_rename2
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
		module procedure Dpermute_rank2 !Dpermute the DTensor whose rank is 2,if rank=1 ,do nothing
		module procedure Dpermute_rank3 !Dpermute the DTensor whose rank is 3
		module procedure Dpermutation	 !Dpermute the DTensor of any rank,give the new order of the dimension
												 !		If operate on a big DTensor,use Dpermute_rank2 or Dpermute_rank3,
												 !	 they are faster.if rank>3,use Dcontract to reshape
		module procedure Dpermutation_name
	end interface
	interface operator(.pf.)
		module procedure Dpermutefo!Dpermute the inde index to the first
										  !T_{1,2,3,..,i,..,n},Dpermutefo(T,i)=_{i,1,2,3,..,i-1,i+1,..,n}
										  ! note the output will Ddecompose all the dimension
		module procedure Dpermutefo_vec
		module procedure Dpermutefo_name
		module procedure Dpermutefo_vec_name
	end interface
	interface operator(.pb.)!	T_{1,2,3,..,i,..,n},permuteback(T,i)=_{1,2,3,..,i-1,i+1,..,n,i}
									! note the output will decompose all the dimension
		module procedure Dpermuteback
		module procedure Dpermuteback_vec
		module procedure Dpermuteback_name
		module procedure Dpermuteback_vec_name
	end interface
	interface operator(.pi.)!T_{1,2,3,..,i,..,n},DpermuteInde(T,i)=_{2,3,..,i,1,i+1,..,n}
									! note the output will Ddecompose all the dimension
		module procedure DpermuteInde
	end interface
	interface operator(.pbi.)!T_{1,2,3,..,i,..,n},permuteInde(T,i)=_{1,2,3,..,i-1,n,i,i+1,..,n-1}
									! note the output will decompose all the dimension
		module procedure DpermutebackInde
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
	
!******************************************************************************************************
!******************************************************************************************************

!								   direct Product
!						There are two type of definition of direct Product
!

!******************************************************************************************************
!The first one,which use in matlab(kron),but may be wrong doing MPS(I do not test)

	interface operator(.kron.)
		module procedure DdirectProduct_Matlab_korn! direct Product [m1,n1] * [m2,n2] => [m2,m1,n2,n1]
												!  or direct Product return a matrix [n1] * [n2] => [ n2,n1 ]
	end interface
		
	interface operator(.mkron.)
		module procedure DdirectProductM_Matlab_korn! direct Product return a matrix [m1,n1] * [m2,n2] => [(m2*m1),(n2*n1)]
												!  or direct Product return a matrix [n1] * [n2] => [ (n2*n1) ]
	end interface
!******************************************************************************************************
!The second one, use for MPS(I have test)	

	interface operator(.xx.)
		module procedure DdirectProduct! direct Product [m1,n1] * [m2,n2] => [m1,m2,n1,n2]
												!  or direct Product return a matrix [n1] * [n2] => [ n1,n2 ]
	end interface
		
	interface operator(.mxx.)
		module procedure DdirectProductM! direct Product return a matrix [m1,n1] * [m2,n2] => [(m1*m2),(n1*n2)]
												!  or direct Product return a matrix [n1] * [n2] => [ (n1*n2) ]
	end interface
	
!******************************************************************************************************
!******************************************************************************************************
		
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
	
	interface operator (.inv.)
		module procedure dinverse! The inverse of a matrix: the input tensor should be a square matrix 
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
	
	interface operator(.iname.)!the index name
		module procedure DoutIndexName
		module procedure DoutAllIndexName
	end interface
	interface operator(.Tname.)!the Tensor name
		module procedure DoutTensorName
		module procedure DoutAllTensorName
	end interface
	
	interface operator(.dim.)
		module procedure DgetTenDim_Namei!output the ith dimension of the Tensor, input the name of the dimension
													!If can not find , output 0
		module procedure DgetTenDim_i!output the i dimension of the DTensor,output an integer
		module procedure DgetTenDim!output the Dimension of Tensor,return type(Dimension)
											!It can use to commit the assignment to a vector
											!vector=type(Dimension),if the vector is allocatable
	end interface
!allocate Tensor according to the dimen		
	interface allocateDTensor
		module procedure  allocateDTensor1
		module procedure  allocateDTensor2
	end interface		
	interface copyDTensor
		module procedure   copyDTensortwoDTensor
		module procedure copyDTensordim1
		module procedure copyDTensordim2
		module procedure copyDTensordim3
	end interface
	interface assignment(=)
		module procedure DassignmentTen!T1=T2 ,both T1 and T2 are DTensor
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
		module procedure DassignmentTenArray!type(Tensor)::T1(len1),T2(len2),
!																			T1=T2,if len1<len2 stop
!																			if len1>len2 T1(1:len2)=T2,and on data in T1(len2+1:)
!																			That means Dgetflag(T1(len2+1))=F
		module procedure DlinkToTenArray!Tensor(len)=link,if len<linklength(link) stop.Do not copy the index in the link
		module procedure DTenArrayTolink!link=Tensor(:),the data in the link will be empty,if there is no data in the Tensor
!													,it will not copy to the link
	end interface
!	other way to build a DTensor:T=DbuildTen(dimen,data),dimen is a vector or type(Dimension),data is any type
!**********************************************************
!	Other Function or Subroutine:
!
!		Dgenerate:generate a DTensor with random number
!
!		setDTensorName(T,TensorName,NameID):set a Name to a Tensor. TensorName is CHARACTER*1 the name, NameID is the Id of the tensor integer,optional
!
!		cleanDTensorName(T):clean the name in the Tensor
!
!		DRNTensorDim:input dimension [1,1,2,1,3,1,1,4,1],output dimenison [2,3,4],if input [1,1,1,1,1],ouput [1] without name
!
!		DTensorindex(T,indexID,indexname,NameID):output the order in the dimension, whose's name is[indexID,NameID(optional),indexname]
!
!		DzeroTen:generate a Tensor with all element being 0d0
!
!		DTensor1(dimen,num):return a DTensor with all the elememt num,dimen is dimension of
!			the DTensor
!
!		cleanDTensor: clean all the data in type(DTensor)
!
!		copydTensor:copy the Tensor data to another,do not allocate new data,it is faster that assignmentTen(=)
!
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
!		subroutine DproductTen(T,T1,T2):
!				productTen regard the last index of T1 and the first index
!	 		of T2 as the dimenion for matrix-product,other index will be see
!	 		as another dimenison.T1 and T2 can be any rank,but the last dimenion
!	 		of T1 and the first diemsion of T2 should be equal.
!			the output store in T,and do not allocate data in T,should input T
!	  		should be allocate before calling.
!				T%totaldata should be equal to T1*T2,the dimension and rank will be
!	 		modify
!
!		Tenproduct(T1_,T2_,i1_,i2_)
!			T1:[i1,i2,i3,i4,i5,i6,i7,i8]
!			T2:[j1,j2,j3,j4,j5,j6,j7,j8,j9,j10]
!			i1=(/5,1,2/)
!			i2=(10,3,5,6)
!			then the result will be T1'*T2'
!			where ,
!			T1'=[i3,i4,i6,i7,i8,(i5*i1*i2)]
!			T2'=[(j10*j3*j5*j6),j1,j2,j4,j7,j8,j9]
!			output T1'*T2
! 			input Tensor should be in its original dimenison,there is no contract on it
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
!		DsubTen(T,inde):inde=[-1,inde_min,inde_max] output data(inde_min:inde_max,:)
!			or[-2,inde_min,inde_max],data(:,inde_min:inde_max)
!			or [-1,inde_row] [-2,inde_col],output row or col
!			or [inde1_min,inde1_max,inde2_min,inde2_max] output data(inde1_min:inde1_max,inde2_min:inde2_max)
!
!		DQRlink:
!			make QR decompostion: A = QR, Q is an orthomomal matrix, and R is an upper triangle
!			The size of matrix A is M x N. the size of Q is M x min(M,N), R is  min(M,N) x N
!  		computes a QR factorization of a complex m by n matrix T  
!			A=(res.i.1)*(res.i.2)
!			ouput is a Tensorlink
!
!		DLQlink:
!			make LQ decompostion: A = LQ, Q is an orthomomal matrix, and L is a lower triangle
!			matrix. The size of matrix A is M x N. the size of L is M x min(M,N), Q is  min(M,N) x N
!			A=(res.i.2)*(res.i.1)
!			ouput is a Tensorlink
!
!		DQRdecompose
!			make QR decompostion: A = QR, Q is an orthomomal matrix, and R is an upper triangle
!			The size of matrix A is M x N. the size of Q is M x min(M,N), R is  min(M,N) x N
!  		computes a QR factorization of a complex m by n matrix T  
!			T=res(1)*res(2)
! 			the output is a array,if it is a array,T should be allocated
!
!
!		DLQdecompose:
!			make LQ decompostion: A = LQ, Q is an orthomomal matrix, and L is a lower triangle
!  		matrix. The size of matrix A is M x N. the size of L is M x min(M,N), Q is  min(M,N) x N
!			T=res(2)*res(1)
!			the output is a array,if it is a array,T should be allocated
!
!		svddecompose(T,Ncut_):SVD decomposion ,output a array,Ncut_ is optional
! 			T should be two dimension		
!			T=U*s*V^T
!			 on output,the order in the array is
! 			[U]->[s]->[V^T]
!			 T=res(1)*res(res(2))*res(3)
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
!		Dcombination
!			T1 :a [...,l,m,n] matrix
!			T2 :a [...,l,m,n] matrix
!			Dcombination(T1,T2):a [...,l,m,n,2] matrix
!			or 
!			T1 :a [...,m,n,l] matrix
!			T2 :a [...,m,n] matrix
!			Dcombination(T1,T2):a [...,m,n,l+1] matrix
!		Dcombinationrow(T1,T2)
!			T1 :a [l,m,n,...] matrix
!			T2 :a [l,m,n,...] matrix
!			Dcombinationrow(T1,T2):a [2,l,m,n,...] matrix
!			or 
!			T1 :a [l,m,n,...] matrix
!			T2 :a [m,n,...] matrix
!			Dcombinationrow(T1,T2):a [l+1,m,n,...] matrix	
!			
!		operator Dvalue:
!			Dvalue(F,Tr):return the Dvalue of <phi|F|phi>,where F is an operator,
!				F=<I',s1',s2',..|F|I,s1,s2,..>
!			Dnorm2(Tr_):return  <phi|phi>
!			Dnorm(Tr):return  sqrt(<phi|phi>)
!
!		DSortTen:Bubble Sort,only rank=1 is allowed
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
!		Dnullifylink:make the head and end of the link point to null and length=0
!		Dnullifylist:make the head and end of the list point to null and length=0
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
!		Some function about random number(in function.f90):
!			randomnumber(): out put a random number from 0~1
!			set_seed(idum):idum a seed of integer,set the seed using for randomnumber,if call set_seed(0) of
!								Do not call this subroutine, the seed will be set by randomly.
!			out_randomseed():output the seed of the computer. One should use it at the begining, and can repete
!								the result by call set_seed(idum) will the same idum at the begining.
!			out_and_set_seed():set a seed randomly and output the seed.
!
!
!
!		MPI function:
!			sent_DTensor(Ten1,Ten2,ID1,ID2,ierr):send the data of Ten1 in ID1 to Ten2 in ID2
!			BCAST_DTensor(Ten1,ID,ierr):BCAST The Ten1 in ID to every cpus
!			sent_DTensorlink(link1,link2,ID1,ID2,ierr):Send TensorLink,form link1 in ID1 to link2 in ID2
!			BCAST_DTensorlink(link1,ID,ierr):BCAST The link1 in ID to every CPUs
!
!			example
!				integer::ierr,proID,proNum
!				type(DTensor)::T1,T2
!				call mpi_init(ierr)
!				call mpi_comm_rank(mpi_comm_world,proID,ierr)
!				call mpi_comm_size(mpi_comm_world,proNum,ierr ) 
!				call sent_DTensor(T1,T2,0,1,ierr)  !T1 in cpu0 sent to cpu2 and store in T2
!				call BCAST_DTensor(T1,0,ierr) 	  !T1 in cpu 0 send to every cpus,store in T1 in other cpus
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
	subroutine DassignmentTen(T,T2)
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
	subroutine DassignmentTenArray(T,T2)
		type(DTensor),intent(inout) ::T(:)
		type(DTensor),intent(in) :: T2(:)
		integer::length,i
		length=size(T2)
		if(size(T).lt.length)then
			write(*,*)"ERROR in assignment of two DTensor array "
			write(*,*)"T1(:)=T2(:),size(T1)<size(T2)"
			write(*,*)size(T),length
			stop
		end if
		do i=1,length
			T(i)=T2(i)
		end do
		return
	end subroutine
	subroutine DlinkToTenArray(T,link)
		type(DTensor),intent(inout) ::T(:)
		type(DTensorlink),intent(in) :: link
		integer::length,i
		length=Dlinklength(link)
		if(size(T).lt.length)then
			write(*,*)"ERROR in assignment oflink to Tensor array "
			write(*,*)"T1(:)=link,size(T1)<linklength(link)"
			write(*,*)size(T),length
			stop
		end if
		do i=1,length
			T(i)=link.i.i
		end do
		return
	end subroutine
	subroutine DTenArrayToLink(link,T)
		type(DTensor),intent(in) ::T(:)
		type(DTensorlink),intent(inout) :: link
		integer::length,i
		length=size(T)
		call Dcleanlink(link)
		do i=1,length
			if(Dgetflag(T(i))) then
				call Dpush_back(link,T(i))
			end if
		end do
		return
	end subroutine	
!***********************************************	
!Do not allocate the Tensor data
!The output Tensor should be allocated
	subroutine copyDTensortwoDTensor(T,T2)
		type(DTensor),intent(inout) ::T
		type(DTensor),intent(in) :: T2
		integer::length
		if(T%TotalData.ne.T2%Totaldata) then
			write(*,*)"ERROR in copyTen"
			stop
		end if
		T%rank=T2%rank
		T%TenDim=T2%TenDim
		call dcopy(length,T2%DTensor_Data,1,T%DTensor_Data,1)
		return
	end subroutine
!Do not allocate the Tensor data
!copy the Tensor_Data to a vector	
	subroutine copyDTensordim1(Vec,T)
		real*8,intent(out) ::Vec(:)
		type(DTensor),intent(in) :: T
		integer::length
		length=T%TotalData
		if(size(Vec).lt.length) then
			write(*,*)"The array can not store the Data"
			write(*,*)"Length of the array",size(Vec)
			write(*,*)"Length of the Tensor",length
			stop
		end if
		call dcopy(length,T%DTensor_Data,1,Vec,1)
		return
	end subroutine
!if rank=2,then copy the Tensor_Data to a mat
!Do not allocate the Tensor data
! mat is a matrix	
	subroutine copyDTensordim2(Mat,T)
		real*8,intent(out) ::Mat(:,:)
		type(DTensor),intent(in) :: T
		integer::m,n
		if(DgetRank(T).ne.2)then
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
		call dcopy(m*n,T%DTensor_Data,1,Mat,1)
		return
	end subroutine
!if rank=3,then copy the Tensor_Data to a mat
!Do not allocate the Tensor data
!mat is a 3 dimension array	
	subroutine copyDTensordim3(Mat,T)
		real*8,intent(out) ::Mat(:,:,:)
		type(DTensor),intent(in) :: T
		integer::m,n,l
		if(DgetRank(T).ne.3)then
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
		call dcopy(m*n*l,T%DTensor_Data,1,Mat,1)
		return
	end subroutine	
!********************************************************
!allocate Tensor according to the dimen
	subroutine allocateDTensor1(T,dimen)
		type(DTensor),intent(inout) ::T
		type(dimension),intent(in)::dimen
		if(T%flag) then
			write(*,*)"T is already allocated"
			stop
		end if
		T%rank=Dimsize(dimen)
		T%TenDim=dimen
		T%totalData=outtotaldata(dimen)
		allocate(T%DTensor_Data(T%totalData))
		T%flag=.true.
		return
	end subroutine
	subroutine allocateDTensor2(T,dimen)
		type(DTensor),intent(inout) ::T
		integer,intent(in)::dimen(:)
		if(T%flag) then
			write(*,*)"T is already allocated"
			stop
		end if
		T%rank=size(dimen)
		T%TenDim=dimen
		T%totalData=product(dimen)
		allocate(T%DTensor_Data(T%totalData))
		T%flag=.true.
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
		real*8,allocatable,intent(inout) ::Vec(:)
		type(DTensor),intent(in) :: T
		integer::length
		length=T%TotalData
		if(allocated(Vec))then
			if(length.ne.size(Vec)) then
				deallocate(Vec)
				allocate(Vec(length))
			end if
		else
			allocate(Vec(length))
		end if
		call dcopy(length,T%DTensor_Data,1,Vec,1)
		return
	end subroutine
!if rank=2,then copy the DTensor_Data to a mat
! mat is a matrix	
	subroutine Dcopy_dim2(Mat,T)
		real*8,allocatable,intent(inout) ::Mat(:,:)
		type(DTensor),intent(in) :: T
		integer::m,n,m2,n2
		if(T%rank.eq.2) then
			m=T.dim.1
			n=T.dim.2
			if(allocated(Mat))then
				m2=size(MAt,1)
				n2=size(MAt,2)
				if((m.ne.m2).or.(n.ne.n2)) then
					deallocate(Mat)
					allocate(Mat(m,n))
				end if
			else
				allocate(Mat(m,n))
			end if
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
		real*8,allocatable,intent(inout) ::Mat(:,:,:)
		type(DTensor),intent(in) :: T
		integer::m,n,l,m2,n2,l2
		if(T%rank.eq.3) then
			m=T.dim.1
			n=T.dim.2
			l=T.dim.3
			if(allocated(Mat))then
				m2=size(MAt,1)
				n2=size(MAt,2)
				l2=size(MAt,3)
				if((m.ne.m2).or.(n.ne.n2).or.(l.ne.l2)) then
					deallocate(Mat)
					allocate(Mat(m,n,l))
				end if
			else
				allocate(Mat(m,n,l))
			end if
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
		if(allocated(T%DTensor_data))then
			if(length.ne.size(T%DTensor_Data)) then
				deallocate(T%DTensor_Data)
				allocate(T%DTensor_Data(length))
			end if
		else
			allocate(T%DTensor_Data(length))
		end if
		call dcopy(length,inData,1,T%DTensor_Data,1)
		return
	end subroutine
	subroutine DstoreTenData_dim2(T,inData)
		real*8,intent(in)::inData(:,:)
		type(DTensor),intent(inout)::T
		integer::length
		length=size(inData)
		if(allocated(T%DTensor_data))then
			if(length.ne.size(T%DTensor_Data)) then
				deallocate(T%DTensor_Data)
				allocate(T%DTensor_Data(length))
			end if
		else
			allocate(T%DTensor_Data(length))
		end if
		call dcopy(length,inData,1,T%DTensor_Data,1)
		return
	end subroutine
	subroutine DstoreTenData_dim3(T,inData)
		real*8,intent(in)::inData(:,:,:)
		type(DTensor),intent(inout)::T
		integer::length
		length=size(inData)
		if(allocated(T%DTensor_data))then
			if(length.ne.size(T%DTensor_Data)) then
				deallocate(T%DTensor_Data)
				allocate(T%DTensor_Data(length))
			end if
		else
			allocate(T%DTensor_Data(length))
		end if
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
	type(DTensor) function Dgenerate_NoName(Tdim) result (T)
		integer,intent(in) :: Tdim(:)
		integer :: rank,totalData,i
		real*8,allocatable :: DTensor_data(:)
		real*8 ::temp_real
		type(Dimension):: TenDim
		rank=size(Tdim)
		totalData=product(Tdim)
		allocate(DTensor_data(totalData))
		do i=1,totalData
			DTensor_data(i)=randomnumber()
		end do
		TenDim=Tdim
		call DstoreTenData(T,DTensor_data)
		T%rank=rank
		T%totalData=totalData
		T%TenDim=TenDim
		T%flag=.true.
		return
	end function
!		generate a Tensor with random number,set the name to the Tensor
	type(DTensor) function Dgenerate_Namechar(Tdim,w) result (T)
		integer,intent(in) :: Tdim(:)
		integer :: rank,totalData,i
		real*8,allocatable :: DTensor_data(:)
		real*8 ::temp_real
		character(len=*),intent(in)::w
		type(Dimension):: TenDim
		rank=size(Tdim)
		totalData=product(Tdim)
		allocate(DTensor_data(totalData))
		do i=1,totalData
			DTensor_data(i)=randomnumber()
		end do
		TenDim=Tdim
		call DimName(TenDim,w)
		call DstoreTenData(T,DTensor_data)
		T%rank=rank
		T%totalData=totalData
		T%TenDim=TenDim
		T%flag=.true.
		return
	end function
!*********************** TensorName   **********************
	subroutine setDTensorName2(T,TensorName)
		type(DTensor),intent(inout)::T
		CHARACTER(len=*),intent(in)::TensorName
		call DimName(T%TenDim,TensorName)
		return
	end subroutine
	subroutine cleanDTensorName(T)
		type(DTensor),intent(inout)::T
		call cleanDimensionName(T%TenDim)
		return
	end subroutine
	subroutine setDTensorName3(T,ith,w)
		type(DTensor),intent(inout)::T
		integer,intent(in)::ith
		CHARACTER(len=*),intent(in)::w
		call resetname(T%TenDim,ith,w)
		return
	end subroutine
	subroutine setDTensorName4(T,oldw,neww)
		type(DTensor),intent(inout)::T
		CHARACTER(len=*),intent(in)::oldw,neww
		call resetname(T%TenDim,oldw,neww)
		return
	end subroutine
	CHARACTER(len=len_of_Name+len_of_Name) function DoutIndexName(T,ith)
		type(DTensor),intent(in)::T
		integer,intent(in)::ith
		DoutIndexName=outName(T%TenDim,ith)
		return
	end function
	function DoutAllIndexName(T)
		CHARACTER(len=len_of_Name+len_of_Name),allocatable::DoutAllIndexName(:)
		type(DTensor),intent(in)::T
		integer::i
		allocate(DoutAllIndexName(T%rank))
		do i=1,T%rank
			DoutAllIndexName(i)=outName(T%TenDim,i)
		end do
		return
	end function
	CHARACTER(len=len_of_Name) function DoutTensorName(T,ith)
		type(DTensor),intent(in)::T
		integer,intent(in)::ith
		DoutTensorName=outNameTen(T%TenDim,ith)
		return
	end function
	function DoutAllTensorName(T)
		CHARACTER(len=len_of_Name),allocatable::DoutAllTensorName(:)
		type(DTensor),intent(in)::T
		integer::i
		allocate(DoutAllTensorName(T%rank))
		do i=1,T%rank
			DoutAllTensorName(i)=outNameTen(T%TenDim,i)
		end do
		return
	end function
!input dimension [1,1,2,1,3,1,1,4,1]
!output dimenison [2,3,4]	
	subroutine DRNTensorDim(T)
		type(DTensor),intent(inout)::T
		type(Dimension)::dimen
		dimen=RNDim(T%TenDim)
		T%TenDim=dimen
		T%rank=Dimsize(dimen)
		return
	end subroutine
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
	integer function DgetTenDim_Namei(T,w)
		type(DTensor),intent(in) :: T
		character(len=*),intent(in)::w
		integer :: inde
		inde=Nameorder(T%TenDim,w)
		if(inde.eq.0)then
			DgetTenDim_Namei=0
			return
		end if
		DgetTenDim_Namei=Dim_i(T%TenDim,inde)
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
				write(*,*) T%DTensor_data
			case (2)
				write(*,'(99999F12.8)') T%DTensor_data
			case (3)
				write(*,'(99999F10.2)') T%DTensor_data
			end 	select
			call Dprint0(T%TenDim)
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
	subroutine DTprint2_file(T,printtype,fileaddress,fileunit,replace)
		type(DTensor),intent(in) :: T
		integer,intent(in)::printtype
		CHARACTER*100,intent(in)::fileaddress
		logical,intent(in)::replace
		integer,intent(in)::fileunit
		if(replace)then
			open(unit=fileunit,file=fileaddress,STATUS='REPLACE',POSITION='APPEND')
		else
			open(unit=fileunit,file=fileaddress,STATUS='old',POSITION='APPEND')
		end if
		if(T%flag) then
			write(fileunit,*) "The rank of the DTensor is"
			write(fileunit,*) T%rank
			write(fileunit,*) "The number of  data of the DTensor is"
			write(fileunit,*) T%totalData
			write(fileunit,*) "The data of the DTensor is"
			select case (printtype)
			case (1)
				write(fileunit,*) T%DTensor_data
			case (2)
				write(fileunit,'(99999F12.8)') T%DTensor_data
			case (3)
				write(fileunit,'(99999F10.2)') T%DTensor_data
			end 	select
			!call Dprint0(T%TenDim)
			write(fileunit,*) "***end***"
			
		else
			write(fileunit,*) "There is no data"
		end if
		close(unit=fileunit)
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
						write(*,*) T%DTensor_data
					case (2)
						write(*,'(99999F12.8)') T%DTensor_data
					case (3)
						write(*,'(99999F10.2)') T%DTensor_data
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
							write (*,*) Tdata(i,:)
							write(*,*) ""
						end do 
					case (2)
						do i=1,Tdim(1)
							write (*,'(99999F12.8)') Tdata(i,:)
							write(*,*) ""
						end do 
					case (3)
						do i=1,Tdim(1)
							write (*,'(99999F10.2)') Tdata(i,:)
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
								write (*,*) Tdata3(i,:,j)
								write(*,*) ""
							end do
						end do 
					case (2)
						do j=1,Tdim(3)
							write(*,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write (*,'(99999F12.8)') Tdata3(i,:,j)
								write(*,*) ""
							end do
						end do 
					case (3)
						do j=1,Tdim(3)
							write(*,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write (*,'(99999F10.2)') Tdata3(i,:,j)
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
								write (*,*) Tdata4(i,:,k,j)
								write(*,*) ""
								end do
							end do 
						end do
					case (2)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(*,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write (*,'(99999F12.8)') Tdata4(i,:,k,j)
								write(*,*) ""
								end do
							end do 
						end do
					case (3)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(*,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write (*,'(99999F10.2)') Tdata4(i,:,k,j)
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

	subroutine DTMprint2_file(T,printtype,fileaddress,fileunit,replace)
		type(DTensor),intent(in) :: T
		integer,intent(in)::printtype
		CHARACTER*100,intent(in)::fileaddress
		logical,intent(in)::replace
		integer,intent(in)::fileunit
		real*8,allocatable :: Tdata(:,:)
		real*8,allocatable :: Tdata3(:,:,:)
		real*8,allocatable :: Tdata4(:,:,:,:)
		integer,allocatable ::Tdim(:)
		integer :: i,j,k
		if(replace)then
			open(unit=fileunit,file=fileaddress,STATUS='REPLACE',POSITION='APPEND')
		else
			open(unit=fileunit,file=fileaddress,STATUS='old',POSITION='APPEND')
		end if
		if(T%flag) then!if 1
			write(fileunit,*) "*** START ***"
			select case(T%rank)
				case(1)
					select case (printtype)
					case (1)
						write(fileunit,*) T%DTensor_data
					case (2)
						write(fileunit,'(99999F12.8)') T%DTensor_data
					case (3)
						write(fileunit,'(99999F10.2)') T%DTensor_data
					end 	select
					write(fileunit,*) "*** END ***"
				case(2)
					allocate(Tdim(2))
					Tdim=T%TenDim
					allocate(Tdata(Tdim(1),Tdim(2)))
					Tdata=reshape(T%DTensor_Data,shape(Tdata))
					select case (printtype)
					case (1)
						do i=1,Tdim(1)
							write(fileunit,*) Tdata(i,:)
							
						end do 
					case (2)
						do i=1,Tdim(1)
							write(fileunit,'(99999F12.8)') Tdata(i,:)
							
						end do 
					case (3)
						do i=1,Tdim(1)
							write(fileunit,'(99999F10.2)') Tdata(i,:)
							
						end do 
					end 	select
						write(fileunit,*) "*** END ***"
				case(3)
					allocate(Tdim(3))
					Tdim=T%TenDim
					allocate(Tdata3(Tdim(1),Tdim(2),Tdim(3)))
					Tdata3=reshape(T%DTensor_Data,shape(Tdata3))
					select case (printtype)
					case (1)
						do j=1,Tdim(3)
							write(fileunit,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write(fileunit,*) Tdata3(i,:,j)
								
							end do
						end do 
					case (2)
						do j=1,Tdim(3)
							write(fileunit,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write(fileunit,'(99999F12.8)') Tdata3(i,:,j)
								
							end do
						end do 
					case (3)
						do j=1,Tdim(3)
							write(fileunit,*) "(*,*,",j,")"
							do i=1,Tdim(1)
								write(fileunit,'(99999F10.2)') Tdata3(i,:,j)
								
							end do
						end do 
					end 	select
					write(fileunit,*) "*** END ***"
				case(4)
					allocate(Tdim(4))
					Tdim=T%TenDim
					allocate(Tdata4(Tdim(1),Tdim(2),Tdim(3),Tdim(4)))
					Tdata4=reshape(T%DTensor_Data,shape(Tdata4))
					select case (printtype)
					case (1)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(fileunit,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write(fileunit,*) Tdata4(i,:,k,j)
								
								end do
							end do 
						end do
					case (2)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(fileunit,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write(fileunit,'(99999F12.8)') Tdata4(i,:,k,j)
								
								end do
							end do 
						end do
					case (3)
						do j=1,Tdim(4)
							do k=1,Tdim(3)
								write(fileunit,*) "(*,*,",k,j,")"
								do i=1,Tdim(1)
								write(fileunit,'(99999F10.2)') Tdata4(i,:,k,j)
								
								end do
							end do 
						end do
					end 	select
					write(fileunit,*) "*** END ***"
				case default
					call DTprint2_file(T,printtype,fileaddress,fileunit,replace)
			end select
		else!if 1
			write(fileunit,*) "There is no data"
		end if!if 1
		close(unit=fileunit)
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
		if(T%rank.eq.1) then
			Dpermute_rank2=T
			return
		end if
		if(T%rank.gt.2) then
			write(*,*)"ERROR in Dpermute_rank2"
			write(*,*)"stop"
			stop
			return
		end if
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
		integer,allocatable ::inde(:)
		integer::lenOrder,i,totaldata,j
		type(Dimension)::dimen		
		lenorder=size(newOrder)-1
		allocate(inde(lenorder))
		totaldata=T%totalData
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		allocate(Dpermutation%DTensor_data(totaldata))
		call DCOPY(totaldata,T%DTensor_data,1,Dpermutation%DTensor_data,1)
		dimen=T%TenDim
		do i=lenorder,1,-1
			call Dpermutefo_data(Dpermutation%DTensor_data,inde(i),dimen,totaldata)
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		call Dresetdim(Dpermutation,dimen)
		Dpermutation%flag=.true.
		Dpermutation%totalData=totaldata
		return
	end function
	type(DTensor) function Dpermutation_name(T,newOrderchar)result(Dpermutation)
		type(DTensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::newOrderchar(:)
		integer,allocatable::newOrder(:)
		integer,allocatable ::inde(:)
		integer::lenOrder,i,totaldata,j
		type(Dimension)::dimen		
		allocate(newOrder(size(newOrderchar)))
		newOrder=Nameorder(T%TenDim,newOrderchar)
		lenorder=size(newOrder)-1
		allocate(inde(lenorder))
		totaldata=T%totalData
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		allocate(Dpermutation%DTensor_data(totaldata))
		call DCOPY(totaldata,T%DTensor_data,1,Dpermutation%DTensor_data,1)
		dimen=T%TenDim
		do i=lenorder,1,-1
			call Dpermutefo_data(Dpermutation%DTensor_data,inde(i),dimen,totaldata)
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		call Dresetdim(Dpermutation,dimen)
		Dpermutation%flag=.true.
		Dpermutation%totalData=totaldata
		return
	end function
	subroutine Dpermutation_data3(T1data,index_not_permute,dimens,totaldata)
		real*8,intent(inout) :: T1data(:)
		integer,intent(in) ::   index_not_permute
		type(Dimension),intent(inout) ::dimens
		integer,allocatable::dimen(:)
		integer ::i,newdimen(3),totaldata
		real*8,allocatable :: newdata(:,:,:),Tdata(:,:,:)
		dimen=dimens
		allocate(Tdata(dimen(1),dimen(2),dimen(3)))
		call DCOPY(totaldata,T1data,1,Tdata,1)
		if(index_not_permute.eq.1) then
			newdimen(1)=dimen(1)
			newdimen(2)=dimen(3)
			newdimen(3)=dimen(2)
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(1)
				newdata(i,:,:)=transpose(Tdata(i,:,:))
			end do
			call DCOPY(totaldata,newdata,1,T1data,1)
			dimens=Dimpermute(dimens,(/1,3,2/))
			return
		end if
		if(index_not_permute.eq.2) then
			newdimen(1)=dimen(3)
			newdimen(2)=dimen(2)
			newdimen(3)=dimen(1)
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(2)
				newdata(:,i,:)=transpose(Tdata(:,i,:))
			end do
			call DCOPY(totaldata,newdata,1,T1data,1)
			dimens=Dimpermute(dimens,(/3,2,1/))
			return
		end if
		if(index_not_permute.eq.3) then
			newdimen(1)=dimen(2)
			newdimen(2)=dimen(1)
			newdimen(3)=dimen(3)
			allocate(newdata(newdimen(1),newdimen(2),newdimen(3)))
			do i=1,dimen(3)
				newdata(:,:,i)=transpose(Tdata(:,:,i))
			end do
			call DCOPY(totaldata,newdata,1,T1data,1)
			dimens=Dimpermute(dimens,(/2,1,3/))
			return
		end if
	 end subroutine
	subroutine Dpermutation_data2 (T1data,dimens,totaldata)
		real*8,intent(inout) :: T1data(:)
		type(Dimension),intent(inout) ::dimens
		integer,intent(inout)::totaldata
		real*8,allocatable::Tdata(:,:),newdata(:,:)
		integer,allocatable::dimen(:)
		if(Dimsize(dimens).eq.1) then
			return
		end if
		dimen=dimens
		allocate(Tdata(dimen(1),dimen(2)))
		call DCOPY(totaldata,T1data,1,Tdata,1)
		allocate(newdata(dimen(2),dimen(1)))
		newdata=transpose(Tdata)
		call DCOPY(totaldata,newdata,1,T1data,1)
		dimens=Dimpermute(dimens,(/2,1/))
	end subroutine
	subroutine Dpermutefo_data(Tdata,inde,dimen,totaldata)
		real*8,allocatable,intent(inout)::Tdata(:)
		type(Dimension),intent(inout) ::dimen
		integer,intent(inout)::totaldata
		integer,intent(in)::inde
		integer::rank,num
		integer::oper(2,3)
		rank=Dimsize(dimen)
		if(inde.gt.rank) then
			write(*,*)"ERROR IN permutefo"
			stop
		end if
		if(inde.eq.1) then
			return
		end if
		if(inde.eq.rank) then
			dimen=DimConstract(dimen,1,rank-2)
			call Dpermutation_data2(Tdata,dimen,totaldata)
			dimen=DimDecomposeAll(dimen)
			return
		end if
		num=inde-2
		dimen=DimConstract(dimen,1,num)
		num=rank-3
		dimen=DimConstract(dimen,3,num)
		call Dpermutation_data3(Tdata,3,dimen,totaldata)
		dimen=DimDecomposeAll(dimen)
		return
	end subroutine
	subroutine Dpermuteback_data(Tdata,inde,dimen,totaldata)
		real*8,allocatable,intent(inout)::Tdata(:)
		type(Dimension),intent(inout) ::dimen
		integer,intent(inout)::totaldata
		integer,intent(in)::inde
		integer::rank,num
		integer::oper(2,3)
		rank=Dimsize(dimen)
		if(inde.eq.rank) then
			return
		end if
		if(inde.eq.1) then
			!permuteback=contracts(T,2,rank)
			!permuteback=.p.permuteback
			!call dimoperation(permuteback,(/3/))
			dimen=DimConstract(dimen,2,rank-2)
			call Dpermutation_data2(Tdata,dimen,totaldata)
			dimen=DimDecomposeAll(dimen)
			return
		end if
		num=inde-2
		dimen=DimConstract(dimen,1,num)
		num=rank-3
		dimen=DimConstract(dimen,3,num)
		call Dpermutation_data3(Tdata,1,dimen,totaldata)
		dimen=DimDecomposeAll(dimen)
		
		!oper(1,:)=(/1,1,inde-1/)
		!oper(2,:)=(/1,3,rank/)
		!permuteback=T.cd.oper
		!permuteback=permuteback.p.1
		!call dimoperation(permuteback,(/3/))
		return
	end subroutine

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
	type(DTensor) function Dpermutefo_name(T,indechar)result(Dpermutefo)
		type(DTensor),intent(in)::T
		character(len=*),intent(in)::indechar
		integer::inde
		integer::rank
		integer::oper(2,3)
		rank=DgetRank(T)
		inde=Nameorder(T%TenDim,indechar)
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
!*****************permutefo******************************
!		T_{1,2,3,..,j,..,i,.,k,...,n},permutefo(T,(/i,j,k/))=_{i,j,k,1,2,3,...,n}
!
	type(DTensor) function Dpermutefo_vec(T,vec_)
		type(DTensor),intent(in)::T
		integer,intent(in)::vec_(:)
		integer::rank,lenVec,totalData,i,j
		integer::oper(2,3)
		integer,allocatable::vec(:)
		type(Dimension)::dimen		
		rank=DgetRank(T)
		lenVec=size(vec_)
		totalData=T%totalData
		dimen=T%TenDim
		allocate(vec(lenVec))
		vec=vec_
		allocate(Dpermutefo_vec%DTensor_Data(totalData))
		call DCOPY(totaldata,T%DTensor_data,1,Dpermutefo_vec%DTensor_Data,1)
		do i=lenVec,1,-1
			call Dpermutefo_data(Dpermutefo_vec%DTensor_data,vec(i),dimen,totaldata)
			do j=1,i-1
				if(vec(j).lt.vec(i))then
					vec(j)=vec(j)+1
				end if
			end do
		end do
		call Dresetdim(Dpermutefo_vec,dimen)
		Dpermutefo_vec%flag=.true.
		Dpermutefo_vec%totalData=totaldata
		return
	end function
	type(DTensor) function Dpermutefo_vec_name(T,indechar)result(Dpermutefo_vec)
		type(DTensor),intent(in)::T
		character(len=*),intent(in)::indechar(:)
		integer::rank,lenVec,totalData,i,j
		integer::oper(2,3)
		integer,allocatable::vec(:)
		type(Dimension)::dimen		
		rank=DgetRank(T)
		lenVec=size(indechar)
		totalData=T%totalData
		dimen=T%TenDim
		allocate(vec(lenVec))
		vec=Nameorder(T%TenDim,indechar)
		allocate(Dpermutefo_vec%DTensor_Data(totalData))
		call DCOPY(totaldata,T%DTensor_data,1,Dpermutefo_vec%DTensor_Data,1)
		do i=lenVec,1,-1
			call Dpermutefo_data(Dpermutefo_vec%DTensor_data,vec(i),dimen,totaldata)
			do j=1,i-1
				if(vec(j).lt.vec(i))then
					vec(j)=vec(j)+1
				end if
			end do
		end do
		call Dresetdim(Dpermutefo_vec,dimen)
		Dpermutefo_vec%flag=.true.
		Dpermutefo_vec%totalData=totaldata
		return
	end function
!*****************permuteback******************************
!		T_{1,2,3,..,i,..,n},permuteback(T,i)=_{1,2,3,..,i-1,i+1,..,n,i}
!
	type(DTensor) function Dpermuteback(T,inde)
		type(DTensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		integer::oper(2,3)
		rank=DgetRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function Dpermuteback"
			write(*,*)"stop"
			stop
		end if
		if(inde.eq.rank) then
			Dpermuteback=.dc.T
			return
		end if
		if(inde.eq.1) then
			Dpermuteback=Dcontracts(T,2,rank)
			Dpermuteback=.p.Dpermuteback
			call Ddimoperation(Dpermuteback,(/3/))
			return
		end if
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,3,rank/)
		Dpermuteback=T.cd.oper
		Dpermuteback=Dpermuteback.p.1
		call Ddimoperation(Dpermuteback,(/3/))
		return
	end function
	type(DTensor) function Dpermuteback_name(T,indechar)result(Dpermuteback)
		type(DTensor),intent(in)::T
		character(len=*),intent(in)::indechar
		integer::inde
		integer::rank
		integer::oper(2,3)
		rank=DgetRank(T)
		inde=Nameorder(T%TenDim,indechar)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function Dpermuteback"
			write(*,*)"stop"
			stop
		end if
		if(inde.eq.rank) then
			Dpermuteback=.dc.T
			return
		end if
		if(inde.eq.1) then
			Dpermuteback=Dcontracts(T,2,rank)
			Dpermuteback=.p.Dpermuteback
			call Ddimoperation(Dpermuteback,(/3/))
			return
		end if
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,3,rank/)
		Dpermuteback=T.cd.oper
		Dpermuteback=Dpermuteback.p.1
		call Ddimoperation(Dpermuteback,(/3/))
		return
	end function
!*****************permuteback******************************
!		T_{1,2,3,..,i,..,n},permuteback(T,i)=_{1,2,3,..,i-1,i+1,..,n,i}
!
	type(DTensor) function Dpermuteback_vec(T,vec_)
		type(DTensor),intent(in)::T
		integer,intent(in)::vec_(:)
		integer::rank,lenVec,totalData,i,j
		integer::oper(2,3)
		integer,allocatable::vec(:)
		type(Dimension)::dimen		
		rank=DgetRank(T)
		lenVec=size(vec_)
		totalData=T%totalData
		dimen=T%TenDim
		allocate(vec(lenVec))
		vec=vec_
		allocate(Dpermuteback_vec%DTensor_Data(totalData))
		call DCOPY(totaldata,T%DTensor_data,1,Dpermuteback_vec%DTensor_Data,1)
		do i=1,lenVec
			call Dpermuteback_data(Dpermuteback_vec%DTensor_data,vec(i),dimen,totaldata)
			do j=lenVec,i+1,-1
				if(vec(j).gt.vec(i))then
					vec(j)=vec(j)-1
				end if
			end do
		end do
		call Dresetdim(Dpermuteback_vec,dimen)
		Dpermuteback_vec%flag=.true.
		Dpermuteback_vec%totalData=totaldata
		return
	end function
	type(DTensor) function Dpermuteback_vec_name(T,indechar)result(Dpermuteback_vec)
		type(DTensor),intent(in)::T
		character(len=*),intent(in)::indechar(:)
		integer::rank,lenVec,totalData,i,j
		integer::oper(2,3)
		integer,allocatable::vec(:)
		type(Dimension)::dimen		
		rank=DgetRank(T)
		lenVec=size(indechar)
		totalData=T%totalData
		dimen=T%TenDim
		allocate(vec(lenVec))
		vec=Nameorder(T%TenDim,indechar)
		allocate(Dpermuteback_vec%DTensor_Data(totalData))
		call DCOPY(totaldata,T%DTensor_data,1,Dpermuteback_vec%DTensor_Data,1)
		do i=1,lenVec
			call Dpermuteback_data(Dpermuteback_vec%DTensor_data,vec(i),dimen,totaldata)
			do j=lenVec,i+1,-1
				if(vec(j).gt.vec(i))then
					vec(j)=vec(j)-1
				end if
			end do
		end do
		call Dresetdim(Dpermuteback_vec,dimen)
		Dpermuteback_vec%flag=.true.
		Dpermuteback_vec%totalData=totaldata
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
!****************** permutebackInde**************************
!		T_{1,2,3,..,i,..,n},permutebackInde(T,i)=_{1,2,3,..,i-1,n,i,i+1,..,n-1}
!
	type(DTensor) function DpermutebackInde(T,inde)
		type(DTensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		integer::oper(2,3)
		rank=DgetRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function DpermutebackInde",inde,rank
			write(*,*)"stop"
			stop
		end if
		if(inde.eq.rank) then
			DpermutebackInde=.dc.T
			return
		end if
		if(inde.eq.1) then
			DpermutebackInde=Dcontracts(T,1,rank-1)
			DpermutebackInde=.p.DpermutebackInde
			call Ddimoperation(DpermutebackInde,(/3/))
			return
		end if
		oper(1,:)=(/1,1,inde-1/)
		oper(2,:)=(/1,2,rank-inde+1/)
		DpermutebackInde=T.cd.oper
		DpermutebackInde=DpermutebackInde.p.1
		call Ddimoperation(DpermutebackInde,(/3/))
		return
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
	type(DTensor) function DTenproduct_noName(T1_,i1,T2_,i2) result(T)
		type(DTensor),intent(in) :: T1_,T2_
		integer,intent(in) :: i1(:),i2(:)
		type(DTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.(if_original_dim(T1_%TenDim).and.if_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		call Ddimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		call Ddimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(DTensor) function DTenproduct_name(T1_,name1,T2_,name2) result(T)
		type(DTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:),name2(:)
		integer :: i1(size(name1)),i2(size(name2))
		type(DTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.(if_original_dim(T1_%TenDim).and.if_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		i1=Nameorder(T1_%TenDim,name1)
		i2=Nameorder(T2_%TenDim,name2)
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		call Ddimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		call Ddimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(DTensor) function DTenproduct_name_rename1(T1_,name1,T2_,name2,rename,whichname) result(T)
		type(DTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:),name2(:),rename(2)
		integer,intent(in)::whichname
		integer :: i1(size(name1)),i2(size(name2))
		type(DTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		character*50::oldName,newName
		if(.not.(if_original_dim(T1_%TenDim).and.if_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		oldName=rename(1)
		newName=rename(2)
		i1=Nameorder(T1_%TenDim,name1)
		i2=Nameorder(T2_%TenDim,name2)
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		if(whichname.eq.1)then
			call resetname(T1%TenDim,oldName,newName)
		end if
		call Ddimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		if(whichname.eq.2)then
			call resetname(T2%TenDim,oldName,newName)
		end if
		call Ddimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(DTensor) function DTenproduct_name_rename2(T1_,name1,rename1,T2_,name2,rename2) result(T)
		type(DTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:),name2(:),rename1(2),rename2(2)
		integer :: i1(size(name1)),i2(size(name2))
		type(DTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		character*50::oldName,newName
		if(.not.(if_original_dim(T1_%TenDim).and.if_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		i1=Nameorder(T1_%TenDim,name1)
		i2=Nameorder(T2_%TenDim,name2)
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		oldName=rename1(1)
		newName=rename1(2)
		call resetname(T1%TenDim,oldName,newName)
		call Ddimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		oldName=rename2(1)
		newName=rename2(2)
		call resetname(T2%TenDim,oldName,newName)
		call Ddimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(DTensor) function DTenproduct_int_name(T1_,i1,T2_,name2) result(T)
		type(DTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name2(:)
		integer,intent(in)::i1(:)
		integer :: i2(size(name2))
		type(DTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.(if_original_dim(T1_%TenDim).and.if_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		i2=Nameorder(T2_%TenDim,name2)
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		call Ddimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		call Ddimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(DTensor) function DTenproduct_name_int(T1_,name1,T2_,i2) result(T)
		type(DTensor),intent(in) :: T1_,T2_
		character(len=*),intent(in)::name1(:)
		integer,intent(in)::i2(:)
		integer :: i1(size(name1))
		type(DTensor) :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.(if_original_dim(T1_%TenDim).and.if_original_dim(T2_%TenDim))) then
			write(*,*)"ERROR in Tenproduct"
			write(*,*)"stop"
			stop
		end if
		i1=Nameorder(T1_%TenDim,name1)
		rank1=T1_%rank
		rank2=T2_%rank
		leni1=size(i1)
		leni2=size(i2)
		T1=T1_.pb.i1
		call Ddimoperation(T1,(/1,rank1-leni1+1,rank1/))
		T2=T2_.pf.i2
		call Ddimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
		return
	end function
	type(DTensor) function DTenproduct_old(T1_,T2_,i1_,i2_) result(T)
		type(DTensor),intent(in) :: T1_,T2_
		integer,intent(in) :: i1_(:),i2_(:)
		integer :: i,j,k,D2,decompoD1,decompoD2!Di use for decompese
		type(DTensor) :: T1,T2
		integer :: decompoD1keep,decompoD2keep,rank1,&
			Tdim1(T1_%rank),Tdim2(T2_%rank),rank2,oper(2,3)
		integer::leni1,leni2,i1(size(i1_)),i2(size(i2_))
		if(.not.(if_original_dim(T1_%TenDim).and.if_original_dim(T2_%TenDim))) then
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
		call Ddimoperation(T1,(/1,rank1-leni1+1,rank1/))
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
		call Ddimoperation(T2,(/1,1,leni2/))
		T=T1 * T2
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
		if(.not.T1%flag) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the first Tensor"
			write(*,*)"stop"
			stop
		end if
		if(.not.T2%flag) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the second Tensor"
			write(*,*)"stop"
			stop
		end if	
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
				if((D1.i.2).ne.(D2.i.1)) then
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
!ccccccccccccccccc   productTen   ccccccccccccccccccc
!		productTen regard the last index of T1 and the first index
!	  of T2 as the dimenion for matrix-product,other index will be see
!	  as another dimenison.T1 and T2 can be any rank,but the last dimenion
!	  of T1 and the first diemsion of T2 should be equal.
!		the output store in T,and do not allocate data in T,should input T
!	  should be allocate before calling.
!		T%totaldata should be equal to T1*T2,the dimension and rank will be
!	  modify	
	subroutine DproductTen(T,T1,T2)
		type(DTensor),intent(inout)::T
		type(DTensor),intent(in) :: T1,T2
		integer::rank1,rank2,flag,T1m,T1n,T2m,T2n,T1l,T2l
		type(Dimension)::D1,D2,newD
		real*8::Ndata
		integer::i
		if(.not.T1%flag) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the first Tensor"
			write(*,*)"stop"
			stop
		end if
		if(.not.T2%flag) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the second Tensor"
			write(*,*)"stop"
			stop
		end if	
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
			write(*,*)"ERROR in DproductTen",rank1,rank2
			stop
		end if
		select case (flag)
			case (1)
				T1m=T1.dim.1
				T2n=T2.dim.1
				if(T1m.ne.T2n) then
					write(*,*)"ERROR in DproductTen,case 1,stop"
					stop
				end if
				call dresetdim(T,(/1/))
				T%DTensor_Data=Ddot(T1m,T1%DTensor_Data,1,T2%DTensor_Data,1)
				return
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
					write(*,*)"ERROR in DproductTen,case 2,stop"
					stop
				end if
				call dresetdim(T,newD)
				call DGEMV('T',T2m,T2n,1D0,T2%DTensor_Data,T2m,T1%DTensor_Data,&
							1,0D0,T%DTensor_Data,1)
				RETURN
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
					write(*,*)"ERROR in DproductTen,case 3,stop"
					stop
				end if
				call dresetdim(T,newD)
				call DGEMV('N',T1m,T1n,1D0,T1%DTensor_Data,T1m,T2%DTensor_Data,&
						1,0D0,T%DTensor_Data,1)
				return
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
				if((D1.i.2).ne.(D2.i.1)) then
					write(*,*)"error DproductTen,dimension"
					stop
				end if
				T1m=D1.i.1
				T1n=D1.i.2
				T2m=D2.i.1
				T2n=D2.i.2
				call dresetdim(T,newD)
				call DGEMM('N','N',T1m,T2n,T1n,1D0,T1%DTensor_Data,T1m,&
					T2%DTensor_Data,T2m,0D0,T%DTensor_data,T1m)
				return
		case default 
			write(*,*) "ERROR in DproductTen,no such data"
			stop
		end 	select
		return
	end subroutine	
	
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
!********************************************************************
!make QR decompostion: A = QR, Q is an orthomomal matrix, and R is an upper triangle
!The size of matrix A is M x N. the size of Q is M x min(M,N), R is  min(M,N) x N
!  computes a QR factorization of a complex m by n matrix T  
!	T=(res.i.1)*(res.i.2)
!********************************************************************

	type(DTensorlink) function DQRlink(T) result(res)
		type(DTensor),intent(in)::T
		real*8,allocatable::  work(:),tau(:),Tdata(:,:)
		Type(DTensor)::v,vv,identity
		type(DTensornode),pointer::respointer
		type(Dimension)::dimen
		integer :: i,j,M,N,INFO,min_MN
		M = T.dim.1
		N = T.dim.2
		if(DgetRank(T).ne.2) then
			write(*,*)"ERROR in DQRdecomposition"
			write(*,*)"input Tensor should be a matrix"
			stop
		endif
		Tdata=T
		min_MN=min(M,N)
		allocate(tau(min_MN))
		allocate(work(N))
		call DGEQR2( M, N, Tdata, M, tau, WORK, INFO )
		if(info.ne.0) then
				write (*,*) "Error in DQRdecomposition ,info=",info
				stop
		end if
		!compute Q
		v=DzeroTen((/M/))
		identity=Deye(M,M)
		allocate(respointer)
		do i=1,min_MN
			if(i.ne.1)then
				call Dmodify(v,i-1,0d0)
			end if
			call Dmodify(v,i,1d0)
			call Dmodify(v,(/i+1,M/),Tdata(i+1:M,i))
			vv=tau(i)*(v.xx.v)
			if(i.eq.1)then
				respointer%Ten=identity-vv
			else
				respointer%Ten=(respointer%Ten)*(identity-vv)
			end if
		end do
		respointer%Ten=DsubTen(respointer%Ten,(/-2,1,min_MN/))
		dimen=(T%TenDim.sub.1)+(/min_MN/)
		respointer%Ten%TenDim=dimen
		call Dpush_back(res,respointer)
		!compute R
		nullify(respointer)
		allocate(respointer)
		do i=1,min_MN
			do j=1,i-1
				Tdata(i,j)=0d0
			end do
		end do
		respointer%Ten=Tdata(1:min_MN,:)
		dimen=(/min_MN/)+(T%TenDim.sub.2)
		respointer%Ten%TenDim=dimen
		call Dpush_back(res,respointer)
		return
	end function
!********************************************************************
!make QR decompostion: A = QR, Q is an orthomomal matrix, and R is an upper triangle
!The size of matrix A is M x N. the size of Q is M x min(M,N), R is  min(M,N) x N
!  computes a QR factorization of a complex m by n matrix T  
!	T=res(1)*res(2),the output is a array of DTensor
!********************************************************************

	function DQRdecompose(T) result(res)
		type(DTensor),allocatable::res(:)
		type(DTensor),intent(in)::T
		real*8,allocatable::  work(:),tau(:),Tdata(:,:)
		Type(DTensor)::v,vv,identity
		type(Dimension)::dimen
		integer :: i,j,M,N,INFO,min_MN
		M = T.dim.1
		N = T.dim.2
		if(DgetRank(T).ne.2) then
			write(*,*)"ERROR in DQRdecomposition"
			write(*,*)"input Tensor should be a matrix"
			stop
		endif
		Tdata=T
		min_MN=min(M,N)
		allocate(tau(min_MN))
		allocate(work(N))
		call DGEQR2( M, N, Tdata, M, tau, WORK, INFO )
		if(info.ne.0) then
				write (*,*) "Error in DQRdecomposition ,info=",info
				stop
		end if
		!compute Q
		v=DzeroTen((/M/))
		identity=Deye(M,M)
		allocate(res(2))
		do i=1,min_MN
			if(i.ne.1)then
				call Dmodify(v,i-1,0d0)
			end if
			call Dmodify(v,i,1d0)
			call Dmodify(v,(/i+1,M/),Tdata(i+1:M,i))
			vv=tau(i)*(v.xx.v)
			if(i.eq.1)then
				res(1)=identity-vv
			else
				res(1)=res(1)*(identity-vv)
			end if
		end do
		res(1)=DsubTen(res(1),(/-2,1,min_MN/))
		dimen=(T%TenDim.sub.1)+(/min_MN/)
		res(1)%TenDim=dimen
		!compute R
		do i=1,min_MN
			do j=1,i-1
				Tdata(i,j)=0d0
			end do
		end do
		res(2)=Tdata(1:min_MN,:)
		dimen=(/min_MN/)+(T%TenDim.sub.2)
		res(2)%TenDim=dimen
		return
	end function
	subroutine DQRdecomposition1(T,Q,R) 
		type(DTensor),intent(in)::T
		type(DTensor),intent(inout)::Q,R
		real*8,allocatable::  work(:),tau(:),Tdata(:,:)
		Type(DTensor)::v,vv,identity
		type(Dimension)::dimen
		integer :: i,j,M,N,INFO,min_MN
		M = T.dim.1
		N = T.dim.2
		if(DgetRank(T).ne.2) then
			write(*,*)"ERROR in DQRdecomposition"
			write(*,*)"input Tensor should be a matrix"
			stop
		endif
		Tdata=T
		min_MN=min(M,N)
		allocate(tau(min_MN))
		allocate(work(N))
		call DGEQR2( M, N, Tdata, M, tau, WORK, INFO )
		if(info.ne.0) then
				write (*,*) "Error in DQRdecomposition ,info=",info
				stop
		end if
		!compute Q
		v=DzeroTen((/M/))
		identity=Deye(M,M)
		do i=1,min_MN
			if(i.ne.1)then
				call Dmodify(v,i-1,0d0)
			end if
			call Dmodify(v,i,1d0)
			call Dmodify(v,(/i+1,M/),Tdata(i+1:M,i))
			vv=tau(i)*(v.xx.v)
			if(i.eq.1)then
				Q=identity-vv
			else
				Q=Q*(identity-vv)
			end if
		end do
		Q=DsubTen(Q,(/-2,1,min_MN/))
		dimen=(T%TenDim.sub.1)+(/min_MN/)
		Q%TenDim=dimen
		!compute R
		do i=1,min_MN
			do j=1,i-1
				Tdata(i,j)=0d0
			end do
		end do
		R=Tdata(1:min_MN,:)
		dimen=(/min_MN/)+(T%TenDim.sub.2)
		R%TenDim=dimen
		return
	end subroutine
	
! Q store in T
	subroutine DQRdecomposition2(T,R) 
		type(DTensor),intent(inout)::T,R
		real*8,allocatable::  work(:),tau(:),Tdata(:,:)
		Type(DTensor)::v,vv,identity
		type(Dimension)::dimen,tempdim
		integer :: i,j,M,N,INFO,min_MN
		M = T.dim.1
		N = T.dim.2
		if(DgetRank(T).ne.2) then
			write(*,*)"ERROR in DQRdecomposition"
			write(*,*)"input Tensor should be a matrix"
			stop
		endif
		Tdata=T
		tempdim=T%TenDim
		min_MN=min(M,N)
		allocate(tau(min_MN))
		allocate(work(N))
		call cleanDTensor(T)
		call DGEQR2( M, N, Tdata, M, tau, WORK, INFO )
		if(info.ne.0) then
				write (*,*) "Error in DQRdecomposition ,info=",info
				stop
		end if
		!compute Q
		v=DzeroTen((/M/))
		identity=Deye(M,M)
		do i=1,min_MN
			if(i.ne.1)then
				call Dmodify(v,i-1,0d0)
			end if
			call Dmodify(v,i,1d0)
			call Dmodify(v,(/i+1,M/),Tdata(i+1:M,i))
			vv=tau(i)*(v.xx.v)
			if(i.eq.1)then
				T=identity-vv
			else
				T=T*(identity-vv)
			end if
		end do
		T=DsubTen(T,(/-2,1,min_MN/))
		dimen=(tempdim.sub.1)+(/min_MN/)
		T%TenDim=dimen
		!compute R
		do i=1,min_MN
			do j=1,i-1
				Tdata(i,j)=0d0
			end do
		end do
		R=Tdata(1:min_MN,:)
		dimen=(/min_MN/)+(tempdim.sub.2)
		R%TenDim=dimen
		return
	end subroutine
!********************************************************************
!!!!! make LQ decompostion: A = LQ, Q is an orthomomal matrix, and L is a lower triangle
!  matrix. The size of matrix A is M x N. the size of L is M x min(M,N), Q is  min(M,N) x N
!	T=(res.i.1)*(res.i.2)
!********************************************************************
	type(DTensorlink) function DLQlink(T) result(res)
		type(DTensor),intent(in)::T
		real*8,allocatable::  work(:),tau(:),Tdata(:,:)
		Type(DTensor)::v,vv,identity
		type(DTensornode),pointer::respointer
		type(Dimension)::dimen
		integer :: i,j,M,N,INFO,min_MN
		M = T.dim.1
		N = T.dim.2
		if(DgetRank(T).ne.2) then
			write(*,*)"ERROR in DQRdecomposition"
			write(*,*)"input Tensor should be a matrix"
			stop
		endif
		Tdata=T
		min_MN=min(M,N)
		allocate(tau(min_MN))
		allocate(work(M))
		call DGELQ2( M, N, Tdata, M, tau, WORK, INFO )
		if(info.ne.0) then
				write (*,*) "Error in DQRdecomposition ,info=",info
				stop
		end if
		!compute Q
		v=DzeroTen((/N/))
		identity=Deye(N,N)
		allocate(respointer)
		do i=1,min_MN
			if(i.ne.1)then
				call Dmodify(v,i-1,0d0)
			end if
			call Dmodify(v,i,1d0)
			call Dmodify(v,(/i+1,N/),Tdata(i,i+1:N))
			vv=tau(i)*(v.xx.v)
			if(i.eq.1)then
				respointer%Ten=identity-vv
			else
				respointer%Ten=(identity-vv)*(respointer%Ten)
			end if
		end do
		respointer%Ten=DsubTen(respointer%Ten,(/-1,1,min_MN/))
		dimen=(/min_MN/)+(T%TenDim.sub.2)
		respointer%Ten%TenDim=dimen
		call 	Dpush_back(res,respointer)
		!compute L
		nullify(respointer)
		allocate(respointer)
		do i=1,M
			do j=i+1,min_MN
				Tdata(i,j)=0d0
			end do
		end do
		respointer%Ten=Tdata(:,1:min_MN)
		dimen=(T%TenDim.sub.1)+(/min_MN/)
		respointer%Ten%TenDim=dimen
		call Dpush_forward(res,respointer)
		return
	end function
!********************************************************************
!!!!! make LQ decompostion: A = LQ, Q is an orthomomal matrix, and L is a lower triangle
!  matrix. The size of matrix A is M x N. the size of L is M x min(M,N), Q is  min(M,N) x N
!	T=res(2)*res(1),the output is a array of DTensor
!********************************************************************
	function DLQdecompose(T) result(res)
		type(DTensor),allocatable::res(:)
		type(DTensor),intent(in)::T
		real*8,allocatable::  work(:),tau(:),Tdata(:,:)
		Type(DTensor)::v,vv,identity
		type(Dimension)::dimen
		integer :: i,j,M,N,INFO,min_MN
		M = T.dim.1
		N = T.dim.2
		if(DgetRank(T).ne.2) then
			write(*,*)"ERROR in DQRdecomposition"
			write(*,*)"input Tensor should be a matrix"
			stop
		endif
		Tdata=T
		min_MN=min(M,N)
		allocate(tau(min_MN))
		allocate(work(M))
		call DGELQ2( M, N, Tdata, M, tau, WORK, INFO )
		if(info.ne.0) then
				write (*,*) "Error in DQRdecomposition ,info=",info
				stop
		end if
		!compute Q
		v=DzeroTen((/N/))
		identity=Deye(N,N)
		allocate(res(2))
		do i=1,min_MN
			if(i.ne.1)then
				call Dmodify(v,i-1,0d0)
			end if
			call Dmodify(v,i,1d0)
			call Dmodify(v,(/i+1,N/),Tdata(i,i+1:N))
			vv=tau(i)*(v.xx.v)
			if(i.eq.1)then
				res(1)=identity-vv
			else
				res(1)=(identity-vv)*res(1)
			end if
		end do
		res(1)=DsubTen(res(1),(/-1,1,min_MN/))
		dimen=(/min_MN/)+(T%TenDim.sub.2)
		res(1)%TenDim=dimen
		!compute L
		do i=1,M
			do j=i+1,min_MN
				Tdata(i,j)=0d0
			end do
		end do
		res(2)=Tdata(:,1:min_MN)
		dimen=(T%TenDim.sub.1)+(/min_MN/)
		res(2)%TenDim=dimen
		return
	end function
	subroutine DLQdecomposition1(T,L,Q)
		type(DTensor),intent(in)::T
		type(DTensor),intent(inout)::L,Q
		real*8,allocatable::  work(:),tau(:),Tdata(:,:)
		Type(DTensor)::v,vv,identity
		type(Dimension)::dimen
		integer :: i,j,M,N,INFO,min_MN
		M = T.dim.1
		N = T.dim.2
		if(DgetRank(T).ne.2) then
			write(*,*)"ERROR in DQRdecomposition"
			write(*,*)"input Tensor should be a matrix"
			stop
		endif
		Tdata=T
		min_MN=min(M,N)
		allocate(tau(min_MN))
		allocate(work(M))
		call DGELQ2( M, N, Tdata, M, tau, WORK, INFO )
		if(info.ne.0) then
				write (*,*) "Error in DQRdecomposition ,info=",info
				stop
		end if
		!compute Q
		v=DzeroTen((/N/))
		identity=Deye(N,N)
		do i=1,min_MN
			if(i.ne.1)then
				call Dmodify(v,i-1,0d0)
			end if
			call Dmodify(v,i,1d0)
			call Dmodify(v,(/i+1,N/),Tdata(i,i+1:N))
			vv=tau(i)*(v.xx.v)
			if(i.eq.1)then
				Q=identity-vv
			else
				Q=(identity-vv)*Q
			end if
		end do
		Q=DsubTen(Q,(/-1,1,min_MN/))
		dimen=(/min_MN/)+(T%TenDim.sub.2)
		Q%TenDim=dimen
		!compute L
		do i=1,M
			do j=i+1,min_MN
				Tdata(i,j)=0d0
			end do
		end do
		L=Tdata(:,1:min_MN)
		dimen=(T%TenDim.sub.1)+(/min_MN/)
		L%TenDim=dimen
		return
	end subroutine
! Q store in T
	subroutine DLQdecomposition2(T,L)
		type(DTensor),intent(inout)::L,T
		real*8,allocatable::  work(:),tau(:),Tdata(:,:)
		Type(DTensor)::v,vv,identity
		type(Dimension)::dimen,tempdim
		integer :: i,j,M,N,INFO,min_MN
		M = T.dim.1
		N = T.dim.2
		if(DgetRank(T).ne.2) then
			write(*,*)"ERROR in DQRdecomposition"
			write(*,*)"input Tensor should be a matrix"
			stop
		endif
		Tdata=T
		tempdim=T%TenDim
		min_MN=min(M,N)
		allocate(tau(min_MN))
		allocate(work(M))
		call DGELQ2( M, N, Tdata, M, tau, WORK, INFO )
		if(info.ne.0) then
				write (*,*) "Error in DQRdecomposition ,info=",info
				stop
		end if
		call cleanDTensor(T)
		!compute Q
		v=DzeroTen((/N/))
		identity=Deye(N,N)
		do i=1,min_MN
			if(i.ne.1)then
				call Dmodify(v,i-1,0d0)
			end if
			call Dmodify(v,i,1d0)
			call Dmodify(v,(/i+1,N/),Tdata(i,i+1:N))
			vv=tau(i)*(v.xx.v)
			if(i.eq.1)then
				T=identity-vv
			else
				T=(identity-vv)*T
			end if
		end do
		T=DsubTen(T,(/-1,1,min_MN/))
		dimen=(/min_MN/)+(tempdim.sub.2)
		T%TenDim=dimen
		!compute L
		do i=1,M
			do j=i+1,min_MN
				Tdata(i,j)=0d0
			end do
		end do
		L=Tdata(:,1:min_MN)
		dimen=(tempdim.sub.1)+(/min_MN/)
		L%TenDim=dimen
		return
	end subroutine
!ccccccccccccccccc svd cccccccccccccccccc	
!ccccccccccccccccc svd cccccccccccccccccc	
! 			T should be two dimension		
!			T=U*s*V^T
! on output,the order in the link is
! [U]->[s]->[V^T]
! T=(link.i.1)*(eye(link.i.2))*(link.i.3)
	type(DTensorlink) function Dsvd_no_cut(T) result(res)
		type(DTensor),intent(in)::T
		type(DTensornode),pointer::Tponter
		real*8,allocatable :: T1data(:)
		real*8,allocatable :: T2data(:)
		real*8,allocatable :: Tdata(:)
		real*8,allocatable :: work(:)
		integer m,n,info,lw,i,j
		integer,allocatable :: iw(:)
		real*8,allocatable ::sdata(:)
		integer :: ms_max,min_MN,max_MN
		type(Dimension) :: T1Dim,T2Dim
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
		m=T.dim.1
		n=T.dim.2
		min_MN=min(M,N)
		max_MN=max(M,N)
		lw=3*min_MN*min_MN+max(max_MN,4*min_MN*min_MN+4*min_MN)
		allocate(iw(8*min_MN))
		allocate(T1data(m*min_MN))
		allocate(T2data(min_MN*n))
		allocate(sdata(min_MN)) 
		allocate(work(lw))
		Tdata=T
		call DGESDD('S',m,n,TData,m,sdata,T1data,m,T2data,min_MN,WORK,lw,iw,INFO)
		if(info.ne.0) then
			write (*,*) "Error in svd ,info=",info
		end if
		T1Dim=(/min_MN/)
		T1Dim=(T%TenDim.sub.1)+T1Dim
		T2Dim=(/min_MN/)
		T2Dim=T2Dim+(T%TenDim.sub.2)
		allocate(Tponter)
		Tponter%Ten=DbuildTen(T1Dim,T1data)
		call Dpush_back(res,Tponter)
		nullify(Tponter)
		
		allocate(Tponter)
		Tponter%Ten=DbuildTen(sdata)
		call Dpush_back(res,Tponter)
		nullify(Tponter)
		
		
		allocate(Tponter)
		Tponter%Ten=DbuildTen(T2Dim,T2data)
		call Dpush_back(res,Tponter)
		nullify(Tponter)
		return
	end function 
	type(DTensorlink) function Dsvd_cut(T,Ncut_) result(res)
		type(DTensor),intent(in)::T
		integer,intent(in)::Ncut_
		integer::Ncut
		type(DTensornode),pointer::Tponter
		real*8,allocatable :: T1data(:,:)
		real*8,allocatable :: T2data(:,:)
		real*8,allocatable :: Tdata(:)
		real*8,allocatable :: work(:)
		integer m,n,info,lw,i,j
		integer,allocatable :: iw(:)
		real*8,allocatable ::sdata(:)
		integer :: ms_max,min_MN,max_MN
		type(Dimension) :: T1Dim,T2Dim
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
		m=T.dim.1
		n=T.dim.2
		min_MN=min(M,N)
		max_MN=max(M,N)
		lw=3*min_MN*min_MN+max(max_MN,4*min_MN*min_MN+4*min_MN)
		allocate(iw(8*min_MN))
		allocate(T1data(m,min_MN))
		allocate(T2data(min_MN,n))
		allocate(sdata(min_MN)) 
		allocate(work(lw))
		if(Ncut_.gt.min_MN) then
			Ncut=min_MN
		else
			Ncut=Ncut_
		end if
		Tdata=T
		call DGESDD('S',m,n,TData,m,sdata,T1data,m,T2data,min_MN,WORK,lw,iw,INFO)
		if(info.ne.0) then
			write (*,*) "Error in svd ,info=",info
		end if
		
		T1Dim=(/Ncut/)
		T1Dim=(T%TenDim.sub.1)+T1Dim
		T2Dim=(/Ncut/)
		T2Dim=T2Dim+(T%TenDim.sub.2)
		
		allocate(Tponter)
		Tponter%Ten=T1data(:,1:Ncut)
		call Dresetdim(Tponter%Ten,T1Dim)
		call Dpush_back(res,Tponter)
		nullify(Tponter)
		
		allocate(Tponter)
		Tponter%Ten=DbuildTen(sdata(1:Ncut))
		call Dpush_back(res,Tponter)
		nullify(Tponter)
		
		
		allocate(Tponter)
		Tponter%Ten=T2data(1:Ncut,:)
		call Dresetdim(Tponter%Ten,T2Dim)
		call Dpush_back(res,Tponter)
		nullify(Tponter)
		return
	end function 
	function Dsvddecompose(T,Ncut_) result(res)
		type(DTensor),allocatable::res(:)
		type(DTensor),intent(in)::T
		integer,optional,intent(in)::Ncut_
		integer::Ncut
		real*8,allocatable :: T1data(:,:)
		real*8,allocatable :: T2data(:,:)
		real*8,allocatable :: Tdata(:)
		real*8,allocatable :: work(:)
		integer m,n,info,lw,i,j
		integer,allocatable :: iw(:)
		real*8,allocatable ::sdata(:)
		integer :: ms_max,min_MN,max_MN
		type(Dimension) :: T1Dim,T2Dim
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
		m=T.dim.1
		n=T.dim.2
		min_MN=min(M,N)
		max_MN=max(M,N)
		lw=3*min_MN*min_MN+max(max_MN,4*min_MN*min_MN+4*min_MN)
		allocate(iw(8*min_MN))
		allocate(T1data(m,min_MN))
		allocate(T2data(min_MN,n))
		allocate(sdata(min_MN)) 
		allocate(work(lw))
		Tdata=T
		call DGESDD('S',m,n,TData,m,sdata,T1data,m,T2data,min_MN,WORK,lw,iw,INFO)
		if(info.ne.0) then
			write (*,*) "Error in svd ,info=",info
		end if
		allocate(res(3))
		if(present(Ncut_))then
			if(Ncut_.gt.min_MN) then
				Ncut=min_MN
			else
				Ncut=Ncut_
			end if
			T1Dim=(/Ncut/)
			T1Dim=(T%TenDim.sub.1)+T1Dim
			T2Dim=(/Ncut/)
			T2Dim=T2Dim+(T%TenDim.sub.2)
		
			res(1)=T1data(:,1:Ncut)
			call Dresetdim(res(1),T1Dim)
			res(2)=DbuildTen(sdata(1:Ncut))
			res(3)=T2data(1:Ncut,:)
			call Dresetdim(res(3),T2Dim)
			return
		else
			T1Dim=(/min_MN/)
			T1Dim=(T%TenDim.sub.1)+T1Dim
			T2Dim=(/min_MN/)
			T2Dim=T2Dim+(T%TenDim.sub.2)
			res(1)=DbuildTen(T1Dim,T1data)
			res(2)=DbuildTen(sdata)
			res(3)=DbuildTen(T2Dim,T2data)
			return
		end if
	end function 
!******************************************************
! 			T should be two dimension		
!			T=U*s*V^T
!			 on output T=U*Deye(s)*V
	subroutine Dsvddecomposition1(T,U,s,V)
		type(DTensor),intent(in)::T
		type(DTensor),intent(out)::U,s,V
		real*8,allocatable :: T1data(:)
		real*8,allocatable :: T2data(:)
		real*8,allocatable :: Tdata(:)
		real*8,allocatable :: work(:)
		integer m,n,info,lw,i,j
		integer,allocatable :: iw(:)
		real*8,allocatable ::sdata(:)
		integer :: ms_max,min_MN,max_MN
		type(Dimension) :: T1Dim,T2Dim
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
		m=T.dim.1
		n=T.dim.2
		min_MN=min(M,N)
		max_MN=max(M,N)
		lw=3*min_MN*min_MN+max(max_MN,4*min_MN*min_MN+4*min_MN)
		allocate(iw(8*min_MN))
		allocate(T1data(m*min_MN))
		allocate(T2data(min_MN*n))
		allocate(sdata(min_MN)) 
		allocate(work(lw))
		Tdata=T
		call DGESDD('S',m,n,TData,m,sdata,T1data,m,T2data,min_MN,WORK,lw,iw,INFO)
		if(info.ne.0) then
			write (*,*) "Error in svd ,info=",info
		end if
		T1Dim=(/min_MN/)
		T1Dim=(T%TenDim.sub.1)+T1Dim
		T2Dim=(/min_MN/)
		T2Dim=T2Dim+(T%TenDim.sub.2)
		U=DbuildTen(T1Dim,T1data)
		s=DbuildTen(sdata)
		V=DbuildTen(T2Dim,T2data)
		return
	end subroutine 
	subroutine DSVDcutoff(T,U,s,V,Ncut_)
		type(DTensor),intent(in)::T
		type(DTensor),intent(out)::U,s,V
		integer,intent(in)::Ncut_
		integer::Ncut
		real*8,allocatable :: T1data(:,:)
		real*8,allocatable :: T2data(:,:)
		real*8,allocatable :: Tdata(:)
		real*8,allocatable :: work(:)
		integer m,n,info,lw,i,j
		integer,allocatable :: iw(:)
		real*8,allocatable ::sdata(:)
		integer :: ms_max,min_MN,max_MN
		type(Dimension) :: T1Dim,T2Dim
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
		m=T.dim.1
		n=T.dim.2
		min_MN=min(M,N)
		max_MN=max(M,N)
		lw=3*min_MN*min_MN+max(max_MN,4*min_MN*min_MN+4*min_MN)
		allocate(iw(8*min_MN))
		allocate(T1data(m,min_MN))
		allocate(T2data(min_MN,n))
		allocate(sdata(min_MN)) 
		allocate(work(lw))
		if(Ncut_.gt.min_MN) then
			Ncut=min_MN
		else
			Ncut=Ncut_
		end if
		Tdata=T
		call DGESDD('S',m,n,TData,m,sdata,T1data,m,T2data,min_MN,WORK,lw,iw,INFO)
		if(info.ne.0) then
			write (*,*) "Error in svd ,info=",info
		end if
		
		T1Dim=(/Ncut/)
		T1Dim=(T%TenDim.sub.1)+T1Dim
		T2Dim=(/Ncut/)
		T2Dim=T2Dim+(T%TenDim.sub.2)
		
		U=T1data(:,1:Ncut)
		call Dresetdim(U,T1Dim)
		
		S=DbuildTen(sdata(1:Ncut))
		
		
		V=T2data(1:Ncut,:)
		call Dresetdim(V,T2Dim)
		return
	end subroutine 
!!! The inverse of a matrix: the input tensor should be a square matrix 
	type(Dtensor) function Dinverse(T)
		type(DTensor),intent(in) :: T
		type(DTensor):: E
		integer :: M,N
		if(DgetRank(T).ne.2) then
			write(*,*)"ERROR in calculating the inverse of a Dtensor"
			write(*,*)"input Tensor should be a square matrix"
			stop
		endif
		M = T.dim.1
		N = T.dim.2
		if(M.ne.N) then
			write(*,*)"ERROR in calculating the inverse of a Dtensor"
			write(*,*)"input Tensor should be a square matrix"
			stop
		endif
		
		E=Deye(M,N)
		Dinverse=Dlinequ(T,E)
		return
	end function
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
	type(DTensor) function deye_real2(s)
		real*8,intent(in) :: s(:)
		integer::n
		real*8,allocatable :: sdata(:,:)
		integer :: i
		n=size(s)
		allocate(sdata(n,n))
		sdata=0d0
		do i=1,n
			sdata(i,i)=s(i)
		end do
		call dstoreTenData(deye_real2,sdata)
		deye_real2%rank=2
		deye_real2%totalData=n*n
		deye_real2%TenDim=(/n,n/)
		deye_real2%flag=.true.
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
!		inde=[-1,inde_min,inde_max] output data(inde_min:inde_max,:)
!	or[-2,inde_min,inde_max],data(:,inde_min:inde_max)
!	or [-1,inde_row] [-2,inde_col],output row or col
!or [inde1_min,inde1_max,inde2_min,inde2_max] output data(inde1_min:inde1_max,inde2_min:inde2_max)

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
		if(size(inde).eq.4) then
			DsubTen=Tdata(inde(1):inde(2),inde(3):inde(4))
			return
		end if
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
	type(DTensor) function DdirectProduct_Matlab_korn(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: T1data(:,:),T2data(:,:)
		real*8,allocatable :: newdata(:,:,:,:)
		integer :: m1,n1,m2,n2,i,j,k,l,rank(2)
		type(Dimension):: D1,Dtemp
		if((DgetRank(T1).eq.1).and.(DgetRank(T2).eq.1)) then
			DdirectProduct_Matlab_korn=DdirectProduct1_Matlab_korn(T1,T2)
			return
		end if
		if((DgetRank(T1).eq.2).and.(DgetRank(T2).eq.2)) then
			D1=T2%TenDim.sub.1
			Dtemp=T1%TenDim.sub.1
			D1=D1+Dtemp
			Dtemp=T2%TenDim.sub.2
			D1=D1+Dtemp
			Dtemp=T1%TenDim.sub.2
			D1=D1+Dtemp
			m1=T1.dim.1
			n1=T1.dim.2
			m2=T2.dim.1
			n2=T2.dim.2
			allocate(T1data(m1,n1))
			allocate(T2data(m2,n2))
			allocate(newdata(m2,m1,n2,n1))
			T1data=reshape(T1%DTensor_Data,(/m1,n1/))
			T2data=reshape(T2%DTensor_Data,(/m2,n2/))
			do l=1,n2
				do k=1,n1
					 do j=1,m2
						 do i=1,m1
							newdata(j,i,l,k)=T1data(i,k)*T2data(j,l)
						end do
					end do
				end do
			end do
			DdirectProduct_Matlab_korn=Dset(D1,reshape(newdata,(/m1*m2*n1*n2/)))
		else
			write(*,*)"The DdirectProduct_Matlab_korn is just for matrix"
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
	type(DTensor) function DdirectProductM_Matlab_korn(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: T1data(:,:),T2data(:,:)
		real*8,allocatable :: newdata(:,:,:,:)
		integer :: m1,n1,m2,n2,i,j,k,l
		type(Dimension):: D1,Dtemp
		if((DgetRank(T1).eq.1).and.(DgetRank(T2).eq.1)) then
			DdirectProductM_Matlab_korn=DdirectProduct1V_Matlab_korn(T1,T2)
			return
		end if
		if((DgetRank(T1).eq.2).and.(DgetRank(T2).eq.2)) then
			D1=T2%TenDim.sub.1
			Dtemp=T1%TenDim.sub.1
			D1=D1+Dtemp
			Dtemp=T2%TenDim.sub.2
			D1=D1+Dtemp
			Dtemp=T1%TenDim.sub.2
			D1=D1+Dtemp
			D1=DimConstract(D1,1,1)
			D1=DimConstract(D1,2,1)
			m1=T1.dim.1
			n1=T1.dim.2
			m2=T2.dim.1
			n2=T2.dim.2
			allocate(T1data(m1,n1))
			allocate(T2data(m2,n2))
			allocate(newdata(m2,m1,n2,n1))
			T1data=reshape(T1%DTensor_Data,(/m1,n1/))
			T2data=reshape(T2%DTensor_Data,(/m2,n2/))
			do l=1,n2
				do k=1,n1
					do j=1,m2
						do i=1,m1
							newdata(j,i,l,k)=T1data(i,k)*T2data(j,l)
						end do
					end do
				end do
			end do
			DdirectProductM_Matlab_korn=Dset(D1,reshape(newdata,(/m1*m2*n1*n2/)))
		else
			write(*,*)"The DdirectProductM_Matlab_korn is just for matrix"
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
		type(DTensor) function DdirectProduct1V_Matlab_korn(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: newdata(:,:)
		integer :: m1,n1,m2,n2,i,j,k,l
		type(Dimension):: D1,Dtemp
		if((DgetRank(T1).eq.1).and.(DgetRank(T2).eq.1)) then
			D1=T2%TenDim
			D1=D1+T1%TenDim
			D1=DimConstract(D1,1,1)
			n1=T1.dim.1
			n2=T2.dim.1
			allocate(newdata(n2,n1))
			do l=1,n2
				do k=1,n1
					newdata(l,k)=T1%DTensor_Data(k)*T2%DTensor_Data(l)
				end do
			end do
			DdirectProduct1V_Matlab_korn=Dset(D1,reshape(newdata,(/n1*n2/)))
		else
			write(*,*)"The DdirectProduct1_Matlab_korn is just for vector"
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
	
	type(DTensor) function DdirectProduct1_Matlab_korn(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: newdata(:,:)
		integer :: m1,n1,m2,n2,i,j,k,l
		type(Dimension):: D1,Dtemp
		if((DgetRank(T1).eq.1).and.(DgetRank(T2).eq.1)) then
			D1=T2%TenDim
			D1=D1+T1%TenDim
			n1=T1.dim.1
			n2=T2.dim.1
			allocate(newdata(n2,n1))
			do l=1,n2
				do k=1,n1
					newdata(l,k)=T1%DTensor_Data(k)*T2%DTensor_Data(l)
				end do
			end do
			DdirectProduct1_Matlab_korn=Dset(D1,reshape(newdata,(/n1*n2/)))
		else
			write(*,*)"The DdirectProduct1_Matlab_korn is just for vector"
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
!**************   Dcombination   ********************
!
!                    / T1 \
!	Dcombination(T1,T2)	=  |----|
!		               \ T2 /
!	T1 :a [...,m,n,l] matrix
!	T2 :a [...,m,n,l] matrix
!	Dcombination(T1,T2):a [...,m,n,l,2] matrix
!or 
!	T1 :a [...,m,n,l] matrix
!	T2 :a [...,m,n] matrix
!	Dcombination(T1,T2):a [...,m,n,l+1] matrix
	type(DTensor) function Dcombination(T1,T2)
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
			call DstoreTenData(Dcombination,newdata)
			Dcombination%Rank=DgetRank(T1)+1
			newDim=.subDim.T1
			newDim=newDim+(/2/)
			Dcombination%TenDim=newDim
			Dcombination%totalData=total
			Dcombination%flag=.true.
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
		call DstoreTenData(Dcombination,newdata)
		Dcombination%Rank=DgetRank(T1)
		newDim=(.subDim.T2)
		newDim=newDim+(/dim_n+1/)
		Dcombination%TenDim=newDim
		Dcombination%totalData=total+total1
		Dcombination%flag=.true.
	end function
	type(DTensor) function Dcombinationrow(T1,T2)
		type(DTensor),intent(in)::T1,T2
		integer,allocatable::dim1(:),dim2(:)
		integer::total,total1,i,dim_n
		real*8,allocatable::newdata(:,:),olddata(:,:)
		type(Dimension)::newDim,dimen1,dimen2
		if(T1%rank.eq.T2%rank) then
			dim1=.subDim.T1
			dim2=.subDim.T2
			if(.not.(dim1.equ.dim2)) then
				write(*,*)"can not Dcombinationrow two DTensor"
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
			call DstoreTenData(Dcombinationrow,newdata)
			Dcombinationrow%Rank=DgetRank(T1)+1
			newDim=.subDim.T1
			newDim=(/2/)+newDim
			Dcombinationrow%TenDim=newDim
			Dcombinationrow%totalData=total
			Dcombinationrow%flag=.true.
			return
		end if
		dimen1=.subDim.T1
		dimen1=DimConstract(dimen1,2,T1%rank)
		dimen2=.subDim.T2
		dimen2=DimConstract(dimen2,1,T2%rank)
		if(.not.((dimen1.sub.2).equ.dimen2)) then
			write(*,*)"can not Dcombinationrow two DTensor"
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
		call DstoreTenData(Dcombinationrow,newdata)
		Dcombinationrow%Rank=DgetRank(T1)
		newDim=(.subDim.T2)
		newDim=(/dim_n+1/)+newDim
		Dcombinationrow%TenDim=newDim
		Dcombinationrow%totalData=total+total1
		Dcombinationrow%flag=.true.
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
! modify some the data in the Tensor,the Tensor should be rank=1
!	Ten(inde(1):inde(2))=val	
	subroutine DmodifyTen_some_data1_real(Ten,inde,val)
		type(DTensor),intent(inout)::Ten
		real*8,intent(in)::val(:)
		integer,intent(in)::inde(:)
		if(Ten%rank.ne.1) then
			write(*,*)"ERROR in Dmodify Tensor data"
			write(*,*)"only rank=1 is allowed,stop"
			stop
		end if
		if((inde(2)-inde(1)+1).ne.size(val)) then
			write(*,*)"ERROR in modify Tensor data"
			write(*,*)"ERROR DmodifyTen_some_data1"
			stop
		end if
		if(inde(2).gt.Ten%totaldata) then
			write(*,*)"ERROR in modify Tensor data"
			write(*,*)"ERROR DmodifyTen_some_data1"
			write(*,*)"ERROR 2"
			stop
		end if
		Ten%DTensor_Data(inde(1):inde(2))=val
		return
	end subroutine	
! Dmodify all the dimension in the DTensor
	subroutine Dresetdim1(Ten,dimen)
		type(DTensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		!type(Dimension)::dimensi
		if(product(dimen).ne.size(Ten%DTensor_Data)) then
			write(*,*)"ERROR in resetdim"
			write(*,*)product(dimen),size(Ten%DTensor_Data)
			write(*,*)"stop"
			stop
		end if
		!dimensi=dimen
		Ten%TenDim=dimen
		Ten%rank=size(dimen)
		return
	end subroutine	
	subroutine Dresetdim2(Ten,dimen)
		type(DTensor),intent(inout)::Ten
		type(dimension),intent(in)::dimen
		!integer,allocatable::dimensi(:)
		!dimensi=dimen
		if(outtotaldata(dimen).ne.size(Ten%DTensor_Data)) then
			write(*,*)"ERROR in Dresetdim"
			write(*,*)outtotaldata(dimen),size(Ten%DTensor_Data)
			write(*,*)"stop"
			stop
		end if
		Ten%TenDim=dimen
		Ten%rank=DimSize(Dimen)
		return
	end subroutine	
!*****************************
!		Bubble Sort
!	only rank=1 is allowed
!			
	type(DTensor) function DSortTen(T,increase)
		type(DTensor),intent(in)::T
		logical,intent(in)::increase
		if(T%rank.ne.1) then
			write(*,*)"Sort,only rank=1 is allowed"
			write(*,*)"stop"
			stop
		end if
		DSortTen%DTensor_Data=T
		call sort(DSortTen%DTensor_Data,increase)
		DSortTen%TenDim=T%TenDim
		DSortTen%rank=T%rank
		DSortTen%flag=T%flag
		DSortTen%totalData=T%totalData
		return
	end function		
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
			write(*,*)"length of index of input"
			write(*,*)size(inde)
			write(*,*)"length of the index in the link is"
			write(*,*)DlenOfIndice(h)
			stop
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
			write(*,*)"length of index of input"
			write(*,*)size(inde)
			write(*,*)"length of the index in the link is"
			write(*,*)DlenOfIndice(h)
			stop
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
		!call Dcleanlink(hout)
		hout%head=>hin%head
		hout%Tend=>hin%Tend
		hout%length=hin%length
		return
	end subroutine
	subroutine Dnullifylink(h)
		type(DTensorlink),intent(inout) :: h
		nullify(h%head)
		nullify(h%Tend)
		h%length=0
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
			!call Dcleanlink(link)
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
!		call Dcleanlink(link)
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
				if(associated(Tp%next)) then
					Tp=>Tp%next
				else
					write(*,*)"ERROR in DLiv",j
					stop
				endif
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
	subroutine Dnullifylist(h)
		type(DTensorlist),intent(inout) :: h
		nullify(h%head)
		nullify(h%Lend)
		h%length=0
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
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!**************  there is no real*8 type function ***************************************
! This part written by Wenyuan Liu
!
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************

!!!!!  normalize: the tensor elements are divided by the maximal abs(elements)
	type(DTensor) function Dnormalize(temp) result(T) 
		type(DTensor),intent(in) :: temp
		real*8	:: maxi
		maxi=maxval(abs(temp%DTensor_data))
		T = temp/maxi
		return
	end function
!!!!  sqrttensor : output a tensor whose element value is corresponding  square root of the input tensor element value
	type(DTensor) function sqrttensor(T) 
		type(DTensor),intent(in) :: T
		real*8,allocatable :: temparray(:)
		integer,allocatable :: Tdim(:)
		integer :: i
		Tdim=T%TenDim
		allocate(temparray(T%totaldata))
		do i=1,T%totaldata
			temparray(i)= sqrt(T%DTensor_data(i))
		end do
		sqrttensor = DbuildTen(Tdim,tempArray)
		return
	end function
!!  give a part of a square matrix to another square matrix:
!! eg: input matrix T1 is size of 10 x 10 , we can cut it, and let T2 = T1(1:5,1:5). 
	type (DTensor) function Tensorcut(T,M)
		implicit none 
		type(DTensor),intent(in)::T
		integer,intent(in) :: M
		real*8,allocatable :: temp(:,:),cut(:,:)
		integer :: i,dim1,dim2
		dim1 = T.dim.1
		dim2 = T.dim.2
		if(T%rank/=2 .or. dim1/=dim2 ) then 
			write(*,*) 'ERROR in Tensorcut'
			stop
		endif
		allocate(temp(dim1,dim2))
		allocate(cut(M,M))
		temp=T
		cut=0d0
		do i=1,M
			cut(i,i)=temp(i,i)
		enddo
		
		Tensorcut=cut

		return
	end function

!BuildTen: this function is used in DQRdecomposition, DLQdecomposition
!and DSVDcutoff
	type(DTensor) function DbuildTen8_2(TenDim,DTensor_data)
		type(Dimension),intent(in)::TenDim
		integer,allocatable::Dimen(:)
		real*8,intent(in)::DTensor_data(:,:)
		Dimen=TenDim
		DbuildTen8_2%rank=size(Dimen)
		DbuildTen8_2%TenDim=TenDim

		if(product(Dimen).ne.product(shape(DTensor_Data))) then
			write(*,*)"ERROR in dimension and input data in function DbuildTen"
			write(*,*)"stop"
			stop
		end if

		call DstoreTenData(DbuildTen8_2,DTensor_Data)
		DbuildTen8_2%totalData=size(DTensor_Data)
		DbuildTen8_2%flag=.true.
		
		return
	end function	
	
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************










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
		
		if((proID.ne.ID1).and.(proID.ne.ID2)) return!The proID do not sent or recv, return
		
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
		!	if(allocated(Ten2%DTensor_data)) then
		!		deallocate(Ten2%DTensor_data)
		!	end if
			if(allocated(Ten2%DTensor_data))then
				if(Ten2%totalData.ne.size(Ten2%DTensor_data)) then
					deallocate(Ten2%DTensor_data)
					allocate(Ten2%DTensor_data(Ten2%totalData))
				end if
			else
				allocate(Ten2%DTensor_data(Ten2%totalData))
			end if
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
		!	if(allocated(Ten1%DTensor_data)) then
		!		deallocate(Ten1%DTensor_data)
		!	end if
			if(allocated(Ten1%DTensor_data)) then
				if(Ten1%totalData.ne.size(Ten1%DTensor_data))then
					deallocate(Ten1%DTensor_data)
					allocate(Ten1%DTensor_data(Ten1%totalData))
				end if
			else
				allocate(Ten1%DTensor_data(Ten1%totalData))
			end if
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
		
		if((proID.ne.ID1).and.(proID.ne.ID2)) return!The proID do not sent or recv, return
		
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
					!if(allocated(indice)) then
					!	deallocate(indice)
					!end if
					if(allocated(indice))then
						if(indeLen.ne.size(indice)) then
							deallocate(indice)
							allocate(indice(indeLen))
						end if
					else
						allocate(indice(indeLen))
					end if
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
				!	if(allocated(indice)) then
				!		deallocate(indice)
				!	end if
					if(allocated(indice))then
						if(indeLen.ne.size(indice))then
							deallocate(indice)
							allocate(indice(indeLen))
						end if
					else
						allocate(indice(indeLen))
					end if
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


			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
