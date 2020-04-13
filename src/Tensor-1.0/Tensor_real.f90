
!		 This file is modify from Tensor.f90 which is used for complex*16 type
!	no test on this code
!		The file Tensor.f90 modify at 2013.10.29.Add some function,but do not
!	add in this file.	

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
!		allocaDTing storage space for A,A will get no data.This happen in 
!		fortran90 but do not in fortran77.The interface assignments is to  
!		solve this program.
!			3.The code is fortran90 version.
!*****************     note      *******************
!			1.The product of DTensor only allow rank<=2,if larger than 2,
!		use .p. and Dcontract(or .c.) to reshape it.
!			2.permutaDTion(.p.) is only allow for rank<=3,if larger then 3,use 
!		Dcontract(or .c.) to reshape it to rank=3 or rank=2,and then do the 
!		permutaDTion,and at last Ddecompose the DTensor back to the original rank.
!			3.to compile the code one should link the files to lapack and blas.
!			4.Send email to tyrants@qq.com to report any bugs.
!***********************************************************

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
		integer,allocatable,private:: indice(:)
		type(DTensor) ::Ten
		type(DTensornode),private,pointer :: next
	end type DTensornode
		
	type DTensorlink
		integer,private :: length=0
		type(DTensornode),pointer,private  :: head
		type(DTensornode),pointer,private  :: Tend
	end type DTensorlink	
!*********		End define		*********************


!****************************************************
!*********		define private data	***************	
	real*8,private :: cpuDTime(2)
	real*8,parameter,private :: II=(0,1),EE=(1,0)! II**2=-1
	real*8,parameter,private :: large_number_warning=1d80!use for checking
	logical,private,save::check_flag=.false.!use for checking
	real*8,private::DDOT!dot product
	real*8,private::DNRM2!norm
	External::DDOT,DNRM2!blas's function
	private:: Dwarning,Dwarning2,DNANwarning,DNANwarning2!They are use for check the input data
	private:: DNANjudge!They are use for check the input data
	private::DstoreTenData
	private::DaddressToIndes
	private::DaddressToIndes2
	private::DIndesToaddress
!*********		End define	***************	
	
!***********   interface for Tensorlink   *************************	
	interface Dpush_back
		module procedure Dpush_backTen
		module procedure Dpush_backnode
		module procedure Dpush_backI
		module procedure Dpush_backnum
		module procedure Dpush_backnumI
	end interface
	
	interface Dpush_forward
		module procedure Dpush_fo
		module procedure Dpush_fonode
		module procedure Dpush_foI
		module procedure Dpush_fonum
		module procedure Dpush_fonumI
	end interface
	
	interface Dnode_i
		module procedure Dnode_i_int
		module procedure Dnode_i_vec
	end interface		
	
	interface Dset!create a DTensor,giving dimension and DTensor_Data
		module procedure Dset1 
		module procedure Dset2 
	end interface
	interface	sqrt
		module procedure DsqrtTen
	end interface
!	DTcopy the array to the T%DTensor_Data,array can be 1,2 or 3 dimension
	interface DstoreTenData
		module procedure DstoreTenData_dim1
		module procedure DstoreTenData_dim2
		module procedure DstoreTenData_dim3
	end interface
	
	interface	DbuildTen!build a DTensor
		module procedure DbuildTen1
		module procedure DbuildTen2
		module procedure DbuildTen3
		module procedure DbuildTen4
		module procedure DbuildTen5
		module procedure DbuildTen6
		module procedure DbuildTen8
	end interface
	
	interface Dmodify
		module procedure Dmodify_link
		module procedure DmodifyTen_val
		module procedure DmodifyTen_dat
		module procedure DmodifyTen_dat2
		module procedure DmodifyTen_dat3
	end interface
	
	interface Deye
!  return a matrix of diag[s1,s2,s3,...s_{max},0,0...]	,a m*n matrix
		module procedure Deye_real!input s(:) is real*8
		module procedure Done_com!idenDTity matrix
		module procedure Done_Ten!return a idenDTity DTensor 
	end interface

	interface operator(.p.)
		module procedure Dpermute_rank2 !permute the DTensor whose rank is 2
		module procedure Dpermute_rank3 !permute the DTensor whose rank is 3
	end interface
	interface operator(.pf.)
		module procedure Dpermutefo!permute the inde index to the first
										  !T_{1,2,3,..,i,..,n},permutefo(T,i)=_{i,1,2,3,..,i-1,i+1,..,n}
	end interface
	interface operator(.pi.)!T_{1,2,3,..,i,..,n},permuteInde(T,i)=_{2,3,..,i,1,i+1,..,n}
		module procedure DpermuteInde
	end interface
	
		
	interface operator(+)
		module procedure Dadd!the DElement of one DTensor push the other one
		module procedure Dadd_num!one DTensor push an idenDTity DTensor 
		module procedure Dadd_num_!num+ DTensor
		module procedure Dconnectlink!combine two DTensor link
	end interface
	
	
	interface operator(.c.)!combine two index of the DTensor,which is con_index and con_index+1
		module procedure Dcontract
		module procedure Dcontractsv
	end interface
	interface operator(.cd.)
		module procedure Dcompose_decompse
	end interface
	interface operator(.dc.)
		module procedure Ddecompose1
		module procedure DdecomposeAll
		module procedure Ddecomposev
	end interface
	
	interface operator(-)
		module procedure Dminus!the DElement of one DTensor minue the other one
		module procedure Dminus_num!one DTensor minue an idenDTity DTensor 
		module procedure Dminus_num_!num - DTensor
	end interface
	
	interface operator(*)
		module procedure DproductTen!Matrixe product,two DTensor intput should be of 1 or 2 dimension
		module procedure DmulDTiply_number!A DTensor DTimes a number,T*num
		module procedure DmulDTiply_number_!num*T
		module procedure DmulDTiply_real4!A DTensor DTimes a real*4 number
	end interface
	
	interface operator(/)
		module procedure Ddivide!the DElement of one DTensor Ddivide the DElement of other one
		module procedure Ddivide_number!A DTensor Ddivide a number of type real*8
	end interface
	
	
	interface operator(.xx.)
		module procedure DdirectProduct! Ddirect Product [m1,n1] * [m2,n2] => [m1,m2,n1,n2]
	end interface
		
	interface operator(.mxx.)
		module procedure DdirectProductM! Ddirect Product return a matrix [m1,n1] * [m2,n2] => [(m1*m2),(n1*n2)]
	end interface
		
	interface operator(.su.)
		module procedure Ddirect_sum! Ddirect sum
	end interface
		
	interface operator(.msu.)
		module procedure Ddirect_sumM! Ddirect sum return a matrix
	end interface
		
	interface operator(.elex.)
		module procedure DmulDTiply_DTensor!the DElement of one DTensor DTimes the DElement of other one
	end interface
	
	interface operator(.numx.)
		module procedure DTensorProduct1Dim!two one-dimension DTensors product,return a number of type real*8
	end interface
	
	interface operator(.H.)
		module procedure Dtranspose! DconjugateTranspose! one or two dimension case
	end interface
	
	interface operator(.equ.)
		module procedure equal_of_DTensor!If two DTensor are equal
	end interface
	
	interface operator(.subDim.)
		module procedure DgetTenSubDim
		module procedure DgetTenDim!output the Dimension,return type(Dimension)
	end interface

	interface operator(.i.)
		module procedure DTi!output the i DTensor in the DTensorlink
		module procedure DElement!output the data in the DTensor of DTim=(i,j)
		module procedure DElement2!output the data in the DTensor of inde
		module procedure DTiI!output the index DTensor in the DTensorlink,index is a vector
	end interface
	
	interface operator(.dim.)
		module procedure DgetTenDim_i!output the i dimension of the DTensor
	end interface
		
	interface assignment(=)
		module procedure DTcopy!T1=T2 ,both T1 and T2 are DTensor
		module procedure DTcopy_dim1!Vec=T ,Vec is a vector and T are DTensor,vec=T%DTensor_Data,used for any case
		module procedure DTcopy_dim2!Mat=T ,Mat is a matrix and T are DTensor,Mat=T%DTensor_Data,used only rank=2
		module procedure DTcopy_dim3!Mat=T ,Mat is a rank=3 and T are DTensor,Mat=T%DTensor_Data,used only rank=3
		module procedure DassignmentTen1
		module procedure DassignmentTen2
		module procedure DassignmentTen3
		module procedure DTcopylink!link1=link2,both link1 and link2 are DTensorLink
		module procedure Dcom_Ten!val=T,val is real*8 ,and T is a DTensor,used only there is one DElement in T
		module procedure DTenLink!Tensor=Tensorlink,if every element(Tensor) in the link is a one-element Tensor
	end interface
!**********************************************************
!	Other function or subroutine:
!
!		DcleanDTensor: clean all the data in type(DTensor)
!
!		Print DTensor: 
!			DTprint(DTensor):print all the informaDTion of the DTensor
!			DTDprint(DTensor):print Dimension
!			DTMprint(Tesnor):Print as matrix 
!
!		get rank or totalData or flag
!			DgetRank
!			DgetTotalData
!			DgetFlag
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
!
!		DTensorProduct1Dim(T1,T2): Dot product,return a comple*16
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
!		diagonal matrix:
!			Deye_com(s,m,n):return a m*n diagonal matrix will diag(s(1),s(2),...s(slen),0,0....),
!				s(:) are real*8
!			Deye_real::return a m*n diagonal matrix will diag(s(1),s(2),...s(slen),0,0....),
!				s(:) are real*8
!			Done_com(m,n)  :return a m*n idenDTity matrix
!
!		Dexpm(H):
!			return a DTensor e^H,
!
!		Deng(H,vec,val):
!			return a eigenvalues matrix of H,H=vec*val*vec^H
!			
!		operator value:
!			Dvalue(F,Tr):return the Dvalue of <phi|F|phi>,where F is an operator,
!				F=<I',s1',s2',..|F|I,s1,s2,..>
!			Dnorm2(Tr_):return  <phi|phi>
!
!		DTensor1(dimen,num):return a DTensor with all the elememt num,dimen is dimension of
!			the DTensor
!		
!		DTen_Ten(T1_,T2_,i1,i2,ifDdecompose):Dcontract The i1 index of T1 and the i2 index of t2
!			T1: (D1,D2,..D_i1,...),T2 :(E1,E2,...,E_i2,....),then the result will be 
!			T:(D1,D2,..D_{i1-1},D_{i1+1},...,D_{rank1},E1,E2,...,E_{i2-1},E_{i2+1},...,E_{rank2})	
!			if if Ddecompose=1, T will be Ddecompose
!
!		Dadd A DTensor to a link
!			Dpush_back(h,T):Dadd a DTensor to the end of the link,there is no vector index in the input DTensor
!			Dpush_backI(h,T,inde):Dadd a DTensor to the begin of the link inde is the index of the DTensor
!			Dpush_fo(h,T):Dadd a DTensor to the begin of the link there is no vector index in the  input DTensor
!			Dpush_foI(h,T):Dadd a DTensor to the begin of the link inde is the index of the DTensor	
!
!		DTensor index
!			DTenIndice(h,inde,indice):return the index of the inde DTensor.T_{1,1}->T_{1,2}->T_{1,3}
!				->T_{2,1}->T_{2,2},if input inde=2 => indice=[1,2].on entry the size of indice should
!				 be equal to the one in h
!			DTenIndiceLog(h,inde,indice):if there no index in h,return .false..other while return inde
!				of the DTensor
!
!		Dchecklength(link):check and output length of the link,check if the length of the link is equal 
!				to the Dvalue in head of the link
!
!		Dlinklength(link):output length of the link
!
!		Dmodify(h,T,inde):Dmodify the inde DElement of the link.The DTensor will be Dmodify
!
!		Dcleanlink(h):clean the link h
!
!		Ddeletelink(link_in,inde1,inde2):
!			link   :	T1->T2->...->T_{inde1-1}->T_inde1->T_{inde1+1}->....->T_{inde2-1}->T_inde2->T_{inde2+1}...->TN 
!		  result:		T1->T2->...->T_{inde1-1}->T_{inde1+1}....->T_{inde2-1}->->T_{inde2+1}...->TN 
!		  delete the inde1 to inde2 DTensor in the link,note that the deleted DTensor are include the DTensors of inde1 and inde2
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
!				A=B before allocaDTing storage space for A,A will get
!			 	no data.This happen fortran90 but do not in fortran77.
!				assignments is used for A=B
!
!**********************************************************	
	contains
	
!***************** assignment *********************
	subroutine DTcopy(T,T2)
		type(DTensor),intent(out) ::T
		type(DTensor),intent(in) :: T2
		if(T2%flag) then
			T%rank=T2%rank
			T%TenDim=T2%TenDim
	!		call assignment_com_dim1(T%DTensor_Data,T2%DTensor_Data)
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
	subroutine DassignmentTen1(T,Tensor_data)
		real*8,intent(in)::Tensor_data(:)
		type(DTensor),intent(inout)::T
		call DcleanDTensor(T)
		T=DbuildTen(Tensor_data)
		return
	end subroutine
	subroutine DassignmentTen2(T,Tensor_data)
		real*8,intent(in)::Tensor_data(:,:)
		type(DTensor),intent(inout)::T
		call DcleanDTensor(T)
		T=DbuildTen(Tensor_data)
		return
	end subroutine
	subroutine DassignmentTen3(T,Tensor_data)
		real*8,intent(in)::Tensor_data(:,:,:)
		type(DTensor),intent(inout)::T
		call DcleanDTensor(T)
		T=DbuildTen(Tensor_data)
		return
	end subroutine
!*************** DcleanDTensor  *****************
	subroutine DcleanDTensor(T)
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
!DTcopy the DTensor_Data to a vector	
	subroutine DTcopy_dim1(Vec,T)
		real*8,allocatable,intent(out) ::Vec(:)
		type(DTensor),intent(in) :: T
		integer::length
		length=T%TotalData
		if(allocated(Vec)) then
			deallocate(Vec)
		end if
		allocate(Vec(length))
		call Dcopy(length,T%DTensor_Data,1,Vec,1)
		return
	end subroutine
!if rank=2,then DTcopy the DTensor_Data to a mat
! mat is a matrix	
	subroutine DTcopy_dim2(Mat,T)
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
			call Dcopy(m*n,T%DTensor_Data,1,Mat,1)
!			Mat=reshape(T%DTensor_Data,(/m,n/))
		else
			write(*,*)"The DTensor is not a matrix"
			write(*,*)"Assignment ERROR,stop"
			stop
		end if
		return
	end subroutine
!if rank=3,then DTcopy the DTensor_Data to a mat
!mat is a 3 dimension array	
	subroutine DTcopy_dim3(Mat,T)
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
			call Dcopy(m*n*l,T%DTensor_Data,1,Mat,1)
!			Mat=reshape(T%DTensor_Data,(/m,n,l/))
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
		call Dcopy(length,inData,1,T%DTensor_Data,1)
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
		call Dcopy(length,inData,1,T%DTensor_Data,1)
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
		call Dcopy(length,inData,1,T%DTensor_Data,1)
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
!		call assignment_com_dim1(DbuildTen1%DTensor_Data,DTensor_Data)
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
!		call assignment_com_dim1(DbuildTen2%DTensor_Data,DTensor_Data)
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
!		T%DTensor_data=Ten_data
!		call assignment_com_dim1(T%DTensor_data,Ten_data)
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
!*********************  DgetFlag	 **********************
	logical function DgetFlag(T)
		type(DTensor),intent(in) :: T
		DgetFlag=T%flag
	end function		
!*********************  Dgenerate *********************
!		Dgenerate a DTensor with random number
	type(DTensor) function Dgenerate(Tdim) result(T)
		integer,intent(in) :: Tdim(:)
		integer :: rank,totalData,i
		real*8,allocatable :: DTensor_data(:)
		real*8 ::temp_real
		type(Dimension):: TenDim
		call Dinit_random_seed()
		rank=size(Tdim)
		totalData=product(Tdim)
		allocate(DTensor_data(totalData))
		do i=1,totalData
			call random_number(temp_real)
			DTensor_data(i)=temp_real
		end do
		TenDim=Tdim
!		T%DTensor_data=DTensor_data
!		call assignment_com_dim1(T%DTensor_data,DTensor_data)
		call DstoreTenData(T,DTensor_data)
		T%rank=rank
		T%totalData=totalData
		T%TenDim=TenDim
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
	subroutine DTprint(T)
		type(DTensor) :: T
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
!********************* print_Matrix *********************
	subroutine DTMprint(T)
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


!cccccccccccccccc Dadd          cccccccccccccccccc
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
!		call assignment_com_dim1(Ddivide_number%DTensor_data,T1%DTensor_data/num)
		call DstoreTenData(Ddivide_number,T1%DTensor_data/num)
		Ddivide_number%totalData=T1%totalData
		Ddivide_number%flag=.true.
		return
	end function
!ccccccccccccccccc Ddivide_number          cccccccccccccccccc	
	type(DTensor) function Ddivide_numberReal(T1,num)
		type(DTensor),intent(in) :: T1
		real*8,intent(in) :: num
		Ddivide_numberReal%rank=T1%rank
		Ddivide_numberReal%TenDim=T1%TenDim
!		call assignment_com_dim1(Ddivide_numberReal%DTensor_data,T1%DTensor_data/num)
		call DstoreTenData(Ddivide_numberReal,T1%DTensor_data/num)
		Ddivide_numberReal%totalData=T1%totalData
		Ddivide_numberReal%flag=.true.
		return
	end function
!ccccccccccccccccc Ddivide_number          cccccccccccccccccc	
	type(DTensor) function Ddivide_numberReal4(T1,num)
		type(DTensor),intent(in) :: T1
		real*4,intent(in) :: num
		Ddivide_numberReal4%rank=T1%rank
		Ddivide_numberReal4%TenDim=T1%TenDim
!		call assignment_com_dim1(Ddivide_numberReal4%DTensor_data,T1%DTensor_data/num)
		call DstoreTenData(Ddivide_numberReal4,T1%DTensor_data/num)
		Ddivide_numberReal4%totalData=T1%totalData
		Ddivide_numberReal4%flag=.true.
		return
	end function			
!cccccccccccccccc DmulDTiply_DTensor          cccccccccccccccccc	
	type(DTensor) function DmulDTiply_DTensor(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		type(Dimension)::dim1,dim2
		dim1=T1%TenDim
		dim2=T2%TenDim
		if(.not.(dim1.equ.dim2)) then
			write(*,*) "The dimension of T1 and T2 are not the same"
			write(*,*)"The program will stop"
			stop
		end if
		DmulDTiply_DTensor=Dset(dim1,T1%DTensor_data*T2%DTensor_data)
		return
	end function
!cccccccccccccccc DmulDTiply_number          cccccccccccccccccc	
	type(DTensor) function DmulDTiply_number(T1,num)
		type(DTensor),intent(in) :: T1
		real*8,intent(in) ::   num
		DmulDTiply_number%rank=T1%rank
		DmulDTiply_number%TenDim=T1%TenDim
!		call assignment_com_dim1(DmulDTiply_number%DTensor_data,T1%DTensor_data*num)
		call DstoreTenData(DmulDTiply_number,T1%DTensor_data*num)
		DmulDTiply_number%DTensor_data=T1%DTensor_data*num
		DmulDTiply_number%totalData=T1%totalData
		DmulDTiply_number%flag=.true.
		return
	end function
		type(DTensor) function DmulDTiply_number_(num,T1)
		type(DTensor),intent(in) :: T1
		real*8,intent(in) ::   num
		DmulDTiply_number_%rank=T1%rank
		DmulDTiply_number_%TenDim=T1%TenDim
!		call assignment_com_dim1(DmulDTiply_number%DTensor_data,T1%DTensor_data*num)
		call DstoreTenData(DmulDTiply_number_,T1%DTensor_data*num)
		DmulDTiply_number_%DTensor_data=T1%DTensor_data*num
		DmulDTiply_number_%totalData=T1%totalData
		DmulDTiply_number_%flag=.true.
		return
	end function
!cccccccccccccccc DmulDTiply_real          cccccccccccccccccc	
	type(DTensor) function DmulDTiply_real(T1,num)
		type(DTensor),intent(in) :: T1
		real*8,intent(in) ::   num
		DmulDTiply_real%rank=T1%rank
		DmulDTiply_real%TenDim=T1%TenDim
!		call assignment_com_dim1(DmulDTiply_real%DTensor_data,T1%DTensor_data*num)
		call DstoreTenData(DmulDTiply_real,T1%DTensor_data*num)
		DmulDTiply_real%totalData=T1%totalData
		DmulDTiply_real%flag=.true.
		return
	end function
	type(DTensor) function DmulDTiply_real_(num,T1)
		type(DTensor),intent(in) :: T1
		real*8,intent(in) ::   num
		DmulDTiply_real_%rank=T1%rank
		DmulDTiply_real_%TenDim=T1%TenDim
!		call assignment_com_dim1(DmulDTiply_real%DTensor_data,T1%DTensor_data*num)
		call DstoreTenData(DmulDTiply_real_,T1%DTensor_data*num)
		DmulDTiply_real_%totalData=T1%totalData
		DmulDTiply_real_%flag=.true.
		return
	end function	
!ccccccccccccccccc DmulDTiply_real4          cccccccccccccccccc	
	type(DTensor) function DmulDTiply_real4(T1,num)
		type(DTensor),intent(in) :: T1
		real*4,intent(in) ::   num
		DmulDTiply_real4%rank=T1%rank
		DmulDTiply_real4%TenDim=T1%TenDim
!		call assignment_com_dim1(DmulDTiply_real4%DTensor_data,T1%DTensor_data*num)
		call DstoreTenData(DmulDTiply_real4,T1%DTensor_data*num)
		DmulDTiply_real4%totalData=T1%totalData
		DmulDTiply_real4%flag=.true.
		return
	end function			
!ccccccccccccccccc Dpermute_rank3   cccccccccccccccccc		
	type(DTensor) function Dpermute_rank3(T1,index_not_permute)
		type(DTensor),intent(in) :: T1
		integer,intent(in) ::   index_not_permute
		integer ::i,newdimen(3),totaldata
		integer,allocatable :: dimen(:)
		type(Dimension) ::newTDim
		real*8,allocatable :: Tdata(:,:,:),newdata(:,:,:)
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
!			call assignment_com_dim1(Dpermute_rank3%DTensor_data,reshape(newdata,(/T1%totalData/)))
			call DstoreTenData(Dpermute_rank3,newdata)
			Dpermute_rank3%rank=3
			Dpermute_rank3%totalData=T1%totalData
			Dpermute_rank3%TenDim=newTDim
			Dpermute_rank3%flag=.true.
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
!			call assignment_com_dim1(Dpermute_rank3%DTensor_data,reshape(newdata,(/T1%totalData/)))
			call DstoreTenData(Dpermute_rank3,newdata)
			Dpermute_rank3%rank=3
			Dpermute_rank3%totalData=T1%totalData
			Dpermute_rank3%TenDim=newTDim
			Dpermute_rank3%flag=.true.
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
!			call assignment_com_dim1(Dpermute_rank3%DTensor_data,reshape(newdata,(/T1%totalData/)))
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
!		call assignment_com_dim1(Dpermute_rank2%DTensor_data,reshape(newdata,(/T%totalData/)))
		call DstoreTenData(Dpermute_rank2,newdata)
		Dpermute_rank2%rank=2
		Dpermute_rank2%totalData=T%totalData
		Dpermute_rank2%TenDim=newTDim
		Dpermute_rank2%flag=.true.
	end function
!cccccccccccccccc permutation   cccccccccccccccccc
	type(DTensor) function Dpermutation(T,newOrder)
		type(DTensor),intent(in) :: T
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
			dpermutation=T
			return
		end if
		total=T%TotalData
		Dpermutation%TenDim=Dimpermute(T%TenDim,newOrder)
		allocate(Dpermutation%DTensor_Data(total))
		maxdim=T%TenDim
		newmaxdim=Dpermutation%TenDim
		lendim=size(maxdim)
		allocate(newinde(lendim))
		i=1
		do i=1,total
			call dIndesToaddress(maxdim,inde,i)
			do j=1,lendim
				newinde(j)=inde(newOrder(j))
			end do
			newaddress=DaddressToIndes2(newmaxdim,newinde)
			Dpermutation%DTensor_Data(newaddress)=T%DTensor_Data(i)
		end do
		Dpermutation%rank=T%rank
		Dpermutation%totalData=total
		Dpermutation%flag=T%flag
	end function
!ccccccccccccccccc   DproductTen   ccccccccccccccccccc
!			the dimension of T1 and T2 should be 2 or 1
!			choose the righ form of
!						DTensorProduct
!					DTensorProduct1Dim
!					DTensorProduct1Dim2
	type(DTensor)	function DproductTen(T1,T2)
		type(DTensor),intent(in) :: T1
		type(DTensor),intent(in) :: T2
		integer :: rank1,rank2
		character*50 ::wt
		rank1=T1%rank
		rank2=T2%rank
		if((rank1.le.2).and.(rank2.le.2)) then
!			wt="T1,in DproductTen"
!			call warning(T1,wt)
!			wt="T2,in DproductTen"
!			call warning(T2,wt)
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
			write(*,*) "error in the dimension of the DTensor"
		end if
		return
	end function
				 	
!cccccccccccccccc DTensorProduct   cccccccccccccccccc			 
!		both T1 and T2 should be a 2 \DTimes 2 DTensor
	type(DTensor) function DTensorProduct(T1,T2)
		type(DTensor),intent(in) :: T1,T2
!		real*8,allocatable :: Ndata(:,:)
!		real*8,allocatable :: T1data(:,:)
!		real*8,allocatable :: T2data(:,:)
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
!			call CPU_DTiME(cpuDTime(1))
			call dGEMM('N','N',T1m,T2n,T1n,dcmplx(1d0,0d0),T1%DTensor_Data,T1m,&
				T2%DTensor_Data,T2m,0d0,DTensorProduct%DTensor_data,T1m)
!			call CPU_DTiME(cpuDTime(2))
!			write(*,*)"DTime2",cpuDTime(2)-cpuDTime(1)
			dim1=T1%TenDim.sub.1
			dim2=T2%TenDim.sub.2
			newdim=dim1+dim2
!			call assignment_com_dim1(DTensorProduct%DTensor_data,reshape(Ndata,(/totalData/)))
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
			Ndata=DDOT(T1n,T1%DTensor_Data,1,T2%DTensor_Data,1)
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
			Ndata=DDOT(T1n,T1%DTensor_Data,1,T2%DTensor_Data,1)
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
			call DGEMV('T',M,N,dcmplx(1d0,0d0),T2%DTensor_Data,M,T1%DTensor_Data,1,0d0,reData,1)
			ResultT=Dset(T2%TenDim.sub.2,reData)
		else if((T1%rank.eq.2).and.(T2%rank.eq.1)) then
			M=T1.dim.1
			N=T1.dim.2
			if((T2.dim.1) .ne. N) then
				write(*,*)"ERROR in DTensorProduct1Dim2,stop"
				stop
			end if
			allocate(reData(M))
			call ZGEMV('N',M,N,dcmplx(1d0,0d0),T1%DTensor_Data,M,T2%DTensor_Data,1,0d0,reData,1)
			ResultT=Dset(T1%TenDim.sub.1,reData)
		else
			write(*,*)"Error in the dimension of the input DTensors"
		end if
		return
	end function
			
!**************** Dcontract   ****************
!		combine two index of the DTensor,which is con_index and con_index+1
	type(DTensor) function Dcontract(T1,con_index)
		integer,intent(in) :: con_index
		type(DTensor),intent(in) :: T1
		type(Dimension) ::newdim
		newdim=DimConstract(T1%TenDim,con_index,1)
!		call assignment_com_dim1(Dcontract%DTensor_data,T1%DTensor_Data)
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
!		call assignment_com_dim1(Dcontracts%DTensor_data,T1%DTensor_Data)
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
!		call assignment_com_dim1(Ddecompose%DTensor_data,T1%DTensor_Data)
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
		newDim=DimDecompose(T1%TenDim,de_index,inde)
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
!		call assignment_com_dim1(Ddecompose1%DTensor_data,T1%DTensor_Data)
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
!		call assignment_com_dim1(DdecomposeAll%DTensor_data,T1%DTensor_Data)
		call DstoreTenData(DdecomposeAll,T1%DTensor_Data)
		DdecomposeAll%rank=DimSize(newDim)
		DdecomposeAll%totalData=T1%totalData
		DdecomposeAll%TenDim=newdim
		DdecomposeAll%flag=.true.
		return
	end function
!	Operar is a matrix,the first element of every row specify
!		1:contracts
!		2:decompose
!		3:decomposeAll
!	other element of every row is input parameter of the function

	type(DTensor) function Dcompose_decompse(T1,Operar)
		integer,intent(in)::Operar(:,:)
		type(DTensor),intent(in) :: T1
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
				Dcompose_decompse%TenDim=DimConstract(T1%TenDim,index1,num)
				Dcompose_decompse%rank=DimSize(T1%TenDim)
			else if(Operar(i,1).eq.2) then
				de_index=Operar(i,2)
				inde=Operar(i,3)
				Dcompose_decompse%TenDim=DimDecompose(T1%TenDim,de_index,inde)
				Dcompose_decompse%rank=DimSize(T1%TenDim)
			else if(Operar(i,1).eq.3) then
				Dcompose_decompse%TenDim=DimDecomposeAll(T1%TenDim)
				Dcompose_decompse%rank=DimSize(T1%TenDim)
			else
				write(*,*) "error in compose_decompse"
				stop
			end if
		end do
		call DstoreTenData(Dcompose_decompse,T1%DTensor_Data)
		Dcompose_decompse%flag=.true.
		Dcompose_decompse%totalData=T1%totalData
		return
	end function
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
	real*8 function maxDElement(T)
		type(DTensor),intent(in) :: T
		maxDElement=maxval(abs(T%DTensor_Data))
		return
	end function
!*****************  maxRealDElement  *****************
	real*8 function maxRealDElement(T)
		type(DTensor),intent(in) :: T
		maxRealDElement=maxval(T%DTensor_Data)
		return
	end function
!*****************  minDElement  *****************
	real*8 function minDElement(T)
		type(DTensor),intent(in) :: T
		minDElement=minval(abs(T%DTensor_Data))
		return
	end function
!*****************  minRealDElement  *****************
	real*8 function minRealDElement(T)
		type(DTensor),intent(in) :: T
		minRealDElement=minval(T%DTensor_Data)
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
		real*8::cpuDTime(2)
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
!		real*8,allocatable :: newT1data(:,:)
		real*8,allocatable :: Tdata(:)
		real*8,allocatable :: work(:)
		integer m,n,info,lw,rwdim,i,j
		integer,allocatable :: iw(:)
		real*8,allocatable :: rw(:),sdata(:)
		integer :: ns_max,inde_ij,inde_ij2
		real*8::cpuDTime(2)
		type(Dimension) :: T1Dim,T2Dim
!		call NANwarning2(T)
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
!			call CPU_DTiME(cpuDTime(1))
			call ZGESDD('S',m,n,Tdata,m,sdata,T1data,m,T2data,n,WORK,lw,rw,iw,INFO)
!			call CPU_DTiME(cpuDTime(2))
!				write(*,*)"DTime",cpuDTime(2)-cpuDTime(1)
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
!			allocate(newT1data(m,ns_max))
			do j=1,ns_max
				do i=1,m
					inde_ij=DaddressToIndes2((/m,n/),(/i,j/))
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
						inde_ij=DaddressToIndes2((/n,n/),(/i,j/))
						inde_ij2=DaddressToIndes2((/ns_max,n/),(/i,j/))
						newT2data(inde_ij)=T2data(inde_ij2)
					end do
				end do
				!newT2data=T2data(1:ns_max,:)
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
!*****************  Deye  *****************
	type(DTensor) function Deye_com(s,m,n)
		real*8,intent(in) :: s(:)
		integer,intent(in) :: m,n
		real*8,allocatable :: sdata(:,:)
		integer :: i,j,lens
		lens=size(s)
		allocate(sdata(m,n))
		do j=1,n
			do i=1,m
				if(i.eq.j) then
					if(lens.lt.i) then
						sdata(i,j)=0D0
					else
						sdata(i,j)=s(i)
					end if
				else
					sdata(i,j)=0
				end if
			end do
		end do
!		call assignment_com_dim1(Deye_com%DTensor_data,reshape(sdata,(/m*n/)))
		call DstoreTenData(Deye_com,sdata)
		Deye_com%rank=2
		Deye_com%totalData=m*n
		Deye_com%TenDim=(/m,n/)
		Deye_com%flag=.true.
		return
	end function
!*****************  Deye_com  *****************
	type(DTensor) function Deye_real(s,m,n)
		real*8,intent(in) :: s(:)
		integer,intent(in) :: m,n
		real*8,allocatable :: sdata(:,:)
		integer :: i,j,lens
		lens=size(s)
		allocate(sdata(m,n))
		do j=1,n
			do i=1,m
				if(i.eq.j) then
					if(lens.lt.i) then
						sdata(i,j)=0d0
					else
						sdata(i,j)=s(i)
					end if
				else
				sdata(i,j)=0d0
				end if
			end do
		end do
!		call assignment_com_dim1(Deye_real%DTensor_data,reshape(sdata,(/m*n/)))
		call DstoreTenData(Deye_real,sdata)
		Deye_real%rank=2
		Deye_real%totalData=m*n
		Deye_real%TenDim=(/m,n/)
		Deye_real%flag=.true.
		return
	end function
!*****************  Done_com  *****************
	type(DTensor) function Done_com(m,n)
		integer,intent(in) :: m,n
		real*8,allocatable :: sdata(:,:)
		integer :: i,j
		allocate(sdata(m,n))
		do j=1,n
			do i=1,m
				if(i.eq.j) then
					sdata(i,j)=1d0
				else
					sdata(i,j)=0d0
				end if
			end do
		end do
!		call assignment_com_dim1(Done_com%DTensor_data,reshape(sdata,(/m*n/)))
		call DstoreTenData(Done_com,sdata)
		Done_com%rank=2
		Done_com%totalData=m*n
		Done_com%TenDim=(/m,n/)
		Done_com%flag=.true.
		return
	end function	
!****************  Done_Ten  *************
	type(DTensor) function Done_Ten(dimen)
		integer,intent(in) :: dimen(:)
		integer :: i,total,lenDim,min_dim,Daddress,dim_i
		real*8,allocatable::TenData(:)
		total=product(dimen)
		lenDim=size(dimen)
		min_dim=minval(dimen)
		allocate(TenData(total))
		TenData=0
		do dim_i=1,min_dim
			Daddress=0
			do i=lenDim,2,-1
				Daddress=Daddress+(dim_i-1)*product(dimen(1:(i-1)))
			end do
			Daddress=Daddress+dim_i
			TenData(Daddress)=1d0
		end do
		Done_Ten=DbuildTen(dimen,TenData)
		return
	end function	
! return exp(H)
	type(DTensor) function Dexpm(H)
		type(DTensor),intent(in) ::H
		real*8,allocatable :: Hdata(:,:),work(:)
		real*8,allocatable :: expHdata(:,:)
		real*8,allocatable :: DengVal(:),DengVec(:,:)
		integer :: hdim(2),i,j,N,lwork,info,SDIM
		real*8,allocatable::rwork(:)
		logical,allocatable::bwork(:)
		type(dimension)::newHdim
		hdim(1)=H.dim.1
		hdim(2)=H.dim.2
		newHdim=H%TenDim
		if(hdim(1).ne.hdim(2)) then
			write(*,*)"ERROR in Dexpm"
			stop
		end if
		N=hdim(1)
		lwork=2*N
		
		allocate(Hdata(N,N))
		allocate(DengVal(N))
		allocate(DengVec(N,N))
		allocate(rwork(N))
		allocate(work(lwork))
		allocate(bwork(N))
		allocate(expHdata(N,N))
		Hdata=H
		call ZGEES('V','N',1,N,Hdata,N,SDIM,DengVal,DengVec,N,WORK,LWORK,RWORK,BWORK,INFO)
		do i=1,N
			DengVal(i)=exp(DengVal(i))
		end do
		do i=1,N
			do j=1,N
				expHdata(i,j)=DengVec(i,j)*DengVal(j)
			end do
		end do
		DengVec=transpose(DengVec)
		expHdata=matmul(expHdata,DengVec)
		Dexpm=Dset(newHdim,reshape(expHdata,(/N*N/)))
		return
	end function
!return eigenDvalue and eigenvectoers		
	subroutine Deng(H,vec,val)
		type(DTensor),intent(in) ::H
		type(DTensor),intent(inout) ::val,vec
		real*8,allocatable :: Hdata(:,:),work(:)
		real*8,allocatable :: DengVal(:),DengVec(:,:)
		integer :: hdim(2),i,j,N,lwork,info,SDIM
		real*8,allocatable::rwork(:)
		logical,allocatable::bwork(:)
		type(dimension)::newHdim
		hdim(1)=H.dim.1
		hdim(2)=H.dim.2
		newHdim=H%TenDim
		if(hdim(1).ne.hdim(2)) then
			write(*,*)"ERROR in Dexpm"
			stop
		end if
		N=hdim(1)
		lwork=2*N
		
		allocate(Hdata(N,N))
		allocate(DengVal(N))
		allocate(DengVec(N,N))
		allocate(rwork(N))
		allocate(work(lwork))
		allocate(bwork(N))
		Hdata=H
		call ZGEES('V','N',1,N,Hdata,N,SDIM,DengVal,DengVec,N,WORK,LWORK,RWORK,BWORK,INFO)
		val=Deye(DengVal,N,N)
		vec=Dset((/N,N/),reshape(DengVec,(/N*N/)))
		return
	end subroutine	
!*****************  sqrt  *********************
! return sqrt(A),A is a matrix

	type(DTensor) function DsqrtTen(h)
		type(DTensor),intent(in) ::H
		type(DTensor)::val,vec
		real*8,allocatable :: Hdata(:,:),work(:)
		real*8,allocatable :: DengVal(:),DengVec(:,:)
		integer :: hdim(2),i,j,N,lwork,info,SDIM
		real*8,allocatable::rwork(:)
		logical,allocatable::bwork(:)
		type(dimension)::newHdim
		hdim(1)=H.dim.1
		hdim(2)=H.dim.2
		newHdim=H%TenDim
		if(hdim(1).ne.hdim(2)) then
			write(*,*)"ERROR in Dexpm"
			stop
		end if
		N=hdim(1)
		lwork=2*N
		
		allocate(Hdata(N,N))
		allocate(DengVal(N))
		allocate(DengVec(N,N))
		allocate(rwork(N))
		allocate(work(lwork))
		allocate(bwork(N))
		Hdata=H
		call ZGEES('V','N',1,N,Hdata,N,SDIM,DengVal,DengVec,N,WORK,LWORK,RWORK,BWORK,INFO)
		do i=1,N
			if(DengVal(i).ge.0) then
				DengVal(i)=dsqrt(DengVal(i))
			else
				write(*,*)"ERROR in Dsqrt of Tensor"
				write(*,*)DengVal(i)
				stop
			end if
		end do
		val=Deye(DengVal,N,N)
		vec=Dset((/N,N/),reshape(DengVec,(/N*N/)))
		DsqrtTen=vec * val * (.H. vec)
	end function	
!*****************  Dtranspose  *****************
	type(DTensor) function Dtranspose(T)
		type(DTensor),intent(in) :: T
		integer::rank
		rank=T%rank
		if(rank.eq.1) then
			Dtranspose=Dconjugate(T)
		else if(rank.eq.2) then
			Dtranspose=DconjugateTranspose(T)
		else
			write(*,*) "DTensor should be 1 or 2 dimension"
		end if
		return
	end function
!*****************  DconjugateTranspose  *****************
	type(DTensor) function DconjugateTranspose(T)
		type(DTensor),intent(in) :: T
		real*8,allocatable :: Tdata(:,:),newdata(:,:)
		integer :: m,n
		type(Dimension) :: Dimen
		if(T%rank.eq.2) then
			m=T.dim.1
			n=T.dim.2
			Dimen=(T%TenDim.sub.2)+(T%TenDim.sub.1)
			allocate(Tdata(m,n))
			allocate(newdata(n,m))
			Tdata=T
			newdata=transpose(Tdata)
!			DconjugateTranspose=Dset(Dimen,reshape(newdata,(/m*n/)))
			DconjugateTranspose%TenDim=Dimen
			DconjugateTranspose%rank=2
			DconjugateTranspose%flag=.true.
			DconjugateTranspose%TotalData=m*n
			call DstoreTenData(DconjugateTranspose,newdata)			
		else
			write(*,*) "The input DTensor should 1 or 2 dimension"
		end if
		return
	end function
!*****************  Dconjugate  *****************
	type(DTensor) function Dconjugate(T)
		type(DTensor),intent(in) :: T
		real*8,allocatable :: newdata(:)
		integer :: m,n
		type(Dimension) :: Dimen
		integer::rank
		m=T.dim. 1
		Dimen=T%TenDim
		allocate(newdata(m))
		newdata=T%DTensor_Data
		Dconjugate=Dset(Dimen,newdata)
		return
	end function
!*****************  Dvalue  *****************
!			return the Dvalue of <phi|F|phi>,where F is an operator
!			F=<I',s1',s2',..|F|I,s1,s2,..>
	real*8 function Dvalue(F_,Tr_)
		type(DTensor),intent(in) :: F_,Tr_
		type(DTensor) :: Tr,Tl,F,midResult
		integer :: i,rank
		rank=Tr_%rank
		F=F_
		Tr=Tr_
		do i=1,rank-1
			Tr=Tr.c.1
			F=F.c.1
		end do
		do i=1,rank-1
			F=F.c.2
		end do
		Tl=.h.Tr
		Dvalue=Tl*F*Tr
		return
	end function				
!*****************  norm   *****************
!			return  <phi|phi>
	real*8 function Dnorm2(Tr_)
		type(DTensor),intent(in) :: Tr_
		type(DTensor) :: Tr
		real*8 :: resul
		Tr=Dcontracts(Tr_,1,Tr_%rank)! Tr is [(I*s1*s2*s3*...*sn),1] dimension
		Dnorm2=DNRM2(Tr%TotalData,Tr%DTensor_Data,1)
		Dnorm2=Dnorm2*Dnorm2
		return
	end function	
!		return  sqrt(<phi|phi>	)
	real*8 function Dnorm(Tr_)
		type(DTensor),intent(in) :: Tr_
		type(DTensor) :: Tr
		Tr=Dcontracts(Tr_,1,Tr_%rank)! Tr is [(I*s1*s2*s3*...*sn),1] dimension
		Dnorm=DNRM2(Tr%TotalData,Tr%DTensor_Data,1)
		return
	end function		
!*****************  DaddressToIndes   *****************
	integer function DaddressToIndes(T,Adim)
		type(DTensor),intent(in) :: T
		integer,intent(in) :: Adim(:)
		integer,allocatable::Tdim(:)
		integer :: i,Dimlen
		
		
		Tdim=T%TenDim
		Dimlen=size(TDim)
		if(Dimlen.eq.1) then
			DaddressToIndes=Adim(1)
			return
		end if
		DaddressToIndes=0
		do i=Dimlen,2,-1
			DaddressToIndes=DaddressToIndes+(Adim(i)-1)*product(TDim(1:(i-1)))
		end do
		DaddressToIndes=DaddressToIndes+Adim(1)
		return 
		
! 		Dmodify the code below in 2013.5.18		
!		DaddressToIndes=Adim(1)
!		Tdim=T%TenDim
!		do i=2,size(Adim)
!			DaddressToIndes=DaddressToIndes+(Adim(i)-1)*Tdim(i-1)
!		end do
!		return
	end function			
	integer function DaddressToIndes2(TDim,Adim)
		integer,intent(in) :: TDim(:)!max dimension
		integer,intent(in) :: Adim(:)
		integer :: i,Dimlen
		Dimlen=size(TDim)
		if(Dimlen.eq.1) then
			DaddressToIndes2=Adim(1)
			return
		end if
		DaddressToIndes2=0
		do i=Dimlen,2,-1
			DaddressToIndes2=DaddressToIndes2+(Adim(i)-1)*product(TDim(1:(i-1)))
		end do
		DaddressToIndes2=DaddressToIndes2+Adim(1)
		return 
! 		Dmodify the code below in 2013.5.18
!		DaddressToIndes2=Adim(1)
!		do i=2,size(Adim)  
!			DaddressToIndes2=DaddressToIndes2+(Adim(i)-1)*Tdim(i-1)
!		end do
		return
	end function
	subroutine DIndesToaddress(TDim,Adim,inde)
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
		inde=DaddressToIndes(T,Tdim)
		DElement=T%DTensor_Data(inde)
		return
	end function		  
	real*8 function DElement2(T,inde)
		type(DTensor),intent(in) ::T
		integer,intent(in)::inde
		DElement2=T%DTensor_Data(inde)
		return
	end function		
!***************  DdirectProduct  *********************
	type(DTensor) function DdirectProduct(T1,T2)
		type(DTensor),intent(in) :: T1,T2
		real*8,allocatable :: T1data(:,:),T2data(:,:)
		real*8,allocatable :: newdata(:,:,:,:)
		integer :: m1,n1,m2,n2,i,j,k,l,rank(2)
		type(Dimension):: D1,Dtemp
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
			write(*,*)"The directProduct is just for matrix"
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
			write(*,*)"The directProductM is just for matrix"
			stop
		end if	
		return
	end function
!************* directProduct1 **************************
!
!       for two Tensor whose ranks are 1		
!
	type(dTensor) function DdirectProduct1V(T1,T2)
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
			write(*,*)"The directProduct1 is just for vector"
			stop
		end if
		return
	end function	
!*******************  Ddirect sum  **********************
	type(DTensor) function  Ddirect_sum(T1,T2) result(T)
		type(DTensor),intent(in) :: T1,T2
		type(DTensor) ::temp1,temp2
		type(Dimension) ::dimT1,dimT2
		integer :: m,n
		dimT1=T1%TenDim
		dimT2=T2%TenDim
		if(dimT1.equ.dimT2) then
			temp1=T1.mxx.Done_com(m,n)
			temp2=Done_com(m,n).mxx.T2
			T=temp1+temp2
			T=.dc.T
		else
			write(*,*) "error in Ddirect_sumM"
		end if
		return 
	end function
!****************  Ddirect sum returning matrix *******************
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
			temp1=T1.mxx.Done_com(m,n)
			temp2=Done_com(m,n).mxx.T2
			T=temp1+temp2
		else
			write(*,*) "error in Ddirect_sumM"
		end if
		return 
	end function
!**************   conbine   ********************
!
!                    / T1 \
!	conbine(T1,T2)	=  |----|
!		               \ T2 /
!	T1 :a [m,n] matrix
!	T2 :a [m,n] matrix
!	conbine(T1,T2):a [m,n,2] matrix
!
	type(DTensor) function Dconbine(T1,T2)
		type(DTensor),intent(in)::T1,T2
		integer,allocatable::dim1(:),dim2(:)
		integer::total,total1,total2,i
		real*8,allocatable::newdata(:)
		type(Dimension)::newDim
		dim1=.subDim.T1
		dim2=.subDim.T2
		if(.not.(dim1.equ.dim2)) then
			write(*,*)"can not conbie two Tensor"
			write(*,*)"stop"
			stop
			return
		end if
		total1=DgettotalData(T1)
		total2=DgettotalData(T2)
		total=total1+total2
		allocate(newdata(total))
		do i=1,total1
			newdata(i)=T1%DTensor_Data(i)
		end do
		do i=1,total2
			newdata(i+total1)=T2%DTensor_Data(i)
		end do
		call DstoreTenData(Dconbine,newdata)
		Dconbine%Rank=DgetRank(T1)+1
		newDim=.subDim.T1
		newDim=newDim+(/2/)
		Dconbine%TenDim=newDim
		Dconbine%totalData=total
		Dconbine%flag=.true.
		return
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
		integer::Daddre
		Daddre=DaddressToIndes(Ten,dimen)
		Ten%DTensor_Data(Daddre)=val
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
		Ten%DTensor_Data=val
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
! modify all the data in the Tensor
	subroutine dresetdim(Ten,dimen)
		type(dTensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		type(Dimension)::dimensi
		if(product(dimen).eq.Ten%TotalData) then
			write(*,*)"ERROR in resetdim"
			write(*,*)"stop"
			stop
		end if
		dimensi=dimen
		Ten%TenDim=dimensi
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
			Tdim1(T1_%rank),Tdim2(T2_%rank),rank2
		rank1=T1_%rank
		rank2=T2_%rank
		T1=T1_
		T2=T2_
		!T1
		if(rank1.ne.1) then
			if(i1.eq.1) then
				T1=Dcontracts(T1,2,rank1)
				T1=.p.T1
			else if (i1.eq.rank1) then
				T1=Dcontracts(T1,1,rank1-1)
			else
				T1=Dcontracts(T1,1,i1-1)
				T1=Dcontracts(T1,3,rank1)
				T1=T1.p.1
				T1=T1.c.1
			end if
		end if
		!T2
		if(rank2.ne.1) then
			if(i2.eq.1) then
				T2=Dcontracts(T2,2,rank2)
			else if (i2.eq.rank2) then
				T2=Dcontracts(T2,1,rank2-1)
				T2=.p.T2
			else
				T2=Dcontracts(T2,1,i2-1)
				T2=Dcontracts(T2,3,rank2)
				T2=T2.p.3
				T2=T2.c.2
			end if
		end if
		T=T1 * T2
		if(ifDdecompose.eq.1) then
			T=.dc.T
		end if
		return
	end function
!*****************permutefo******************************
!		T_{1,2,3,..,i,..,n},permutefo(T,i)=_{i,1,2,3,..,i-1,i+1,..,n}
!
	type(DTensor) function Dpermutefo(T,inde)
		type(DTensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		rank=DgetRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permutefo"
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
			Dpermutefo=.dc.Dpermutefo
			return
		end if
		Dpermutefo=Dcontracts(T,1,inde-1)
		Dpermutefo=Dcontracts(Dpermutefo,3,rank)
		Dpermutefo=Dpermutefo.p.3
		Dpermutefo=.dc.Dpermutefo
		return
	end function
!****************** permuteInde**************************
!		T_{1,2,3,..,i,..,n},permuteInde(T,i)=_{2,3,..,i,1,i+1,..,n}
!
	type(DTensor) function DpermuteInde(T,inde)
		type(DTensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		rank=DgetRank(T)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permuteInde"
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
			DpermuteInde=.dc.DpermuteInde
			return
		end if
		DpermuteInde=Dcontracts(T,2,inde)
		DpermuteInde=Dcontracts(DpermuteInde,3,rank)
		DpermuteInde=DpermuteInde.p.3
		DpermuteInde=.dc.DpermuteInde
		return
	end function

!cccccccccccccccc	pauli_matrix  ccccccccccccccccccccccc
!	the output is sigma_i*num,where i=x,y,z
	subroutine Dpauli_matrix(Tx,Ty,Tz,num)
		type(DTensor) :: Tx,Ty,Tz
		real*8,intent(in) :: num
		REAL*8 :: x(2,2)
		REAL*8 :: y(2,2)
		REAL*8 :: z(2,2)
			x(1,1)=0.
			x(1,2)=1D0
			x(2,1)=1D0
			x(2,2)=0.

			y(1,1)=0.
			y(1,2)=-1D0
			y(2,1)=1D0
			y(2,2)=0.

			z(1,1)=1D0
			z(1,2)=0.
			z(2,1)=0.
			z(2,2)=1D0
			Tx=Dset((/2,2/),reshape(x,(/4/)))
			Ty=Dset((/2,2/),reshape(y,(/4/)))*num
			Tz=Dset((/2,2/),reshape(z,(/4/)))*num
		return
	end subroutine




	 
! 	DTensorlink and DTensornode 	  





!	Dadd a DTensor to the end of the link
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
!	Dadd a DTensor to the begin of the link	
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
!	Dadd a DTensor to the end of the link	
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
!	Dadd a DTensor to the begin of the link	
!	inde is the index of the DTensor	
	subroutine Dpush_foI(h,T,inde)
		type(DTensorlink),intent(inout) ::h
		type(DTensor),intent(in):: T
		integer,intent(in)::inde(:)
		type(DTensornode),pointer ::node,p
		allocate(node)
		node%Ten=T
		node%indice=inde
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
!	add a Tensornode to the begin of the link	
!	there is no vector index in the  input Tensor
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
	subroutine Dpush_backnum(h,num)
		type(DTensorlink),intent(inout) ::h
		real*8,intent(in) ::num
		type(DTensor)::Ten
		Ten=(/num/)
		call Dpush_back(h,Ten)
		return
	end subroutine
!	add a compelx munber to the end of the link
!	there is vector index in the  input Tensor
	subroutine Dpush_backnumI(h,num,inde)
		type(DTensorlink),intent(inout) ::h
		real*8,intent(in) ::num
		integer,intent(in) :: inde(:)
		type(DTensor)::Ten
		Ten=(/num/)
		call Dpush_backI(h,Ten,inde)
		return
	end subroutine
	subroutine Dpush_fonum(h,num)
		type(DTensorlink),intent(inout) ::h
		real*8,intent(in) ::num
		type(DTensor)::Ten
		Ten=(/num/)
		call Dpush_fo(h,Ten)
		return
	end subroutine
!	add a compelx munber to the begin of the link	
!	there is vector index in the  input Tensor
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
			write(*,*)"ERROR in DTi"
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
		logical::conDTinu=.true.
		if(DlenOfIndice(h).ne.size(inde)) then
			write(*,*)"ERROR in DTiI"
		else
			p=>h%head
			i=1
			do while (conDTinu)
				if(inde.equ.p%indice) then
					DTiI=p%Ten
					conDTinu=.false.
				end if
				p=>p%next
				i=i+1
				if(i.eq.h%length+1) then
					conDTinu=.false.
				end if
			end do
		end if
		return
	end  function
!	return the inde Tensornode's address in the link
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
!	return the inde Tensornode's address in the link
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
					write(*,*)"ERROR in node_i_vec"
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
	
	logical function Difnextnode(p)
		type(DTensornode),pointer,intent(inout) ::p
		if(associated(p%next)) then
			Difnextnode=.true.
			p=>p%next
		else
			Difnextnode=.false.
		end if
	end function
!	p is a pointer of Tensornode,p points to next
	subroutine Dnextnode(p)
		type(DTensornode),pointer,intent(inout) ::p
		if(associated(p%next)) then
			p=>p%next
		else
			write(*,*)"Error in nextnode,p is pointing to the end of the link"
			stop
		end if
	end subroutine
!	p is a pointer of Tensornode,p points to head of the h	
	subroutine Dheadnode(h,p)
		type(DTensorlink),intent(in)::h
		type(DTensornode),pointer,intent(inout) ::p
		p=>h%head
	end subroutine
!	p is a pointer of Tensornode,p points to end of the h	
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
!check and output length of the link		
!check if the length of the link is equal to the Dvalue in head of the link
	integer function Dchecklength(link)
		type(DTensorlink),intent(in) :: link
		type(DTensornode),pointer ::p
		integer ::num
		Dchecklength=link%length
		p=>link%head
		num=0
		do while(associated(p))
			p=>p%next
			num=num+1
		end do
		if(num.ne.Dchecklength)then
			write(*,*)"The length of link is",num
			write(*,*)"The Dvalue of length in link is" ,Dchecklength
			write(*,*)"Error in length of link"
			write(*,*)"stop !"
			stop
		end if
		return
	end function
!	output length of the link		
	integer function Dlinklength(link)
		type(DTensorlink),intent(in) :: link
		Dlinklength=link%length
		return
	end function	
	
!	Dmodify the inde DElement of the link
!	The DTensor will be Dmodify
	subroutine Dmodify_link(h,T,inde)
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
			deallocate(p)
			p=>p1
		end do
		h%length=0
		nullify(h%head)
		nullify(h%Tend)
		return
	end subroutine
!	DTcopy the link hin to hout
!	hin will not change
	subroutine DTcopylink(hout,hin)
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
!	connect two links:
!			link1 :T1->T2->T4->T5...->TN
!			link2: TT1->TT2->..->TTM
!			result link will b2 T1->T2->..->TN->TT1->TT2->...->TTM
	type(DTensorlink) function Dconnectlink(link1,link2)
		type(DTensorlink),intent(in) :: link1,link2
		type(DTensorlink) :: Templink
		type(DTensornode),pointer ::p
		Dconnectlink=link1
		Templink=link2
		p=>Dconnectlink%Tend
		p%next=>Templink%head
		Dconnectlink%Tend=>Templink%Tend
		Dconnectlink%length=Dconnectlink%length+Templink%length
		return
	end function
!***************  Ddeletelink  **********************
!	link   :	T1->T2->...->T_{inde1-1}->T_inde1->T_{inde1+1}->....->T_{inde2-1}->T_inde2->T_{inde2+1}...->TN 
!  result:		T1->T2->...->T_{inde1-1}->T_{inde1+1}....->T_{inde2-1}->->T_{inde2+1}...->TN 
!  delete the inde1 to inde2 DTensor in the link,note that the deleted DTensor are include the DTensors of inde1 and inde2
	type(DTensorlink) function Ddeletelink(link_in,inde1,inde2) result(link)			
		type(DTensorlink),intent(in) :: link_in
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
		link=link_in
		if(inde1.eq.1) then!if[4]
			if(inde2.eq.link%length) then!if[5]
				call Dcleanlink(link)
				return
			end if!if[5]
			p2=>link%head
			do i=1,inde2
				p=>p2
				p2=>p2%next
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
			deallocate(p)
		end do
!		call DTDprint(p2%Ten)
		p1%next=>p2
		link%length=link%length-(inde2-inde1)-1
		return
	end function	
	
!Tensor=Tensorlink,if every element(Tensor) in the link is a one-element Tensor
	subroutine DTenLink(Ten,link)
		type(DTensor),intent(inout)::Ten
		type(DTensorlink),intent(in)::link
		real*8,allocatable::Tensordata(:)
		integer::Tenlen,i
		logical::goon
		type(DTensornode),pointer::p
		Tenlen=dlinklength(link)
		allocate(Tensordata(Tenlen))
		call Dheadnode(link,p)
		do i=1,Tenlen
			if(DgetTotaldata(p%Ten).eq.1)then
				Tensordata(i)=p%Ten.i.1
			else
				write(*,*)"error in assignment of a link to Tensor"
			end if
			goon=Difnextnode(p)
		end do
		Ten=Tensordata
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
			write(*,*)"There is no Tensor in the link"
			return
		end if
		write(*,*) "length of the link is",lenlink
		goon=DTenIndiceLog(link,1,indice)
		if(goon) then
			write(*,*)"indice of every Tensor in the link are"
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
					write(*,*)"error in the link,there are some Tensor with indice,some are not"
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











!******************  init_random_seed *****************
	subroutine Dinit_random_seed()
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
		maxV=maxDElement(T)
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
		maxV=maxDElement(T)
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
			if(isnan(val)) then
				DNANjudge=.true.
			return
			end if
		end do
		return
	end function			


	
end module
!***************************************************
!***************** END OF TNESOR *************
!***************************************************


			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
