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
!
!************************************************************
!************* START OF Tensor **********************
!************************************************************
!*****************     worning    *******************
!			1.The code is fortran90 version, gfortran4.8.4.
!			2.to compile the code one should link the files to lapack and blas.
!			3.Send email to sj.dong@outlook.com to report any bugs.
!
!
!
!

module Tensor_type
	use TData_module
	use Dimension_typede
	use Contract_real8_Tools
	use mpi
	use memory_type
	use Tools
	implicit none
	private
	character(len=1)::array_character_divider='|'

	public::Tensor
	type,extends (Tdata) :: Tensor
	contains
		procedure::outAllName,outAllName2
		generic,public::getAllName=>outAllName,outAllName2
		procedure::reset_dim_no_check1,reset_dim_no_check2
		generic,public::reset_dim_no_check => reset_dim_no_check1,reset_dim_no_check2
		!resetdim define in Dimension.f90
		procedure::splitDimension2,splitDimension3,splitDimensionAll!overwrite in dimension.f90
		procedure::fuseDimension_val,fuseDimension_vec!overwrite in dimension.f90

		procedure::allocateTensor1,allocateTensor2,allocateTensor3,allocateTensor4,&
				allocateTensor5,allocateTensor6,allocateTensor7,allocateTensor8,allocateTensor9
		generic,public::allocate=>allocateTensor1,allocateTensor2,allocateTensor3,allocateTensor4,&
				allocateTensor5,allocateTensor6,allocateTensor7,allocateTensor8,allocateTensor9
		procedure::TElement,TElement2
		procedure::iElement,iElement2,ielementAll
		procedure::sElement,sElement2,selementAll
		procedure::dElement,dElement2,delementAll
		procedure::cElement,cElement2,celementAll
		procedure::zElement,zElement2,zelementAll
		procedure::lElement,lElement2,lelementAll
		procedure::aElement,aElement2,aelementAll
		generic,public::i => TElement,TElement2
		generic,public::ii => iElement,iElement2,ielementAll
		generic,public::si => sElement,sElement2,selementAll
		generic,public::di => dElement,dElement2,delementAll
		generic,public::ci => cElement,cElement2,celementAll
		generic,public::zi => zElement,zElement2,zelementAll
		generic,public::li => lElement,lElement2,lelementAll
		generic,public::ai => aElement,aElement2,aelementAll

		procedure::TmaxElement,TmaxminElement,TminElement
		procedure::intmaxElement,intmaxminElement,intminElement
		procedure::realmaxElement,realmaxminElement,realminElement
		procedure::dblemaxElement,dblemaxminElement,dbleminElement
		procedure::cmplxmaxElement,cmplxmaxminElement,cmplxminElement
		procedure::dcmplxmaxElement,dcmplxmaxminElement,dcmplxminElement
		generic,public::max => TmaxElement,TmaxminElement
		generic,public::imax => intmaxElement,intmaxminElement
		generic,public::smax => realmaxElement,realmaxminElement
		generic,public::dmax => dblemaxElement,dblemaxminElement
		generic,public::cmax => cmplxmaxElement,cmplxmaxminElement
		generic,public::zmax => dcmplxmaxElement,dcmplxmaxminElement

		generic,public::min => TminElement,TmaxminElement
		generic,public::imin => intminElement,intmaxminElement
		generic,public::smin => realminElement,realmaxminElement
		generic,public::dmin => dbleminElement,dblemaxminElement
		generic,public::cmin => cmplxminElement,cmplxmaxminElement
		generic,public::zmin => dcmplxminElement,dcmplxmaxminElement
		!maxminflag=
		!'maxa': max abs 
		!'mina': min abs 
		!'maxr': max real
		!'minr': min real
		!'maxi': 0(not com) or max imag
		!'mini': 0(not com) or min imag


		procedure,public::maxmin => dcmplxmaxminElement
		procedure,public::imaxmin => intmaxminElement
		procedure,public::smaxmin => realmaxminElement
		procedure,public::dmaxmin => dblemaxminElement
		procedure,public::cmaxmin => cmplxmaxminElement
		procedure,public::zmaxmin => dcmplxmaxminElement
		
		
		procedure,public::sum =>sumTensorT
		procedure,public::isum =>sumTensori
		procedure,public::ssum =>sumTensors
		procedure,public::dsum =>sumTensord
		procedure,public::csum =>sumTensorc
		procedure,public::zsum =>sumTensorz

		procedure::Dprint2!overwrite in dimension.f90
		procedure::Dprint!overwrite in dimension.f90
		procedure::Tprint2,Tprint3,Tprint4
		procedure::TMprint2,TMprint3,TMprint4
		procedure::Tprint_file1,Tprint_file2
		procedure::TMprint_file1,TMprint_file2
		procedure::TDprint1,TDprint2,TDprint3
		procedure::Tread_data
		generic,public::readData =>Tread_data
		generic,public::info =>Tprint2,Tprint3,Tprint4 !Dprint2 define in diemsnion.f90
		generic,public::print =>TMprint2,TMprint3,TMprint4!Dprint  define in diemsnion.f90
		generic,public::diminfo=>TDprint1,TDprint2,TDprint3
		generic,public::writeinfo =>Tprint_file1,Tprint_file2
		generic,public::write =>Tprint_file1,Tprint_file2
		procedure::readdimension!overwrite in dimension.f90, read
		generic,public::printFile =>TMprint_file1,TMprint_file2
		procedure::subTen,subTen2,subTen2_Name,subTen3,subTen3_Name
		generic,public::subTensor=>subTen,subTen2,subTen2_Name,subTen3,subTen3_Name
		generic,public::subBlock=>subTen,subTen2,subTen2_Name,subTen3,subTen3_Name
		procedure,public::paste=>pasteTensorSubroutine
		!call T1%paste(T2,.true.)
		!			T1 :a [l,m,n,...] matrix
		!			T2 :a [l,m,n,...] matrix
		!			pasteTensor(T1,T2):a [2*l,m,n,...] matrix
		!        [1--->l, m, n,...] is T1
		!			[l+1--->2*l, m, n,...] is T2
		!        /   \
		!        | T1 |
		!        |----|
		!        | T2 |
		!        \    /
		!call T1%paste(T2,.false.)
		!			T1 :a [...,m,n,l] matrix
		!			T2 :a [...,m,n,l] matrix
		!			pasteTensor(T1,T2):a [...,m,n,2*l] matrix
		!        [... , m, n, 1--->l] is T1
		!			[... , m, n, l+1--->2*l] is T2
		!        /         \
		!        | T1 | T2 |
		!        \         /

		!Tdata(:)=a            setvalue(a)                -----set_value_Tensor_*
		!Tdata(i)=a            setvalue(i,a)              -----modifyTen_val_*
		!Tdata((/i,j,k/))=a    setvalue((/i,j,k/),a)      -----modifyTen_val1_*
		!Tdata(i)=T.i.1            setvalue(i,T)   T is a Tensor             -----modifyTen_val_*
		!Tdata((/i,j,k/))=T.i.1    setvalue((/i,j,k/),T)  T is a Tensor      -----modifyTen_val1_*
		!Tdata(:)=a(*)         setvalue(a) , a is a array -----store_*
		!Tdata(i1,i2)=a(*)     setvalue(i1,i2,a), a is a array -----storeSome_*
		!Tdata(i1,i2)=T        setvalue(i1,i2,T), T is a Tensor -----storeSome_Tensor*
		!Tdata(i1:i2)=a(i1':i2') a is a Tensor or array  --- modify_Some_Data1_*
		!Tdata(i1:i2,j1:j2)=a(i1':i2',j1':j2') a is a Tensor or array  --- modify_Some_Data2_*
		!To be finished
		!Tdata(i1:i2,j1:j2,k1:k2)=a(i1':i2',j1':j2',k1':k2')
		!Tdata(i1:i2,j1:j2,k1:k2,l1:l2)=a(i1':i2',j1':j2',k1':k2',l1':l2')

		procedure::modifyTen_val_int,modifyTen_val_real4,modifyTen_val_real8,modifyTen_val_com4,&
											modifyTen_val_com8,modifyTen_val_logi,modifyTen_val_char,&
											modifyTen_val1_int,modifyTen_val1_real4,modifyTen_val1_real8,modifyTen_val1_com4,&
											modifyTen_val1_com8,modifyTen_val1_logi,modifyTen_val1_char,modifyTen_val_Tensor,&
											modifyTen_val1_Tensor,&
											set_value_Tensor_int,set_value_Tensor_real4,set_value_Tensor_real8,set_value_Tensor_com4,&
											set_value_Tensor_com8,set_value_Tensor_logi,set_value_Tensor_char, store_T,&
											store_int,store_real4,store_real8,store_com4,store_com8,store_logi,store_char,&
											storeSome_int,storeSome_real4,storeSome_real8,storeSome_com4,storeSome_com8,storeSome_logi,storeSome_char,&
											storeSome_Tensor,&
											modify_Some_Data1_int,modify_Some_Data1_real4,modify_Some_Data1_real8,modify_Some_Data1_com4,&
											modify_Some_Data1_com8,modify_Some_Data1_char,modify_Some_Data1_logi,modify_Some_Data1_Tensor,&
											modify_Some_Data2_int,modify_Some_Data2_real4,modify_Some_Data2_real8,modify_Some_Data2_com4,&
											modify_Some_Data2_com8,modify_Some_Data2_char,modify_Some_Data2_logi,modify_Some_Data2_Tensor,&
											modify_Some_Data2_Tensor2
		generic,public::setValue =>modifyTen_val_int,modifyTen_val_real4,modifyTen_val_real8,modifyTen_val_com4,&
											modifyTen_val_com8,modifyTen_val_logi,modifyTen_val_char,&
											modifyTen_val1_int,modifyTen_val1_real4,modifyTen_val1_real8,modifyTen_val1_com4,&
											modifyTen_val1_com8,modifyTen_val1_logi,modifyTen_val1_char,modifyTen_val_Tensor,&
											modifyTen_val1_Tensor,&
											set_value_Tensor_int,set_value_Tensor_real4,set_value_Tensor_real8,set_value_Tensor_com4,&
											set_value_Tensor_com8,set_value_Tensor_logi,set_value_Tensor_char, store_T,&
											store_int,store_real4,store_real8,store_com4,store_com8,store_logi,store_char,&
											storeSome_int,storeSome_real4,storeSome_real8,storeSome_com4,storeSome_com8,storeSome_logi,storeSome_char,&
											storeSome_Tensor,&
											modify_Some_Data1_int,modify_Some_Data1_real4,modify_Some_Data1_real8,modify_Some_Data1_com4,&
											modify_Some_Data1_com8,modify_Some_Data1_char,modify_Some_Data1_logi,modify_Some_Data1_Tensor,&
											modify_Some_Data2_int,modify_Some_Data2_real4,modify_Some_Data2_real8,modify_Some_Data2_com4,&
											modify_Some_Data2_com8,modify_Some_Data2_char,modify_Some_Data2_logi,modify_Some_Data2_Tensor,&
											modify_Some_Data2_Tensor2
		procedure,public::zero=>set_zero_Tensor
		procedure::eye_Tensor2,eye_Tensor3,eye_Tensor4,eye_Tensor5,eye_Tensor6
		generic,public::eye=> eye_Tensor2,eye_Tensor3,eye_Tensor4,eye_Tensor5,eye_Tensor6
		procedure::iwhichindex,swhichindex,dwhichindex,cwhichindex,zwhichindex,awhichindex,awhichindex2
		procedure::sortTensor1,sortTensor2,sortTensor3,sortTensor4,sortTensor5,sortTensor6
		generic,public::which=>iwhichindex,swhichindex,dwhichindex,cwhichindex,zwhichindex,awhichindex,awhichindex2
		generic,public::sort=>sortTensor1,sortTensor2,sortTensor3,sortTensor4,sortTensor5,sortTensor6

		procedure::enlargeTensorReal8,enlargeTensorReal4,enlargeTensorInt
		procedure::enlargeTensorReal8_name,enlargeTensorReal4_name,enlargeTensorInt_name
		procedure::enlargeTensorReal8Routine,enlargeTensorReal4Routine,enlargeTensorIntRoutine
		procedure::enlargeTensorReal8_nameRoutine,enlargeTensorReal4_nameRoutine,enlargeTensorInt_nameRoutine
		procedure::enlargeTensorAllReal8Routine,enlargeTensorAllReal4Routine,enlargeTensorAllIntRoutine
		procedure::enlargeTensorAllReal8,enlargeTensorAllReal4,enlargeTensorAllInt
		generic,public::enLargeTensor=>enlargeTensorReal8,enlargeTensorReal4,enlargeTensorInt,enlargeTensorAllReal8&
							,enlargeTensorAllReal4,enlargeTensorAllInt,&
							enlargeTensorReal8_name,enlargeTensorReal4_name,enlargeTensorInt_name
		generic,public::GetenLarge=>enlargeTensorReal8,enlargeTensorReal4,enlargeTensorInt,enlargeTensorAllReal8&
							,enlargeTensorAllReal4,enlargeTensorAllInt,&
							enlargeTensorReal8_name,enlargeTensorReal4_name,enlargeTensorInt_name
		generic,public::enLarge=>enlargeTensorReal8Routine,enlargeTensorReal4Routine,enlargeTensorIntRoutine,&
								enlargeTensorAllReal8Routine,enlargeTensorAllReal4Routine,enlargeTensorAllIntRoutine,&
								enlargeTensorReal8_nameRoutine,enlargeTensorReal4_nameRoutine,enlargeTensorInt_nameRoutine

		procedure::set_random_Tensor,set_random_Tensor_region4,set_random_Tensor_regioni
		generic,public::random =>set_random_Tensor,set_random_Tensor_region4,set_random_Tensor_regioni

		procedure::permutation_routine,permutation_Name_routine,permutation_Tensor_routine
		procedure::permutefo_routine,permutefo_name_routine,permutefo_vec_routine,permutefo_vec_Name_routine&
										,permutefo_Tensor_routine
		procedure::permuteback_routine,permuteback_name_routine,permuteback_vec_routine&
										,permuteback_vec_name_routine,permuteback_Tensor_routine
		generic,public::permute=>permutation_routine,permutation_Name_routine,permutation_Tensor_routine
		generic,public::forward=>permutefo_routine,permutefo_name_routine,permutefo_vec_routine,permutefo_vec_Name_routine&
										,permutefo_Tensor_routine
		generic,public::backward=>permuteback_routine,permuteback_name_routine,permuteback_vec_routine&
										,permuteback_vec_name_routine,permuteback_Tensor_routine

		procedure,public::norm => TnormTensor
		procedure,public::inorm => inormTensor
		procedure,public::snorm => snormTensor
		procedure,public::dnorm => dnormTensor
		procedure,public::cnorm => cnormTensor
		procedure,public::znorm => znormTensor
		procedure,public::norm2 => Tnorm2Tensor
		procedure,public::inorm2 => inorm2Tensor
		procedure,public::snorm2 => snorm2Tensor
		procedure,public::dnorm2 => dnorm2Tensor
		procedure,public::cnorm2 => cnorm2Tensor
		procedure,public::znorm2 => znorm2Tensor
		
		procedure,public::trace  => TtraceTensor
		procedure,public::itrace => itraceTensor
		procedure,public::strace => straceTensor
		procedure,public::dtrace => dtraceTensor
		procedure,public::ctrace => ctraceTensor
		procedure,public::ztrace => ztraceTensor

		procedure,public::eigvalue!check,only ok on sysmetry martix
		procedure,public::eig =>eigTensor!check,only ok on sysmetry martix
		generic,public::inv=>inverse,inverseTen
		generic,public::invTensor=>inverse,inverseTen
		procedure,public::linequ
		!linear equations
		!linequ:X=A%linequ(B),solve A*X=B,X could be a matrix

		procedure::SVDcutoff,SVDcutoff_name,SVDNameLeftRoutine,SVDNameLeftRoutine_no_cut,SVDTensor_Leg2_Routine
		procedure::SVDcutoff_kill_inData,SVDcutoff_name_kill_inData,SVDNameLeftRoutine_kill_indata,&
					SVDNameLeftRoutine_kill_indata_no_cut,SVDTensor_Leg2_Routine_kill_inData
		procedure::SVDTensor_noname,SVDTensor_name,SVDNameLeft,SVDNameLeft_no_cut,SVDTensor_Leg2
		generic,public::SVD=>SVDcutoff,SVDcutoff_name,SVDNameLeftRoutine,SVDNameLeftRoutine_no_cut,SVDTensor_Leg2_Routine
		generic,public::SVDroutine=>SVDcutoff,SVDcutoff_name,SVDNameLeftRoutine,SVDNameLeftRoutine_no_cut,SVDTensor_Leg2_Routine
		generic,public::SVDKill=>SVDcutoff_kill_inData,SVDcutoff_name_kill_inData,SVDNameLeftRoutine_kill_indata,&
					SVDNameLeftRoutine_kill_indata_no_cut,SVDTensor_Leg2_Routine_kill_inData
		generic,public::getSVD=>SVDTensor_noname,SVDTensor_name,SVDNameLeft,SVDNameLeft_no_cut,SVDTensor_Leg2
		generic,public::SVDTensor=>SVDTensor_noname,SVDTensor_name,SVDNameLeft,SVDNameLeft_no_cut,SVDTensor_Leg2

		procedure::LQdecomposition1,LQTensorNameLeftRoutine,LQTensorLegNameRoutine
		procedure::LQdecomposition_kill_inData,LQTensorNameLeftRoutine_kill_inData,LQTensorLegNameRoutine_kill_inData
		procedure::LQTensor_noName,LQTensor_name,LQTensorNameLeft,LQTensorLegName
		generic,public::LQ =>LQdecomposition1,LQTensorNameLeftRoutine,LQTensorLegNameRoutine
		generic,public::LQKill =>LQdecomposition_kill_inData,LQTensorNameLeftRoutine_kill_inData,LQTensorLegNameRoutine_kill_inData
		generic,public::LQTensor=>LQTensor_noName,LQTensor_name,LQTensorNameLeft,LQTensorLegName
		generic,public::getLQ=>LQTensor_noName,LQTensor_name,LQTensorNameLeft,LQTensorLegName

		procedure::QRdecomposition1,QRTensorNameLeftRoutine,QRTensorLegNameRoutine
		procedure::QRdecomposition_kill_inData,QRTensorNameLeftRoutine_kill_inData,QRTensorLegNameRoutine_kill_inData
		procedure::QRTensor_noName,QRTensor_name,QRTensorNameLeft,QRTensorLegName
		generic,public::QR =>QRdecomposition1,QRTensorNameLeftRoutine,QRTensorLegNameRoutine
		generic,public::QRKill =>QRdecomposition_kill_inData,QRTensorNameLeftRoutine_kill_inData,QRTensorLegNameRoutine_kill_inData
		generic,public::QRTensor=>QRTensor_noName,QRTensor_name,QRTensorNameLeft,QRTensorLegName
		generic,public::getQR=>QRTensor_noName,QRTensor_name,QRTensorNameLeft,QRTensorLegName

		procedure::ProductTensorRoutine1,ProductTensorRoutine2,ProductTensorRoutine3
		generic,public::ProductTensorRoutine=>ProductTensorRoutine1,ProductTensorRoutine2,ProductTensorRoutine3

		procedure::contract_name_routine,contract_name_routine1,contract_name_routine2
		procedure::contract_name_routine4,contract_name_routine5,contract_name_routine6
		procedure::contract_name_ownlegs_routine,contract_ownlegs_routine
		generic,public::contract=>contract_name_routine,contract_name_routine1,contract_name_routine2,&
											contract_name_routine4,contract_name_routine5,contract_name_routine6,&
											contract_name_ownlegs_routine,contract_ownlegs_routine

		!********* assignment **************************
 		procedure::DimInitialization!overwrite in dimension.f90
 		procedure::DimInitialization2!overwrite in dimension.f90
 		procedure,pass(Dimen)::copyDimToVec!overwrite in dimension.f90
 		procedure::assignmentTenarray_int2
 		procedure::assignmentTenarray_int3
 		procedure::assignmentTenarray_int4
 		procedure::assignmentTenarray_real4_1
 		procedure::assignmentTenarray_real4_2
 		procedure::assignmentTenarray_real4_3
 		procedure::assignmentTenarray_real4_4
 		procedure::assignmentTenarray_real8_1
 		procedure::assignmentTenarray_real8_2
 		procedure::assignmentTenarray_real8_3
 		procedure::assignmentTenarray_real8_4
 		procedure::assignmentTenarray_com4_1
 		procedure::assignmentTenarray_com4_2
 		procedure::assignmentTenarray_com4_3
 		procedure::assignmentTenarray_com4_4
 		procedure::assignmentTenarray_com8_1
 		procedure::assignmentTenarray_com8_2
 		procedure::assignmentTenarray_com8_3
 		procedure::assignmentTenarray_com8_4
 		procedure::assignmentTenarray_logi1
 		procedure::assignmentTenarray_logi2
 		procedure::assignmentTenarray_logi3
 		procedure::assignmentTenarray_logi4
 		procedure::assignmentTenarray_char1
 		procedure::assignmentTenarray_char2
 		procedure::assignmentTenarray_char3
 		procedure::assignmentTenarray_char4
 		procedure::assignmentTenNum_int
 		procedure::assignmentTenNum_real4
 		procedure::assignmentTenNum_real8
 		procedure::assignmentTenNum_com4
 		procedure::assignmentTenNum_com8
 		procedure::assignmentTenNum_logi
 		procedure::assignmentTenNum_char
 		procedure,pass(T)::assignmentNumTen_int
 		procedure,pass(T)::assignmentNumTen_real4
 		procedure,pass(T)::assignmentNumTen_real8
 		procedure,pass(T)::assignmentNumTen_com4
 		procedure,pass(T)::assignmentNumTen_com8
 		procedure,pass(T)::assignmentNumTen_logi
 		procedure,pass(T)::assignmentNumTen_char
		procedure,pass(T)::copyTensor_int_dim2
		procedure,pass(T)::copyTensor_int_dim3
		procedure,pass(T)::copyTensor_int_dim4
		procedure,pass(T)::copyTensor_real4_dim1
		procedure,pass(T)::copyTensor_real4_dim2
		procedure,pass(T)::copyTensor_real4_dim3
		procedure,pass(T)::copyTensor_real4_dim4
		procedure,pass(T)::copyTensor_real8_dim1
		procedure,pass(T)::copyTensor_real8_dim2
		procedure,pass(T)::copyTensor_real8_dim3
		procedure,pass(T)::copyTensor_real8_dim4
		procedure,pass(T)::copyTensor_com4_dim1
		procedure,pass(T)::copyTensor_com4_dim2
		procedure,pass(T)::copyTensor_com4_dim3
		procedure,pass(T)::copyTensor_com4_dim4
		procedure,pass(T)::copyTensor_com8_dim1
		procedure,pass(T)::copyTensor_com8_dim2
		procedure,pass(T)::copyTensor_com8_dim3
		procedure,pass(T)::copyTensor_com8_dim4
		procedure,pass(T)::copyTensor_logi_dim1
		procedure,pass(T)::copyTensor_logi_dim2
		procedure,pass(T)::copyTensor_logi_dim3
		procedure,pass(T)::copyTensor_logi_dim4
		procedure,pass(T)::copyTensor_char_dim1
		procedure,pass(T)::copyTensor_char_dim2
		procedure,pass(T)::copyTensor_char_dim3
		procedure,pass(T)::copyTensor_char_dim4
 		generic,public :: assignment(=) =>  assignmentTenarray_int2,&
									 		assignmentTenarray_int3,&
									 		assignmentTenarray_int4,&
									 		assignmentTenarray_real4_1,&
									 		assignmentTenarray_real4_2,&
									 		assignmentTenarray_real4_3,&
									 		assignmentTenarray_real4_4,&
									 		assignmentTenarray_real8_1,&
									 		assignmentTenarray_real8_2,&
									 		assignmentTenarray_real8_3,&
									 		assignmentTenarray_real8_4,&
									 		assignmentTenarray_com4_1,&
									 		assignmentTenarray_com4_2,&
									 		assignmentTenarray_com4_3,&
									 		assignmentTenarray_com4_4,&
									 		assignmentTenarray_com8_1,&
									 		assignmentTenarray_com8_2,&
									 		assignmentTenarray_com8_3,&
									 		assignmentTenarray_com8_4,&
									 		assignmentTenarray_logi1,&
									 		assignmentTenarray_logi2,&
									 		assignmentTenarray_logi3,&
									 		assignmentTenarray_logi4,&
									 		assignmentTenarray_char1,&
									 		assignmentTenarray_char2,&
									 		assignmentTenarray_char3,&
									 		assignmentTenarray_char4,&
									 		assignmentTenNum_int,&
									 		assignmentTenNum_real4,&
									 		assignmentTenNum_real8,&
									 		assignmentTenNum_com4,&
									 		assignmentTenNum_com8,&
									 		assignmentTenNum_logi,&
									 		assignmentTenNum_char,&
									 		assignmentNumTen_int,&
									 		assignmentNumTen_real4,&
									 		assignmentNumTen_real8,&
									 		assignmentNumTen_com4,&
									 		assignmentNumTen_com8,&
									 		assignmentNumTen_logi,&
									 		assignmentNumTen_char,&
									 		copyTensor_int_dim2,&
											copyTensor_int_dim3,&
											copyTensor_int_dim4,&
											copyTensor_real4_dim1,&
											copyTensor_real4_dim2,&
											copyTensor_real4_dim3,&
											copyTensor_real4_dim4,&
											copyTensor_real8_dim1,&
											copyTensor_real8_dim2,&
											copyTensor_real8_dim3,&
											copyTensor_real8_dim4,&
											copyTensor_com4_dim1,&
											copyTensor_com4_dim2,&
											copyTensor_com4_dim3,&
											copyTensor_com4_dim4,&
											copyTensor_com8_dim1,&
											copyTensor_com8_dim2,&
											copyTensor_com8_dim3,&
											copyTensor_com8_dim4,&
											copyTensor_logi_dim1,&
											copyTensor_logi_dim2,&
											copyTensor_logi_dim3,&
											copyTensor_logi_dim4,&
											copyTensor_char_dim1,&
											copyTensor_char_dim2,&
											copyTensor_char_dim3,&
											copyTensor_char_dim4
		procedure::AddFunc!over write in Dimension.f90
		procedure::AddFunc2!over write in Dimension.f90
		procedure,pass(Dimen)::AddFunc3!over write in Dimension.f90
		procedure::add_int
		procedure::add_real4
		procedure::add_real8
		procedure::add_com4
		procedure::add_com8
		procedure::add_char
		procedure,pass(T1)::add_int_
		procedure,pass(T1)::add_real4_
		procedure,pass(T1)::add_real8_
		procedure,pass(T1)::add_com4_
		procedure,pass(T1)::add_com8_
		procedure,pass(T1)::add_char_
 		generic,public :: operator(+) => add_int,&
										add_real4,&
										add_real8,&
										add_com4,&
										add_com8,&
										add_char,&
										add_int_,&
										add_real4_,&
										add_real8_,&
										add_com4_,&
										add_com8_,&
										add_char_
		procedure::minus
		procedure::minus_int
		procedure::minus_real4
		procedure::minus_real8
		procedure::minus_com4
		procedure::minus_com8
		procedure,pass(T1)::minus_int_
		procedure,pass(T1)::minus_real4_
		procedure,pass(T1)::minus_real8_
		procedure,pass(T1)::minus_com4_
		procedure,pass(T1)::minus_com8_
		generic,public :: operator(-) => minus,&
										minus_int,&
										minus_real4,&
										minus_real8,&
										minus_com4,&
										minus_com8,&
										minus_int_,&
										minus_real4_,&
										minus_real8_,&
										minus_com4_,&
										minus_com8_
		procedure::multiply_number_int
		procedure,pass(T1)::multiply_number_int_
		procedure::multiply_number_real4
		procedure,pass(T1)::multiply_number_real4_
		procedure::multiply_number_real8
		procedure,pass(T1)::multiply_number_real8_
		procedure::multiply_number_com4
		procedure,pass(T1)::multiply_number_com4_
		procedure::multiply_number_com8
		procedure,pass(T1)::multiply_number_com8_
		procedure::ProductTensor
		generic,public :: operator(*) =>multiply_number_int,&
										multiply_number_int_,&
										multiply_number_real4,&
										multiply_number_real4_,&
										multiply_number_real8,&
										multiply_number_real8_,&
										multiply_number_com4,&
										multiply_number_com4_,&
										multiply_number_com8,&
										multiply_number_com8_,&
										ProductTensor
		procedure::divide_Tensor
		procedure,pass(T)::int_divide_Tensor
		procedure,pass(T)::real4_divide_Tensor
		procedure,pass(T)::real8_divide_Tensor
		procedure,pass(T)::com4_divide_Tensor
		procedure,pass(T)::com8_divide_Tensor
		procedure::divide_num_int
		procedure::divide_num_real4
		procedure::divide_num_real8
		procedure::divide_num_com4
		procedure::divide_num_com8
		generic,public :: operator(/) =>divide_Tensor,&
										int_divide_Tensor,&
										real4_divide_Tensor,&
										real8_divide_Tensor,&
										com4_divide_Tensor,&
										com8_divide_Tensor,&
										divide_num_int,&
										divide_num_real4,&
										divide_num_real8,&
										divide_num_com4,&
										divide_num_com8
		procedure::equal_of_dim!overwrite in Dimension.f90
		procedure::T_eq_int
		procedure::T_eq_real4
		procedure::T_eq_real8
		procedure,pass(T)::int_eq_T
		procedure,pass(T)::real4_eq_T
		procedure,pass(T)::real8_eq_T
		generic,public :: operator(.equ.) =>T_eq_int,&
											T_eq_real4,&
											T_eq_real8,&
											int_eq_T,&
											real4_eq_T,&
											real8_eq_T
		generic,public :: operator(.eq.) =>T_eq_int,&
											T_eq_real4,&
											T_eq_real8,&
											int_eq_T,&
											real4_eq_T,&
											real8_eq_T
		procedure::gt_of_Tensor
		procedure::T_gt_int
		procedure::T_gt_real4
		procedure::T_gt_real8
		procedure,pass(T)::int_gt_T
		procedure,pass(T)::real4_gt_T
		procedure,pass(T)::real8_gt_T
		generic,public :: operator(.gt.) =>gt_of_Tensor,&
											T_gt_int,&
											T_gt_real4,&
											T_gt_real8,&
											int_gt_T,&
											real4_gt_T,&
											real8_gt_T
		procedure::ge_of_Tensor
		procedure::T_ge_int
		procedure::T_ge_real4
		procedure::T_ge_real8
		procedure,pass(T)::int_ge_T
		procedure,pass(T)::real4_ge_T
		procedure,pass(T)::real8_ge_T	
		generic,public :: operator(.ge.) =>ge_of_Tensor,&
											T_ge_int,&
											T_ge_real4,&
											T_ge_real8,&
											int_ge_T,&
											real4_ge_T,&
											real8_ge_T	
		procedure::lt_of_Tensor
		procedure::T_lt_int
		procedure::T_lt_real4
		procedure::T_lt_real8
		procedure,pass(T)::int_lt_T
		procedure,pass(T)::real4_lt_T
		procedure,pass(T)::real8_lt_T
		generic,public :: operator(.lt.) =>lt_of_Tensor,&
											T_lt_int,&
											T_lt_real4,&
											T_lt_real8,&
											int_lt_T,&
											real4_lt_T,&
											real8_lt_T
		procedure::le_of_Tensor
		procedure::T_le_int
		procedure::T_le_real4
		procedure::T_le_real8
		procedure,pass(T)::int_le_T
		procedure,pass(T)::real4_le_T
		procedure,pass(T)::real8_le_T
		generic,public :: operator(.le.) =>le_of_Tensor,&
											T_le_int,&
											T_le_real4,&
											T_le_real8,&
											int_le_T,&
											real4_le_T,&
											real8_le_T
		procedure::permute_rank2 !permute the tensor whose rank is 2,if rank=1 ,do nothing
		procedure::permute_rank3 !permute the tensor whose rank is 3
		procedure::permutation	 !permute the tensor of any rank,give the new order of the dimension
								 !		If operate on a big Tensor,use permute_rank2 or permute_rank3,
								 !	 they are faster.if rank>3,use contract to reshape
		procedure::permutation_name!input a character(:) as the new order
		procedure::permutation_Tensor ! if input 'A', goes worg, to be modify
		generic,public :: operator(.p.) =>permute_rank2,&
											permute_rank3,&
											permutation,&
											permutation_name,&
											permutation_Tensor
		procedure::permutefo!permute the inde index to the first
										  !T_{1,2,3,..,i,..,n},permutefo(T,i)=_{i,1,2,3,..,i-1,i+1,..,n}
										  !T_{1,2,3,..,j,..,i,.,k,...,n},permutefo_vec(T,(/i,j,k/))=_{i,j,k,1,2,3,...,n}
		procedure::permutefo_vec
		procedure::permutefo_name
		procedure::permutefo_vec_name
		procedure::permutefo_Tensor
		generic,public :: operator(.pf.) =>permutefo,&
											permutefo_vec,&
											permutefo_name,&
											permutefo_vec_name,&
											permutefo_Tensor
		procedure::permuteback
		procedure::permuteback_vec
		procedure::permuteback_name
		procedure::permuteback_vec_name
		procedure::permuteback_Tensor
		generic,public :: operator(.pb.) =>permuteback,&
											permuteback_vec,&
											permuteback_name,&
											permuteback_vec_name,&
											permuteback_Tensor
		procedure::permuteInde,permuteInde_name
		generic,public :: operator(.pi.) =>permuteInde,permuteInde_name
		procedure::permutebackInde,permutebackInde_name
		generic,public :: operator(.pbi.) =>permutebackInde,permutebackInde_name
		procedure::combinationCol,combinationRow
		generic,public :: operator(.AddCol.)=>combinationCol
		generic,public :: operator(.AddRow.)=>combinationRow
		!combinationCol:
		!			T1 :a [...,l,m,n] matrix
		!			T2 :a [...,l,m,n] matrix
		!			combination(T1,T2):a [...,l,m,n,2] matrix
		!			or 
		!			T1 :a [...,m,n,l] matrix
		!			T2 :a [...,m,n] matrix
		!			combination(T1,T2):a [...,m,n,l+1] matrix
		!combinationrow
		!			T1 :a [l,m,n,...] matrix
		!			T2 :a [l,m,n,...] matrix
		!			combinationrow(T1,T2):a [2,l,m,n,...] matrix
		!			or 
		!			T1 :a [l,m,n,...] matrix
		!			T2 :a [m,n,...] matrix
		!			combinationrow(T1,T2):a [l+1,m,n,...] matrix	

		procedure::pasteTensorRow,pasteTensorCol
		generic,public :: operator(.RPaste.)=>pasteTensorRow
		generic,public :: operator(.CPaste.)=>pasteTensorCol
		!			T=T1.Rpaste.T2
		!			T1 :a [l,m,n,...] matrix
		!			T2 :a [l,m,n,...] matrix
		!			pasteTensorRow(T1,T2):a [2*l,m,n,...] matrix
		!        [1--->l, m, n,...] is T1
		!			[l+1--->2*l, m, n,...] is T2
		!        /   \
		!        | T1 |
		!        |----|
		!        | T2 |
		!        \    /
		!			T=T1.Cpaste.T2
		!pasteTensorCol(T1,T2,.false.)
		!			T1 :a [...,m,n,l] matrix
		!			T2 :a [...,m,n,l] matrix
		!			pasteTensor(T1,T2):a [...,m,n,2*l] matrix
		!        [... , m, n, 1--->l] is T1
		!			[... , m, n, l+1--->2*l] is T2
		!        /         \
		!        | T1 | T2 |
		!        \         /
		generic,public :: operator(.i.) => TElement,TElement2
		generic,public :: operator(.ii.) => iElement,iElement2,ielementAll
		generic,public :: operator(.si.) => sElement,sElement2,selementAll
		generic,public :: operator(.di.) => dElement,dElement2,delementAll
		generic,public :: operator(.ci.) => cElement,cElement2,celementAll
		generic,public :: operator(.zi.) => zElement,zElement2,zelementAll
		generic,public :: operator(.li.) => lElement,lElement2,lelementAll
		generic,public :: operator(.ai.) => aElement,aElement2,aelementAll
		procedure::conjugate,Htranspose2,Htranspose,transposeTensor
		generic,public :: operator(.con.) => conjugate
		generic,public :: operator(.H.) => Htranspose
		generic,public :: operator(.Hn.) => Htranspose2
		generic,public :: operator(.T.) => transposeTensor
		procedure::directProductM,directProduct,directProductTensor
		generic,public :: operator(.mxx.) => directProductM!direct Product,input two rank<=2 Tensor only
		generic,public :: operator(.xx.) => directProduct!direct Product,input two rank<=2 Tensor only
		generic,public :: operator(.kron.) => directProductTensor!direct Product,support any rank tensor and keep their TensorName
		procedure::TdotTensor,idotTensor,sdotTensor,ddotTensor,cdotTensor,zdotTensor
		! dot product conjugating the first vector,The Tensor will be regard as a vector

		generic,public :: operator(.x.) => TdotTensor
		generic,public :: operator(.ix.) => idotTensor
		generic,public :: operator(.sx.) => sdotTensor
		generic,public :: operator(.dx.) => ddotTensor
		generic,public :: operator(.cx.) => cdotTensor
		generic,public :: operator(.zx.) => zdotTensor
		procedure::TdotUTensor,idotUTensor,sdotUTensor,ddotUTensor,cdotUTensor,zdotUTensor
		! dot product DO NOT conjugating the first vector,The Tensor will be regard as a vector

		generic,public :: operator(.dot.) => TdotUTensor
		generic,public :: operator(.idot.) => idotUTensor
		generic,public :: operator(.sdot.) => sdotUTensor
		generic,public :: operator(.ddot.) => ddotUTensor
		generic,public :: operator(.cdot.) => cdotUTensor
		generic,public :: operator(.zdot.) => zdotUTensor
		procedure::inverse,inverseTen
		generic,public :: operator(.inv.) =>inverse,inverseTen
	end type

	interface Tensor
		procedure constructor_int_scal
		procedure constructor_int_dim1
		procedure constructor_int_dim2
		procedure constructor_int_dim3
		procedure constructor_int_dim4
		
		procedure constructor_real4_scal
		procedure constructor_real4_dim1
		procedure constructor_real4_dim2
		procedure constructor_real4_dim3
		procedure constructor_real4_dim4
		
		procedure constructor_real8_scal
		procedure constructor_real8_dim1
		procedure constructor_real8_dim2
		procedure constructor_real8_dim3
		procedure constructor_real8_dim4
		
		procedure constructor_com4_scal
		procedure constructor_com4_dim1
		procedure constructor_com4_dim2
		procedure constructor_com4_dim3
		procedure constructor_com4_dim4
		
		procedure constructor_com8_scal
		procedure constructor_com8_dim1
		procedure constructor_com8_dim2
		procedure constructor_com8_dim3
		procedure constructor_com8_dim4
		
		procedure constructor_logi_scal
		procedure constructor_logi_dim1
		procedure constructor_logi_dim2
		procedure constructor_logi_dim3
		procedure constructor_logi_dim4
		
		procedure constructor_char_scal
		procedure constructor_char_dim1
		procedure constructor_char_dim2
		procedure constructor_char_dim3
		procedure constructor_char_dim4
	end interface

	public::writemess
	interface writemess
		module procedure writemess_Tensor
		module procedure writemess_Tensor_form
	end interface

	public::eye
	interface eye! some subroutine will conflict 
		module procedure  eye_int
		module procedure  eye_real4
		module procedure  eye_real8
		module procedure  eye_com4
		module procedure  eye_com8
		module procedure  eye_logi
		module procedure  eye_char
		module procedure  eye_Tensor
		module procedure  diag_int
		module procedure  identity_matrix_Tensor
		module procedure  diag_real4
		module procedure  diag_real8
		module procedure  diag_com4
		module procedure  diag_com8
		module procedure  diag_logi
		module procedure  diag_char
		module procedure  diag_type
	end interface

	public::pauli_matrix
	!		pauli_matrix(Tx,Ty,Tz,num)
	public::expm
	interface expm
		module procedure  expmTensor
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
		module procedure contract_name_ownlegs
		module procedure contract_ownlegs
	end interface

	interface modifyTen_val_class
		module procedure modifyTen_val_class_i
		module procedure modifyTen_val_class_s
		module procedure modifyTen_val_class_d
		module procedure modifyTen_val_class_c
		module procedure modifyTen_val_class_z
		module procedure modifyTen_val_class_l
		module procedure modifyTen_val_class_a
	end interface

	public::assignment(=)
	interface assignment(=)
		module procedure assignmentTenArray
	end interface

	public::set_xgesdd_subroutine,set_xgesvd_subroutine
	public::set_array_character_divider,get_array_character_divider
	
	!**************************************************************
	!
	!           WorkingMemory
	!
	!**************************************************************
	type(memory),save::WorkingMemory
	type(Tensor),target,save::WorkingTensor1,WorkingTensor2!use in contract
	type(dimension),target,save::Workingdimension4,Workingdimension5!!use in contract
	type(dimension),target,save::Workingdimension1,Workingdimension2,Workingdimension3!use in (*)

contains

	function outAllName(T,w)
		class(Tensor),allocatable::outAllName
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::w
		logical::goon
		character(len=max_len_of_char_in_TData),allocatable::outchar(:)
		allocate(Tensor::outAllName)
		goon=T%Dimension%outSomedimensionName(outchar,w)
		if(goon)then
			outAllName=outchar
		else
			call outAllName%empty()
		end if
		return
	end function
	function outAllName2(T)
		class(Tensor),allocatable::outAllName2
		class(Tensor),intent(in)::T
		allocate(Tensor::outAllName2)
		outAllName2=T%Dimension%getName()
		return
	end function


	!********************************************************************************************
	!********************************************************************************************
	!                                Tensor's  constructor    
	!
	!*********************************************************************************************
	!*********************************************************************************************

	function constructor_int_scal(val,dimen)result(Res)
		integer,intent(in)::val
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=val
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_int_dim1(vec,dimen)result(Res)
		integer,intent(in)::vec(:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_int_dim2(vec,dimen)result(Res)
		integer,intent(in)::vec(:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_int_dim3(vec,dimen)result(Res)
		integer,intent(in)::vec(:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_int_dim4(vec,dimen)result(Res)
		integer,intent(in)::vec(:,:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_real4_scal(val,dimen)result(Res)
		real(kind=4),intent(in)::val
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=val
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_real4_dim1(vec,dimen)result(Res)
		real(kind=4),intent(in)::vec(:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_real4_dim2(vec,dimen)result(Res)
		real(kind=4),intent(in)::vec(:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_real4_dim3(vec,dimen)result(Res)
		real(kind=4),intent(in)::vec(:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_real4_dim4(vec,dimen)result(Res)
		real(kind=4),intent(in)::vec(:,:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_real8_scal(val,dimen)result(Res)
		real(kind=8),intent(in)::val
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=val
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_real8_dim1(vec,dimen)result(Res)
		real(kind=8),intent(in)::vec(:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_real8_dim2(vec,dimen)result(Res)
		real(kind=8),intent(in)::vec(:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_real8_dim3(vec,dimen)result(Res)
		real(kind=8),intent(in)::vec(:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_real8_dim4(vec,dimen)result(Res)
		real(kind=8),intent(in)::vec(:,:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_com4_scal(val,dimen)result(Res)
		complex(kind=4),intent(in)::val
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=val
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_com4_dim1(vec,dimen)result(Res)
		complex(kind=4),intent(in)::vec(:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_com4_dim2(vec,dimen)result(Res)
		complex(kind=4),intent(in)::vec(:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_com4_dim3(vec,dimen)result(Res)
		complex(kind=4),intent(in)::vec(:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_com4_dim4(vec,dimen)result(Res)
		complex(kind=4),intent(in)::vec(:,:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_com8_scal(val,dimen)result(Res)
		complex(kind=8),intent(in)::val
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=val
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_com8_dim1(vec,dimen)result(Res)
		complex(kind=8),intent(in)::vec(:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_com8_dim2(vec,dimen)result(Res)
		complex(kind=8),intent(in)::vec(:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_com8_dim3(vec,dimen)result(Res)
		complex(kind=8),intent(in)::vec(:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_com8_dim4(vec,dimen)result(Res)
		complex(kind=8),intent(in)::vec(:,:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_logi_scal(val,dimen)result(Res)
		logical,intent(in)::val
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=val
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_logi_dim1(vec,dimen)result(Res)
		logical,intent(in)::vec(:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_logi_dim2(vec,dimen)result(Res)
		logical,intent(in)::vec(:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_logi_dim3(vec,dimen)result(Res)
		logical,intent(in)::vec(:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_logi_dim4(vec,dimen)result(Res)
		logical,intent(in)::vec(:,:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_char_scal0(val,dimen)result(Res)
		character(len=*),intent(in)::val
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=val
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	subroutine set_array_character_divider(cha)
		character(len=1),intent(in)::cha
		array_character_divider=cha
		return
	end subroutine
	function get_array_character_divider()
		character(len=1)::get_array_character_divider
		get_array_character_divider=array_character_divider
		return
	end function
	function constructor_char_scal(cha_,dimen)result(res)
		character(len=*),intent(in)::cha_
		integer,optional,intent(in)::dimen(:)
		character(len=len(trim(adjustl(cha_))))::cha
		integer::TotalData,lenCha,i,j,ith,jth
		character(len=1)::w,divider,TensorType
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex(kind=4),pointer::cp(:)
		complex(kind=8),pointer::zp(:)
		logical,pointer::lp(:)
		integer,allocatable::location(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		cha=trim(adjustl(cha_))
		lenCha=len(cha)
 		TotalData=1
 		divider=array_character_divider
 		TensorType=cha(1:1)
 		w=cha(2:2)
 		if(w.nequ.'=')then
 			Res=cha
			if(present(dimen))call Res%resetdim(dimen)
 			return
 		end if
 		do i=3,lenCha
 			w=cha(i:i)
 			if(w.equ.divider)TotalData=TotalData+1
 		end do
		allocate(location(TotalData+1))
		j=1
		location(1)=2
		do i=3,lenCha
			w=cha(i:i)
			if(w.equ.divider)then
				if(j.ge.(TotalData+1))then
					call writemess('ERROR in converting character for Tensor=character')
					call error_stop
				end if
				j=j+1
				location(j)=i
			end if
		end do
		if(j.ge.(TotalData+1))then
			call writemess('ERROR in converting character for Tensor=character')
			call error_stop
		end if
		j=j+1
		location(j)=lenCha+1
		select case(TensorType)
			case ('i')
				call Res%allocate([TotalData],'integer')
				call Res%pointer(ip)
				do i=1,TotalData
					ith=location(i)+1
					jth=location(i+1)-1
					read(cha(ith:jth),*)ip(i)
				end do
			case ('s')
				call Res%allocate([TotalData],'real*4')
				call Res%pointer(sp)
				do i=1,TotalData
					ith=location(i)+1
					jth=location(i+1)-1
					read(cha(ith:jth),*)sp(i)
				end do
			case ('d')
				call Res%allocate([TotalData],'real*8')
				call Res%pointer(dp)
				do i=1,TotalData
					ith=location(i)+1
					jth=location(i+1)-1
					read(cha(ith:jth),*)dp(i)
				end do
			case ('c')
				call writemess('DO NO finished this type, in Tensor=character')
				call error_stop
			case ('z')
				call writemess('DO NO finished this type, in Tensor=character')
				call error_stop
			case ('l')
				call Res%allocate([TotalData],'logical')
				call Res%pointer(lp)
				do i=1,TotalData
					ith=location(i)+1
					jth=location(i+1)-1
					read(cha(ith:jth),*)lp(i)
				end do
			case ('a')
				call Res%allocate([TotalData],'character')
				do i=1,TotalData
					ith=location(i)+1
					jth=location(i+1)-1
					call Res%setValue(i,cha(ith:jth))
				end do
			case default
				call writemess('ERROR in converting character for Tensor=character')
				call error_stop
		end select
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_char_dim1(vec,dimen)result(Res)
		character(len=*),intent(in)::vec(:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_char_dim2(vec,dimen)result(Res)
		character(len=*),intent(in)::vec(:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_char_dim3(vec,dimen)result(Res)
		character(len=*),intent(in)::vec(:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function
	function constructor_char_dim4(vec,dimen)result(Res)
		character(len=*),intent(in)::vec(:,:,:,:)
		integer,optional,intent(in)::dimen(:)
		class(Tensor),allocatable::Res
		allocate(Tensor::Res)
		Res=Vec
		if(present(dimen))call Res%resetdim(dimen)
		return
	end function

	!********************************************************************************************
	!********************************************************************************************
	!                                    allocateTensor    
	!
	!*********************************************************************************************
	!*********************************************************************************************

	!allocate Tensor according to the dimen

	subroutine allocateTensor1(T,dimen,typede)
		class(Tensor),intent(inout) ::T
		type(dimension),intent(in)::dimen
		integer,intent(in)::typede
		integer::length
		if(T%getflag())then
			call writemess('Can not allocate to a allocated Tensor')
			call error_stop
		end if
		T%Dimension=dimen
		length=dimen%size()
		if(length.le.0) then
			write(*,*)"ERROR in allocateTensor"
			call dimen%print()
			call error_stop()
		end if
		call allocateData(T%TData,length,typede,.true.)
		return
	end subroutine

	!allocate Tensor according to the dimen

	subroutine allocateTensor2(T,dimen,typede)
		class(Tensor),intent(inout) ::T
		integer,intent(in)::dimen(:)
		integer,intent(in)::typede
		integer::length
		if(T%getflag())then
			call writemess('Can not allocate to a allocated Tensor')
			call error_stop
		end if
		T%Dimension=dimen
		length=product(dimen)
		if(length.le.0) then
			write(*,*)"ERROR in allocateTensor"
			write(*,*)dimen
			call error_stop()
		end if
		call allocateData(T%TData,length,typede,.true.)
		return
	end subroutine
	subroutine allocateTensor3(T,T2,typede)
		class(Tensor),intent(inout) ::T
		class(Tensor),intent(in) ::T2
		integer,intent(in)::typede
		integer::length
		if(T%getflag())then
			call writemess('Can not allocate to a allocated Tensor')
			call error_stop
		end if
		if(.not.T2%getflag())then
			call T%empty()
			return
		end if
		T%Dimension=T2%Dimension
		length=T2%TData%gettotalData()
		call allocateData(T%TData,length,typede,.true.)
		return
	end subroutine
	subroutine allocateTensor4(T,dimen,typede)
		class(Tensor),intent(inout) ::T
		type(dimension),intent(in)::dimen
		character(len=*),intent(in)::typede
		logical::deallocateflag
		integer::length
		if(T%getflag())then
			call writemess('Can not allocate to a allocated Tensor')
			call error_stop
		end if
		T%Dimension=dimen
		length=dimen%size()
		if(length.le.0) then
			write(*,*)"ERROR in allocateTensor"
			call dimen%print()
			call error_stop()
		end if
		call allocateData(T%TData,length,typede,.true.)
		return
	end subroutine

	!allocate Tensor according to the dimen

	subroutine allocateTensor5(T,dimen,typede)
		class(Tensor),intent(inout) ::T
		integer,intent(in)::dimen(:)
		character(len=*),intent(in)::typede
		logical::deallocateflag
		integer::length
		if(T%getflag())then
			call writemess('Can not allocate to a allocated Tensor')
			call error_stop
		end if
		T%Dimension=dimen
		length=product(dimen)
		if(length.le.0) then
			write(*,*)"ERROR in allocateTensor"
			write(*,*)dimen
			call error_stop()
		end if
		call allocateData(T%TData,length,typede,.true.)
		return
	end subroutine
	subroutine allocateTensor6(T,T2,typede)
		class(Tensor),intent(inout) ::T
		class(Tensor),intent(in) ::T2
		character(len=*),intent(in)::typede
		integer::length
		if(T%getflag())then
			call writemess('Can not allocate to a allocated Tensor')
			call error_stop
		end if
		if(.not.T2%getflag())then
			call T%empty()
			return
		end if
		T%Dimension=T2%Dimension
		length=T2%TData%gettotalData()
		call allocateData(T%TData,length,typede,.true.)
		return
	end subroutine
	subroutine allocateTensor7(T,dimen)
		class(Tensor),intent(inout) ::T
		type(dimension),intent(in)::dimen
		logical::deallocateflag
		integer::length
		if(T%getflag())then
			call writemess('Can not allocate to a allocated Tensor')
			call error_stop
		end if
		T%Dimension=dimen
		length=dimen%size()
		if(length.le.0) then
			write(*,*)"ERROR in allocateTensor"
			call dimen%print()
			call error_stop()
		end if
		call allocateData(T%TData,length,.false.)
		return
	end subroutine
	subroutine allocateTensor8(T,dimen)
		class(Tensor),intent(inout) ::T
		integer,intent(in)::dimen(:)
		logical::deallocateflag
		integer::length
		if(T%getflag())then
			call writemess('Can not allocate to a allocated Tensor')
			call error_stop
		end if
		T%Dimension=dimen
		length=product(dimen)
		if(length.le.0) then
			write(*,*)"ERROR in allocateTensor"
			write(*,*)dimen
			call error_stop()
		end if
		call allocateData(T%TData,length,.false.)
		return
	end subroutine
	subroutine allocateTensor9(T,T2)
		class(Tensor),intent(inout) ::T
		class(Tensor),intent(in) ::T2
		integer::length
		if(T%getflag())then
			call writemess('Can not allocate to a allocated Tensor')
			call error_stop
		end if
		if(.not.T2%getflag())then
			call T%empty()
			return
		end if
		T%Dimension=T2%dimension
		length=T2%TData%gettotalData()
		call allocateData(T%TData,length,T2%GetType(),.false.)
		return
	end subroutine


	!***********************************************
	!********* assignment **************************
	!***********************************************

	subroutine assignmentTenArray(T,T2)
		class(Tensor),intent(inout) ::T(:)
		class(Tensor),intent(in) :: T2(:)
		integer::length,i
		length=size(T2)
		if(size(T).lt.length)then
			write(*,*)"ERROR in assignment of two Tensor array "
			write(*,*)"T1(:)=T2(:),size(T1)<size(T2)"
			write(*,*)size(T),length
			call error_stop()
		end if
		do i=1,length
			T(i)=T2(i)
		end do
		return
	end subroutine

	subroutine DimInitialization2(Dimen,Dimen2)
		class(Tensor),intent(inout) ::Dimen
		class(Dimension),intent(in) ::Dimen2
		select type(Dimen2)
		class is (Tensor)
			if(.not.Dimen2%getFlag())then
				call Dimen%empty()
				return
			end if
			call Dimen%empty()
			call Dimen%allocate(dimen2)
			call assignmentTData_routine(Dimen%TData,Dimen2%TData)
		type is (TData)
			if(.not.Dimen2%getFlag())then
				call Dimen%empty()
				return
			end if
			call Dimen%empty()
			call Dimen%allocate([dimen2%getTotalData()],dimen2%getType())
			call assignmentTData_routine(Dimen%TData,Dimen2)
		type is (Dimension)
			if(.not.Dimen2%outlenDimData().eq.0)then
				call Dimen%empty()
				return
			end if
			call Dimen%empty()
			call Dimen%allocate(dimen2)
			call assignmentTData_int(Dimen%TData,Dimen2%dim(),Dimen2%size())
		class default
			call writemess('ERROR in assignment for tensor',-1)
			call error_stop
		end select
		return
	end subroutine
	subroutine DimInitialization(Dimen,DimData)
		class(Tensor),intent(inout) ::Dimen
		integer,intent(in) :: DimData(:)
		integer::length
		length=size(DimData)
		call Dimen%allocate([length],1)
		call assignmentTData_int(Dimen%TData,DimData,length)
		return
	end subroutine
	subroutine assignmentTenarray_int2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		integer,intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,1)
		call assignmentTData_int(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_int3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		integer,intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,1)
		call assignmentTData_int(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_int4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		integer,intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,1)
		call assignmentTData_int(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real4_1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call T%empty()
		call T%allocate([length],2)
		call assignmentTData_real4(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real4_2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,2)
		call assignmentTData_real4(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real4_3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,2)
		call assignmentTData_real4(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real4_4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,2)
		call assignmentTData_real4(T%TData,Tensor_data,length)
		return
	end subroutine
	
	subroutine assignmentTenarray_real8_1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call T%empty()
		call T%allocate([length],3)
		call assignmentTData_real8(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real8_2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,3)
		call assignmentTData_real8(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real8_3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,3)
		call assignmentTData_real8(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_real8_4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,3)
		call assignmentTData_real8(T%TData,Tensor_data,length)
		return
	end subroutine
	

	subroutine assignmentTenarray_com4_1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call T%empty()
		call T%allocate([length],4)
		call assignmentTData_com4(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com4_2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,4)
		call assignmentTData_com4(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com4_3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,4)
		call assignmentTData_com4(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com4_4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty()
		call T%allocate(dimen,4)
		call assignmentTData_com4(T%TData,Tensor_data,length)
		return
	end subroutine
	
	subroutine assignmentTenarray_com8_1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call T%empty()
		call T%allocate([length],5)
		call assignmentTData_com8(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com8_2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty() 
		call T%allocate(dimen,5)
		call assignmentTData_com8(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com8_3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty() 
		call T%allocate(dimen,5)
		call assignmentTData_com8(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_com8_4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty() 
		call T%allocate(dimen,5)
		call assignmentTData_com8(T%TData,Tensor_data,length)
		return
	end subroutine	

	

	subroutine assignmentTenarray_logi1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		logical,intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call T%empty() 
		call T%allocate([length],6)
		call assignmentTData_logi(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_logi2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		logical,intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty() 
		call T%allocate(dimen,6)
		call assignmentTData_logi(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_logi3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		logical,intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty() 
		call T%allocate(dimen,6)
		call assignmentTData_logi(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_logi4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		logical,intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty() 
		call T%allocate(dimen,6)
		call assignmentTData_logi(T%TData,Tensor_data,length)
		return
	end subroutine
	
	subroutine assignmentTenarray_char1(T,Tensor_data)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::Tensor_data(:)
		integer::length
		length=size(Tensor_data)
		call T%empty() 
		call T%allocate([length],7)
		call assignmentTData_char(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_char2(T,Tensor_data)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::Tensor_data(:,:)
		integer::length
		integer::dimen(2)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty() 
		call T%allocate(dimen,7)
		call assignmentTData_char(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_char3(T,Tensor_data)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::Tensor_data(:,:,:)
		integer::length
		integer::dimen(3)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty() 
		call T%allocate(dimen,7)
		call assignmentTData_char(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenarray_char4(T,Tensor_data)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::Tensor_data(:,:,:,:)
		integer::length
		integer::dimen(4)
		length=size(Tensor_data)
		dimen=shape(Tensor_data)
		call T%empty() 
		call T%allocate(dimen,7)
		call assignmentTData_char(T%TData,Tensor_data,length)
		return
	end subroutine
	subroutine assignmentTenNum_int(T,num)
		class(Tensor),intent(inout)::T
		integer,intent(in)::num
		call DimInitialization(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_real4(T,num)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::num
		call assignmentTenarray_real4_1(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_real8(T,num)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::num
		call assignmentTenarray_real8_1(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_com4(T,num)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::num
		call assignmentTenarray_com4_1(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_com8(T,num)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::num
		call assignmentTenarray_com8_1(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_logi(T,num)
		class(Tensor),intent(inout)::T
		logical,intent(in)::num
		call assignmentTenarray_logi1(T,[num])
		return
	end subroutine
	subroutine assignmentTenNum_char(T,num)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::num
		call assignmentTenarray_char1(T,[num])
		return
	end subroutine
	subroutine assignmentNumTen_int(val,T)
		integer,intent(inout)::val
		class(Tensor),intent(in)::T
		if(T%TData%getTotalData().eq.1) then
			call assignment_int_Tdata_value(val,T%TData)
		else 
			write(*,*)"ERROR in assignment for Tensor to integer"
			call error_stop()
		end if
		return
	end subroutine
	subroutine assignmentNumTen_real4(val,T)
		real(kind=4),intent(inout)::val
		class(Tensor),intent(in)::T
		if(T%TData%getTotalData().eq.1) then
			call assignment_real4_Tdata_value(val,T%TData)
		else if(.not.T%getFlag())then
			call writemess("ERROR in assignment for Tensor to real",-1)
			call writemess("The Tensor is a empty Tensor",-1)
			call error_stop()
		else
			call writemess("ERROR in assignment for Tensor to real",-1)
			call T%print()
			call error_stop()
		end if
		return
	end 	subroutine 
	subroutine assignmentNumTen_real8(val,T)
		real(kind=8),intent(inout)::val
		class(Tensor),intent(in)::T
		if(T%TData%getTotalData().eq.1) then
			call assignment_real8_Tdata_value(val,T%TData)
		else if(.not.T%getFlag())then
			call writemess("ERROR in assignment for Tensor to real",-1)
			call writemess("The Tensor is a empty Tensor",-1)
			call error_stop()
		else
			call writemess("ERROR in assignment for Tensor to real",-1)
			call T%print()
			call error_stop()
		end if
		return
	end subroutine
	subroutine assignmentNumTen_com4(val,T)
		complex(kind=4),intent(inout)::val
		class(Tensor),intent(in)::T
		if(T%TData%getTotalData().eq.1) then
			call assignment_com4_Tdata_value(val,T%TData)
		else
			call writemess("ERROR in assignment for Tensor to complex",-1)
			call error_stop()
		end if
		return
	end 	subroutine
	subroutine  assignmentNumTen_com8(val,T)
		complex(kind=8),intent(inout)::val
		class(Tensor),intent(in)::T
		if(T%TData%getTotalData().eq.1) then
			call assignment_com8_Tdata_value(val,T%TData)
		else
			write(*,*)"ERROR in assignment for Tensor to complex"
			call error_stop()
		end if
		return
	end subroutine
	subroutine  assignmentNumTen_logi(val,T)
		logical,intent(inout)::val
		class(Tensor),intent(in)::T
		if(T%TData%getTotalData().eq.1) then
			call assignment_logi_Tdata_value(val,T%TData)
		else
			write(*,*)"ERROR in assignment for Tensor to complex"
			call error_stop()
		end if
		return
	end subroutine
	subroutine  assignmentNumTen_char(val,T)
		character(len=*),intent(inout)::val
		class(Tensor),intent(in)::T
		if(T%TData%getTotalData().eq.1) then
			call assignment_char_Tdata_value(val,T%TData)
		else
			write(*,*)"ERROR in assignment for Tensor to complex"
			call error_stop()
		end if
		return
	end subroutine
	subroutine copyDimToVec(dimenVec,Dimen)
		integer,intent(inout) :: dimenVec(:)
		class(Tensor),intent(in) :: Dimen
		integer::length
		length=Dimen%TData%getTotalData()
		if(size(dimenVec).lt.length) then
			write(*,*)"The array can not store the Data"
			write(*,*)"Length of the array",size(dimenVec)
			write(*,*)"Length of the Tensor",length
			call error_stop()
		end if
		call assignment_int_Tdata(dimenVec,Dimen%TData,length)
		return
	end subroutine
	subroutine copyTensor_int_dim2(Mat,T)
		integer,intent(inout) ::Mat(:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,length
		if(T%Dimension%getRank().ne.2)then
			write(*,*)"T is not a rank-2 Tensor"
			write(*,*)"Can not copy to a matrix"
			call error_stop()
		end if
		length=T%TData%getTotalData()
		m=T.dim.1
		n=T.dim.2
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2)
			write(*,*)"size of the Tensor",m,n
			call error_stop()
		end if
		call assignment_int_Tdata(Mat,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_int_dim3(Mat,T)
		integer,intent(inout) ::Mat(:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,length
		if(T%Dimension%getRank().ne.3)then
			write(*,*)"T is not a rank-3 Tensor"
			write(*,*)"Can not copy to the array of dimension 3"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3)
			write(*,*)"size of the Tensor",m,n,l
			call error_stop()
		end if
		call assignment_int_Tdata(Mat,T%TData,length)
		return
	end subroutine	
	subroutine copyTensor_int_dim4(Mat,T)
		integer,intent(inout) ::Mat(:,:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,k,length
		if(T%Dimension%getRank().ne.4)then
			write(*,*)"T is not a rank-4 Tensor"
			write(*,*)"Can not copy to the array of dimension 4"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		k=T.dim.4
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l).or.(size(Mat,4).ne.k)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3),size(Mat,4)
			write(*,*)"size of the Tensor",m,n,l,k
			call error_stop()
		end if
		call assignment_int_Tdata(Mat,T%TData,length)
		return
	end subroutine		
	subroutine copyTensor_real4_dim1(Vec,T)
		real(kind=4),intent(inout) ::Vec(:)
		class(Tensor),intent(in) :: T
		integer::length
		length=T%TData%getTotalData()
		if(size(Vec).lt.length) then
			write(*,*)"The array can not store the Data"
			write(*,*)"Length of the array",size(Vec)
			write(*,*)"Length of the Tensor",length
			call error_stop()
		end if
		call assignment_real4_Tdata(vec,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_real4_dim2(Mat,T)
		real(kind=4),intent(inout) ::Mat(:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,length
		if(T%Dimension%getRank().ne.2)then
			write(*,*)"T is not a rank-2 Tensor"
			write(*,*)"Can not copy to a matrix"
			call error_stop()
		end if
		length=T%TData%getTotalData()
		m=T.dim.1
		n=T.dim.2
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2)
			write(*,*)"size of the Tensor",m,n
			call error_stop()
		end if
		call assignment_real4_Tdata(Mat,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_real4_dim3(Mat,T)
		real(kind=4),intent(inout) ::Mat(:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,length
		if(T%Dimension%getRank().ne.3)then
			write(*,*)"T is not a rank-3 Tensor"
			write(*,*)"Can not copy to the array of dimension 3"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3)
			write(*,*)"size of the Tensor",m,n,l
			call error_stop()
		end if
		call assignment_real4_Tdata(Mat,T%TData,length)
		return
	end subroutine	
	subroutine copyTensor_real4_dim4(Mat,T)
		real(kind=4),intent(inout) ::Mat(:,:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,k,length
		if(T%Dimension%getRank().ne.4)then
			write(*,*)"T is not a rank-4 Tensor"
			write(*,*)"Can not copy to the array of dimension 4"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		k=T.dim.4
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l).or.(size(Mat,4).ne.k)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3),size(Mat,4)
			write(*,*)"size of the Tensor",m,n,l,k
			call error_stop()
		end if
		call assignment_real4_Tdata(Mat,T%TData,length)
		return
	end subroutine		
	subroutine copyTensor_real8_dim1(Vec,T)
		real(kind=8),intent(inout) ::Vec(:)
		class(Tensor),intent(in) :: T
		integer::length
		length=T%TData%getTotalData()
		if(size(Vec).lt.length) then
			write(*,*)"The array can not store the Data"
			write(*,*)"Length of the array",size(Vec)
			write(*,*)"Length of the Tensor",length
			call error_stop()
		end if
		call assignment_real8_Tdata(vec,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_real8_dim2(Mat,T)
		real(kind=8),intent(inout) ::Mat(:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,length
		if(T%Dimension%getRank().ne.2)then
			write(*,*)"T is not a rank-2 Tensor"
			write(*,*)"Can not copy to a matrix"
			call error_stop()
		end if
		length=T%TData%getTotalData()
		m=T.dim.1
		n=T.dim.2
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2)
			write(*,*)"size of the Tensor",m,n
			call error_stop()
		end if
		call assignment_real8_Tdata(Mat,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_real8_dim3(Mat,T)
		real(kind=8),intent(inout) ::Mat(:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,length
		if(T%Dimension%getRank().ne.3)then
			write(*,*)"T is not a rank-3 Tensor"
			write(*,*)"Can not copy to the array of dimension 3"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3)
			write(*,*)"size of the Tensor",m,n,l
			call error_stop()
		end if
		call assignment_real8_Tdata(Mat,T%TData,length)
		return
	end subroutine	
	subroutine copyTensor_real8_dim4(Mat,T)
		real(kind=8),intent(inout) ::Mat(:,:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,k,length
		if(T%Dimension%getRank().ne.4)then
			write(*,*)"T is not a rank-4 Tensor"
			write(*,*)"Can not copy to the array of dimension 4"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		k=T.dim.4
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l).or.(size(Mat,4).ne.k)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3),size(Mat,4)
			write(*,*)"size of the Tensor",m,n,l,k
			call error_stop()
		end if
		call assignment_real8_Tdata(Mat,T%TData,length)
		return
	end subroutine		
	subroutine copyTensor_com4_dim1(Vec,T)
		complex(kind=4),intent(inout) ::Vec(:)
		class(Tensor),intent(in) :: T
		integer::length
		length=T%TData%getTotalData()
		if(size(Vec).lt.length) then
			write(*,*)"The array can not store the Data"
			write(*,*)"Length of the array",size(Vec)
			write(*,*)"Length of the Tensor",length
			call error_stop()
		end if
		call assignment_com4_Tdata(vec,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_com4_dim2(Mat,T)
		complex(kind=4),intent(inout) ::Mat(:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,length
		if(T%Dimension%getRank().ne.2)then
			write(*,*)"T is not a rank-2 Tensor"
			write(*,*)"Can not copy to a matrix"
			call error_stop()
		end if
		length=T%TData%getTotalData()
		m=T.dim.1
		n=T.dim.2
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2)
			write(*,*)"size of the Tensor",m,n
			call error_stop()
		end if
		call assignment_com4_Tdata(Mat,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_com4_dim3(Mat,T)
		complex(kind=4),intent(inout) ::Mat(:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,length
		if(T%Dimension%getRank().ne.3)then
			write(*,*)"T is not a rank-3 Tensor"
			write(*,*)"Can not copy to the array of dimension 3"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3)
			write(*,*)"size of the Tensor",m,n,l
			call error_stop()
		end if
		call assignment_com4_Tdata(Mat,T%TData,length)
		return
	end subroutine	
	subroutine copyTensor_com4_dim4(Mat,T)
		complex(kind=4),intent(inout) ::Mat(:,:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,k,length
		if(T%Dimension%getRank().ne.4)then
			write(*,*)"T is not a rank-4 Tensor"
			write(*,*)"Can not copy to the array of dimension 4"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		k=T.dim.4
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l).or.(size(Mat,4).ne.k)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3),size(Mat,4)
			write(*,*)"size of the Tensor",m,n,l,k
			call error_stop()
		end if
		call assignment_com4_Tdata(Mat,T%TData,length)
		return
	end subroutine		
	subroutine copyTensor_com8_dim1(Vec,T)
		complex(kind=8),intent(inout) ::Vec(:)
		class(Tensor),intent(in) :: T
		integer::length
		length=T%TData%getTotalData()
		if(size(Vec).lt.length) then
			write(*,*)"The array can not store the Data"
			write(*,*)"Length of the array",size(Vec)
			write(*,*)"Length of the Tensor",length
			call error_stop()
		end if
		call assignment_com8_Tdata(vec,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_com8_dim2(Mat,T)
		complex(kind=8),intent(inout) ::Mat(:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,length
		if(T%Dimension%getRank().ne.2)then
			write(*,*)"T is not a rank-2 Tensor"
			write(*,*)"Can not copy to a matrix"
			call error_stop()
		end if
		length=T%TData%getTotalData()
		m=T.dim.1
		n=T.dim.2
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2)
			write(*,*)"size of the Tensor",m,n
			call error_stop()
		end if
		call assignment_com8_Tdata(Mat,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_com8_dim3(Mat,T)
		complex(kind=8),intent(inout) ::Mat(:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,length
		if(T%Dimension%getRank().ne.3)then
			write(*,*)"T is not a rank-3 Tensor"
			write(*,*)"Can not copy to the array of dimension 3"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3)
			write(*,*)"size of the Tensor",m,n,l
			call error_stop()
		end if
		call assignment_com8_Tdata(Mat,T%TData,length)
		return
	end subroutine	
	subroutine copyTensor_com8_dim4(Mat,T)
		complex(kind=8),intent(inout) ::Mat(:,:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,k,length
		if(T%Dimension%getRank().ne.4)then
			write(*,*)"T is not a rank-4 Tensor"
			write(*,*)"Can not copy to the array of dimension 4"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		k=T.dim.4
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l).or.(size(Mat,4).ne.k)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3),size(Mat,4)
			write(*,*)"size of the Tensor",m,n,l,k
			call error_stop()
		end if
		call assignment_com8_Tdata(Mat,T%TData,length)
		return
	end subroutine		
	subroutine copyTensor_logi_dim1(Vec,T)
		logical,intent(inout) ::Vec(:)
		class(Tensor),intent(in) :: T
		integer::length
		length=T%TData%getTotalData()
		if(size(Vec).lt.length) then
			write(*,*)"The array can not store the Data"
			write(*,*)"Length of the array",size(Vec)
			write(*,*)"Length of the Tensor",length
			call error_stop()
		end if
		call assignment_logi_Tdata(vec,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_logi_dim2(Mat,T)
		logical,intent(inout) ::Mat(:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,length
		if(T%Dimension%getRank().ne.2)then
			write(*,*)"T is not a rank-2 Tensor"
			write(*,*)"Can not copy to a matrix"
			call error_stop()
		end if
		length=T%TData%getTotalData()
		m=T.dim.1
		n=T.dim.2
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2)
			write(*,*)"size of the Tensor",m,n
			call error_stop()
		end if
		call assignment_logi_Tdata(Mat,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_logi_dim3(Mat,T)
		logical,intent(inout) ::Mat(:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,length
		if(T%Dimension%getRank().ne.3)then
			write(*,*)"T is not a rank-3 Tensor"
			write(*,*)"Can not copy to the array of dimension 3"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3)
			write(*,*)"size of the Tensor",m,n,l
			call error_stop()
		end if
		call assignment_logi_Tdata(Mat,T%TData,length)
		return
	end subroutine	
	subroutine copyTensor_logi_dim4(Mat,T)
		logical,intent(inout) ::Mat(:,:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,k,length
		if(T%Dimension%getRank().ne.4)then
			write(*,*)"T is not a rank-4 Tensor"
			write(*,*)"Can not copy to the array of dimension 4"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		k=T.dim.4
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l).or.(size(Mat,4).ne.k)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3),size(Mat,4)
			write(*,*)"size of the Tensor",m,n,l,k
			call error_stop()
		end if
		call assignment_logi_Tdata(Mat,T%TData,length)
		return
	end subroutine		
	subroutine copyTensor_char_dim1(Vec,T)
		character(len=*),intent(inout) ::Vec(:)
		class(Tensor),intent(in) :: T
		integer::length
		length=T%TData%getTotalData()
		if(size(Vec).lt.length) then
			write(*,*)"The array can not store the Data"
			write(*,*)"Length of the array",size(Vec)
			write(*,*)"Length of the Tensor",length
			call error_stop()
		end if
		call assignment_char_Tdata(vec,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_char_dim2(Mat,T)
		character(len=*),intent(inout) ::Mat(:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,length
		if(T%Dimension%getRank().ne.2)then
			write(*,*)"T is not a rank-2 Tensor"
			write(*,*)"Can not copy to a matrix"
			call error_stop()
		end if
		length=T%TData%getTotalData()
		m=T.dim.1
		n=T.dim.2
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2)
			write(*,*)"size of the Tensor",m,n
			call error_stop()
		end if
		call assignment_char_Tdata(Mat,T%TData,length)
		return
	end subroutine
	subroutine copyTensor_char_dim3(Mat,T)
		character(len=*),intent(inout) ::Mat(:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,length
		if(T%Dimension%getRank().ne.3)then
			write(*,*)"T is not a rank-3 Tensor"
			write(*,*)"Can not copy to the array of dimension 3"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3)
			write(*,*)"size of the Tensor",m,n,l
			call error_stop()
		end if
		call assignment_char_Tdata(Mat,T%TData,length)
		return
	end subroutine	
	subroutine copyTensor_char_dim4(Mat,T)
		character(len=*),intent(inout) ::Mat(:,:,:,:)
		class(Tensor),intent(in) :: T
		integer::m,n,l,k,length
		if(T%Dimension%getRank().ne.4)then
			write(*,*)"T is not a rank-4 Tensor"
			write(*,*)"Can not copy to the array of dimension 4"
			call error_stop()
		end if
		m=T.dim.1
		n=T.dim.2
		l=T.dim.3
		k=T.dim.4
		length=T%TData%getTotalData()
		if((size(Mat,1).ne.m).or.(size(Mat,2).ne.n).or.(size(Mat,3).ne.l).or.(size(Mat,4).ne.k)) then
			write(*,*)"The array can not store the Data"
			write(*,*)"size of the array",size(Mat,1),size(Mat,2),size(Mat,3),size(Mat,4)
			write(*,*)"size of the Tensor",m,n,l,k
			call error_stop()
		end if
		call assignment_char_Tdata(Mat,T%TData,length)
		return
	end subroutine	

	!******************************************************************************
	!******************************************************************************
	!
	!                                  modify element in Tensor
	!
	!******************************************************************************
	!******************************************************************************

	subroutine modifyTen_val_class_i(Ten,dimen,val)!it will go wrong for character
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		integer,intent(in)::val
		integer::addre
		integer,allocatable::LDdimen(:)
		if(Ten%Dimension%getRank().eq.1) then
			if(dimen(1).gt.Ten%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(1))
			LDdimen=Ten%dim(1)
			call modify_TData_class(Ten,dimen,LDdimen,1,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.2) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt.Ten%Dim(2)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(2))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,2,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.3) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(3))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,3,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.4) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3))&
				.or.(dimen(4).gt.Ten%Dim(4)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(4))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,4,val)
			return
		end if
		addre=addressToIndes(Ten,dimen)
		if(addre.gt.Ten%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
			call error_stop()
		end if
		call modify_TData_class(Ten,(/addre/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine

	subroutine modifyTen_val_class_s(Ten,dimen,val)!it will go wrong for character
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		real*4,intent(in)::val
		integer::addre
		integer,allocatable::LDdimen(:)
		if(Ten%Dimension%getRank().eq.1) then
			if(dimen(1).gt.Ten%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(1))
			LDdimen=Ten%dim(1)
			call modify_TData_class(Ten,dimen,LDdimen,1,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.2) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt.Ten%Dim(2)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(2))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,2,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.3) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(3))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,3,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.4) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3))&
				.or.(dimen(4).gt.Ten%Dim(4)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(4))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,4,val)
			return
		end if
		addre=addressToIndes(Ten,dimen)
		if(addre.gt.Ten%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
			call error_stop()
		end if
		call modify_TData_class(Ten,(/addre/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine

	subroutine modifyTen_val_class_d(Ten,dimen,val)!it will go wrong for character
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		real*8,intent(in)::val
		integer::addre
		integer,allocatable::LDdimen(:)
		if(Ten%Dimension%getRank().eq.1) then
			if(dimen(1).gt.Ten%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(1))
			LDdimen=Ten%dim(1)
			call modify_TData_class(Ten,dimen,LDdimen,1,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.2) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt.Ten%Dim(2)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(2))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,2,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.3) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(3))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,3,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.4) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3))&
				.or.(dimen(4).gt.Ten%Dim(4)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(4))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,4,val)
			return
		end if
		addre=addressToIndes(Ten,dimen)
		if(addre.gt.Ten%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
			call error_stop()
		end if
		call modify_TData_class(Ten,(/addre/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine

	subroutine modifyTen_val_class_c(Ten,dimen,val)!it will go wrong for character
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		complex(kind=4),intent(in)::val
		integer::addre
		integer,allocatable::LDdimen(:)
		if(Ten%Dimension%getRank().eq.1) then
			if(dimen(1).gt.Ten%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(1))
			LDdimen=Ten%dim(1)
			call modify_TData_class(Ten,dimen,LDdimen,1,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.2) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt.Ten%Dim(2)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(2))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,2,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.3) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(3))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,3,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.4) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3))&
				.or.(dimen(4).gt.Ten%Dim(4)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(4))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,4,val)
			return
		end if
		addre=addressToIndes(Ten,dimen)
		if(addre.gt.Ten%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
			call error_stop()
		end if
		call modify_TData_class(Ten,(/addre/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine

	subroutine modifyTen_val_class_z(Ten,dimen,val)!it will go wrong for character
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		complex(kind=8),intent(in)::val
		integer::addre
		integer,allocatable::LDdimen(:)
		if(Ten%Dimension%getRank().eq.1) then
			if(dimen(1).gt.Ten%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(1))
			LDdimen=Ten%dim(1)
			call modify_TData_class(Ten,dimen,LDdimen,1,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.2) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt.Ten%Dim(2)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(2))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,2,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.3) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(3))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,3,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.4) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3))&
				.or.(dimen(4).gt.Ten%Dim(4)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(4))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,4,val)
			return
		end if
		addre=addressToIndes(Ten,dimen)
		if(addre.gt.Ten%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
			call error_stop()
		end if
		call modify_TData_class(Ten,(/addre/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine

	subroutine modifyTen_val_class_l(Ten,dimen,val)!it will go wrong for character
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		logical,intent(in)::val
		integer::addre
		integer,allocatable::LDdimen(:)
		if(Ten%Dimension%getRank().eq.1) then
			if(dimen(1).gt.Ten%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(1))
			LDdimen=Ten%dim(1)
			call modify_TData_class(Ten,dimen,LDdimen,1,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.2) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt.Ten%Dim(2)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(2))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,2,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.3) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(3))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,3,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.4) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3))&
				.or.(dimen(4).gt.Ten%Dim(4)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(4))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,4,val)
			return
		end if
		addre=addressToIndes(Ten,dimen)
		if(addre.gt.Ten%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
			call error_stop()
		end if
		call modify_TData_class(Ten,(/addre/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine

	subroutine modifyTen_val_class_a(Ten,dimen,val)!it will go wrong for character
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		character(len=*),intent(in)::val
		integer::addre
		integer,allocatable::LDdimen(:)
		if(Ten%Dimension%getRank().eq.1) then
			if(dimen(1).gt.Ten%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(1))
			LDdimen=Ten%dim(1)
			call modify_TData_class(Ten,dimen,LDdimen,1,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.2) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt.Ten%Dim(2)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(2))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,2,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.3) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(3))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,3,val)
			return
		end if
		if(Ten%Dimension%getRank().eq.4) then
			if( (dimen(1).gt.Ten%Dim(1)) .or. (dimen(2).gt. Ten%Dim(2)).or. (dimen(3).gt.Ten%dim(3))&
				.or.(dimen(4).gt.Ten%Dim(4)))Then
				call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
				call error_stop()
			end if
			allocate(LDdimen(4))
			LDdimen=Ten%dim()
			call modify_TData_class(Ten,dimen,LDdimen,4,val)
			return
		end if
		addre=addressToIndes(Ten,dimen)
		if(addre.gt.Ten%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+Ten%TData%getTotalData(),-1)
			call error_stop()
		end if
		call modify_TData_class(Ten,(/addre/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine

	subroutine modifyTen_val_Tensor(Ten,dimen,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		class(Tensor),intent(in)::val
		if(val%TData%getTotalData().ne.1)then
			call writemess("Do no finished this case, in modifyTen_val_Tensor",-1)
			call error_stop()
		end if
		if(.not.val%getFlag())then
			call writemess("There is no data in input element, setValue(element)",-1)
			call error_stop()
		end if	
		select case(val%getType())
			case (1)
				call modifyTen_val_class_i(Ten,dimen,val%ii(1))
			case (2)
				call modifyTen_val_class_s(Ten,dimen,val%si(1))	
			case (3)
				call modifyTen_val_class_d(Ten,dimen,val%di(1))
			case (4)
				call modifyTen_val_class_c(Ten,dimen,val%ci(1))
			case (5)
				call modifyTen_val_class_z(Ten,dimen,val%zi(1))
			case (6)
				call modifyTen_val_class_l(Ten,dimen,val%li(1))
			case (7)
				call modifyTen_val_class_a(Ten,dimen,val%ai(1))
			case default
				call writemess("ERROR in modifyTen_val_Tensor",-1)
				call error_stop()
		end select
		return
	end subroutine
	subroutine modifyTen_val_int(Ten,dimen,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		integer,intent(in)::val
		call modifyTen_val_class(Ten,dimen,val)
		return
	end subroutine
	subroutine modifyTen_val_real4(Ten,dimen,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		real(kind=4),intent(in)::val
		call modifyTen_val_class(Ten,dimen,val)
		return
	end subroutine
	subroutine modifyTen_val_real8(Ten,dimen,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		real(kind=8),intent(in)::val
		call modifyTen_val_class(Ten,dimen,val)
		return
	end subroutine
	subroutine modifyTen_val_com4(Ten,dimen,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		complex(kind=4),intent(in)::val
		call modifyTen_val_class(Ten,dimen,val)
		return
	end subroutine
	subroutine modifyTen_val_com8(Ten,dimen,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		complex(kind=8),intent(in)::val
		call modifyTen_val_class(Ten,dimen,val)
		return
	end subroutine
	subroutine modifyTen_val_logi(Ten,dimen,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		logical,intent(in)::val
		call modifyTen_val_class(Ten,dimen,val)
		return
	end subroutine
	subroutine modifyTen_val_char(Ten,dimen,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		character(len=*),intent(in)::val
		call modifyTen_val_class(Ten,dimen,val)
		return
	end subroutine

	
	subroutine modifyTen_val1_int(Ten,ith,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::ith
		integer,intent(in)::val
		call modify_TData_class(Ten,(/ith/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine
	subroutine modifyTen_val1_real4(Ten,ith,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::ith
		real(kind=4),intent(in)::val
		call modify_TData_class(Ten,(/ith/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine
	subroutine modifyTen_val1_real8(Ten,ith,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::ith
		real(kind=8),intent(in)::val
		call modify_TData_class(Ten,(/ith/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine
	subroutine modifyTen_val1_com4(Ten,ith,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::ith
		complex(kind=4),intent(in)::val
		call modify_TData_class(Ten,(/ith/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine
	subroutine modifyTen_val1_com8(Ten,ith,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::ith
		complex(kind=8),intent(in)::val
		call modify_TData_class(Ten,(/ith/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine
	subroutine modifyTen_val1_logi(Ten,ith,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::ith
		logical,intent(in)::val
		call modify_TData_class(Ten,(/ith/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine

	subroutine modifyTen_val1_char(Ten,ith,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::ith
		character(len=*),intent(in)::val
			call modify_TData_class(Ten,(/ith/),(/Ten%TData%getTotalData()/),1,val)
		return
	end subroutine
	subroutine modifyTen_val1_Tensor(Ten,dimen,val)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen
		class(Tensor),intent(in)::val
		if(val%TData%getTotalData().ne.1)then
			call writemess("Do no finished this case, in modifyTen_val_Tensor",-1)
			call error_stop()
		end if
		select case(val%getType())
			case (1)
				call modifyTen_val_class(Ten,(/dimen/),val%ii(1))
			case (2)
				call modifyTen_val_class(Ten,(/dimen/),val%si(1))	
			case (3)
				call modifyTen_val_class(Ten,(/dimen/),val%di(1))
			case (4)
				call modifyTen_val_class(Ten,(/dimen/),val%ci(1))
			case (5)
				call modifyTen_val_class(Ten,(/dimen/),val%zi(1))
			case (6)
				call modifyTen_val_class(Ten,(/dimen/),val%li(1))
			case (7)
				call modifyTen_val_class(Ten,(/dimen/),val%ai(1))
			case default
				call writemess("ERROR in modifyTen_val_Tensor",-1)
				call error_stop()
		end select
		return
	end subroutine
	subroutine store_T(T,Vec)
		class(Tensor),intent(inout) :: T
		class(Tensor),intent(in) ::Vec
		T%TData=Vec%TData
		return
	end subroutine	 
	subroutine store_int(T,Vec)
		class(Tensor),intent(inout) :: T
		integer,intent(in) ::Vec(*)
		integer::length
		length=T%TData%getTotalData()
		call assignmentTData_int(T%TData,Vec,length)
		return
	end subroutine	
	
	subroutine storeSome_int(T,i1,i2,Vec)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::i1,i2
		integer,intent(in) ::Vec(*)
		integer::length
		length=i2-i1+1
		call assignmentSomeTData_int(T%TData,vec,length,i1,i2)
		return
	end subroutine	
	subroutine storeSome_Tensor(T,i1,i2,Vec)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::i1,i2
		class(Tensor),intent(in) ::Vec
		integer::length
		if(.not.Vec%getFlag())then
			call writemess("There is no data in input element, setValue(element)",-1)
			call error_stop()
		end if	
		length=i2-i1+1
		call assignmentSomeTData_T(T%TData,Vec%TData,length,i1,i2)
		return
	end subroutine	
	
	subroutine store_real4(T,Vec)
		class(Tensor),intent(inout) :: T
		real(kind=4),intent(in) ::Vec(*)
		integer::length
		length=T%TData%getTotalData()
		call assignmentTData_real4(T%TData,Vec,length)
		return
	end subroutine	
	subroutine storeSome_real4(T,i1,i2,Vec)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::i1,i2
		real*4,intent(in) ::Vec(*)
		integer::length
		length=i2-i1+1
		call assignmentSomeTData_real4(T%TData,vec,length,i1,i2)
		return
	end subroutine	
	subroutine store_real8(T,Vec)
		class(Tensor),intent(inout) :: T
		real(kind=8),intent(in) ::Vec(*)
		integer::length
		length=T%TData%getTotalData()
		call assignmentTData_real8(T%TData,Vec,length)
		return
	end subroutine	
	subroutine storeSome_real8(T,i1,i2,Vec)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::i1,i2
		real*8,intent(in) ::Vec(*)
		integer::length
		length=i2-i1+1
		call assignmentSomeTData_real8(T%TData,vec,length,i1,i2)
		return
	end subroutine	
	subroutine store_com4(T,Vec)
		class(Tensor),intent(inout) :: T
		complex(kind=4),intent(in) ::Vec(*)
		integer::length
		length=T%TData%getTotalData()
		call assignmentTData_com4(T%TData,Vec,length)
		return
	end subroutine	
	subroutine storeSome_com4(T,i1,i2,Vec)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::i1,i2
		complex(kind=4),intent(in) ::Vec(*)
		integer::length
		length=i2-i1+1
		call assignmentSomeTData_com4(T%TData,vec,length,i1,i2)
		return
	end subroutine	
	subroutine store_com8(T,Vec)
		class(Tensor),intent(inout) :: T
		complex(kind=8),intent(in) ::Vec(*)
		integer::length
		length=T%TData%getTotalData()
		call assignmentTData_com8(T%TData,Vec,length)
		return
	end subroutine	
	subroutine storeSome_com8(T,i1,i2,Vec)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::i1,i2
		complex(kind=8),intent(in) ::Vec(*)
		integer::length
		length=i2-i1+1
		call assignmentSomeTData_com8(T%TData,vec,length,i1,i2)
		return
	end subroutine	
	subroutine store_logi(T,Vec)
		class(Tensor),intent(inout) :: T
		logical,intent(in) ::Vec(*)
		integer::length
		length=T%TData%getTotalData()
		call assignmentTData_logi(T%TData,Vec,length)
		return
	end subroutine	
	subroutine storeSome_logi(T,i1,i2,Vec)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::i1,i2
		logical,intent(in) ::Vec(*)
		integer::length
		length=i2-i1+1
		call assignmentSomeTData_logi(T%TData,vec,length,i1,i2)
		return
	end subroutine	
	subroutine store_char(T,Vec)
		class(Tensor),intent(inout) :: T
		character(len=*),intent(in) ::Vec(*)
		integer::length
		length=T%TData%getTotalData()
		call assignmentTData_char(T%TData,Vec,length)
		return
	end subroutine	
	subroutine storeSome_char(T,i1,i2,Vec)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::i1,i2
		character(len=*),intent(in) ::Vec(*)
		integer::length
		length=i2-i1+1
		call assignmentSomeTData_char(T%TData,vec,length,i1,i2)
		return
	end subroutine	
	
	subroutine modify_chech_input(A,ia)
		class(Tensor)::A
		integer::ia(2)
		if(ia(1).gt.ia(2))then
			call writemess('ERROR in input parameter, when setvalue in Tensor.f90',-1)
			call writemess('ERROR:ia(1)>ia(2)'+',ia(1)='+ia(1)+',ia(2)='+ia(2),-1)
			call error_stop
		end if
		if(ia(1).le.0)then
			call writemess('ERROR in input parameter, when setvalue in Tensor.f90',-1)
			call writemess('ERROR:ia(1)<0'+',ia(1)='+ia(1))
			call error_stop
		end if
		if(ia(2).gt.A%TData%getToTalData())then
			call writemess('ERROR in input parameter, when setvalue in Tensor.f90',-1)
			call writemess('ERROR:ia(2)>len of Tensor Data'+',ia(2)='+ia(2)+',TotalData='+A%TData%getToTalData(),-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine modify_chech_inputParameter(ia,lenA)
		integer::ia(2),lenA
		if(ia(1).gt.ia(2))then
			call writemess('ERROR in input parameter, when setvalue in Tensor.f90',-1)
			call writemess('ERROR:ia(1)>ia(2)'+',ia(1)='+ia(1)+',ia(2)='+ia(2),-1)
			call error_stop
		end if
		if(ia(1).le.0)then
			call writemess('ERROR in input parameter, when setvalue in Tensor.f90',-1)
			call writemess('ERROR:ia(1)<0'+',ia(1)='+ia(1),-1)
			call error_stop
		end if
		if(ia(2).gt.lenA)then
			call writemess('ERROR in input parameter, when setvalue in Tensor.f90',-1)
			call writemess('ERROR:ia(2)>len of lenA Data'+',ia(2)='+ia(2)+',TotalData='+lenA,-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine modify_chech_inputindex(ia,ib)
		integer::ia(2),ib(2)
		if((ia(2)-ia(1)).ne.(ib(2)-ib(1)))then
			call writemess('ERROR in input parameter, when setvalue in Tensor.f90',-1)
			call writemess('ERROR: index do not match,ia(2)-ia(1)!=ib(2)-ib(1)',-1)
			call writemess('ia(1)='+ia(1),-1)
			call writemess('ia(2)='+ia(2),-1)
			call writemess('ib(1)='+ib(1),-1)
			call writemess('ib(2)='+ib(2),-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine modify_chech_inputindex2(ia,ja,total)
		integer::ia(2),ja(2),total
		integer::m,n
		m=ia(2)-ia(1)+1
		n=ja(2)-ja(1)+1
		if((m*n).ne.total)then
			call writemess('ERROR in input parameter, when setvalue in Tensor.f90',-1)
			call writemess('ERROR: index do not match,(ia(2)-ia(1))*(ja(2)-ja(1))!=total',-1)
			call writemess('ia(1)='+ia(1),-1)
			call writemess('ia(2)='+ia(2),-1)
			call writemess('ja(1)='+ja(1),-1)
			call writemess('ja(2)='+ja(2),-1)
			call writemess('total='+total,-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine modify_chech_inputindex1(ia,total)
		integer::ia(2),total
		integer::m
		m=ia(2)-ia(1)
		if(m.ne.total)then
			call writemess('ERROR in input parameter, when setvalue in Tensor.f90',-1)
			call writemess('ERROR: index do not match,(ia(2)-ia(1))!=total',-1)
			call writemess('ia(1)='+ia(1),-1)
			call writemess('ia(2)='+ia(2),-1)
			call writemess('total='+total,-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine modify_Some_Data1_int(A,ia,B,ib)	
		class(Tensor),intent(inout) :: A
		integer,intent(in)::B(:)
		integer,intent(in)::ia(2),ib(2)
		integer::lenB
		call modify_chech_input(A,ia)
		lenB=size(B)
		call modify_chech_inputParameter(ib,lenB)
		call modify_chech_inputindex(ia,ib)
		call modify_Some_TData_class1(A%TData,A%TData%GetTotalData(),ia,B,lenB,ib)
		return
	end subroutine
	subroutine modify_Some_Data1_real4(A,ia,B,ib)	
		class(Tensor),intent(inout) :: A
		real*4,intent(in)::B(:)
		integer,intent(in)::ia(2),ib(2)
		integer::lenB
		call modify_chech_input(A,ia)
		lenB=size(B)
		call modify_chech_inputParameter(ib,lenB)
		call modify_chech_inputindex(ia,ib)
		call modify_Some_TData_class1(A%TData,A%TData%GetTotalData(),ia,B,lenB,ib)
		return
	end subroutine
	subroutine modify_Some_Data1_real8(A,ia,B,ib)	
		class(Tensor),intent(inout) :: A
		real*8,intent(in)::B(:)
		integer,intent(in)::ia(2),ib(2)
		integer::lenB
		call modify_chech_input(A,ia)
		lenB=size(B)
		call modify_chech_inputParameter(ib,lenB)
		call modify_chech_inputindex(ia,ib)
		call modify_Some_TData_class1(A%TData,A%TData%GetTotalData(),ia,B,lenB,ib)
		return
	end subroutine
	subroutine modify_Some_Data1_com4(A,ia,B,ib)	
		class(Tensor),intent(inout) :: A
		complex*8,intent(in)::B(:)
		integer,intent(in)::ia(2),ib(2)
		integer::lenB
		call modify_chech_input(A,ia)
		lenB=size(B)
		call modify_chech_inputParameter(ib,lenB)
		call modify_chech_inputindex(ia,ib)
		call modify_Some_TData_class1(A%TData,A%TData%GetTotalData(),ia,B,lenB,ib)
		return
	end subroutine
	subroutine modify_Some_Data1_com8(A,ia,B,ib)	
		class(Tensor),intent(inout) :: A
		complex*16,intent(in)::B(:)
		integer,intent(in)::ia(2),ib(2)
		integer::lenB
		call modify_chech_input(A,ia)
		lenB=size(B)
		call modify_chech_inputParameter(ib,lenB)
		call modify_chech_inputindex(ia,ib)
		call modify_Some_TData_class1(A%TData,A%TData%GetTotalData(),ia,B,lenB,ib)
		return
	end subroutine
	subroutine modify_Some_Data1_logi(A,ia,B,ib)	
		class(Tensor),intent(inout) :: A
		logical,intent(in)::B(:)
		integer,intent(in)::ia(2),ib(2)
		integer::lenB
		call modify_chech_input(A,ia)
		lenB=size(B)
		call modify_chech_inputParameter(ib,lenB)
		call modify_chech_inputindex(ia,ib)
		call modify_Some_TData_class1(A%TData,A%TData%GetTotalData(),ia,B,lenB,ib)
		return
	end subroutine
	subroutine modify_Some_Data1_char(A,ia,B,ib)	
		class(Tensor),intent(inout) :: A
		character(len=*),intent(in)::B(:)
		integer,intent(in)::ia(2),ib(2)
		integer::lenB
		call modify_chech_input(A,ia)
		lenB=size(B)
		call modify_chech_inputParameter(ib,lenB)
		call modify_chech_inputindex(ia,ib)
		call modify_Some_TData_class1(A%TData,A%TData%GetTotalData(),ia,B,lenB,ib)
		return
	end subroutine
	subroutine modify_Some_Data1_Tensor(A,ia,B,ib)	
		class(Tensor),intent(inout) :: A
		class(Tensor),intent(in)::B
		integer,intent(in)::ia(2),ib(2)
		call modify_chech_input(A,ia)
		call modify_chech_input(B,ib)
		call modify_chech_inputindex(ia,ib)
		call modify_Some_TData_TData1(A%TData,A%TData%GetTotalData(),ia,B%TData,A%TData%GetTotalData(),ib)
		return
	end subroutine
	
	subroutine modify_Some_Data2_int(A,ia,ja,B,ib,jb)	
		class(Tensor),intent(inout) :: A
		integer,intent(in)::B(:,:)
		integer,intent(in)::ia(2),ja(2),ib(2),jb(2)
		integer::Bdim(2),Adim(2)
		if(A%Dimension%getRank().ne.2)then
			call writemess('ERROR in SETTing value in Tensor.f90',-1)
			call writemess('The rank of the Tensor should be 2',-1)
			call error_stop
		end if
		Adim=A%dim()
		Bdim(1)=size(B,1)
		Bdim(2)=size(B,2)
		call modify_chech_inputParameter(ia,Adim(1))
		call modify_chech_inputParameter(ja,Adim(2))
		call modify_chech_inputParameter(ib,Bdim(1))
		call modify_chech_inputParameter(jb,Bdim(2))
		call modify_chech_inputindex(ia,ib)
		call modify_chech_inputindex(ja,jb)
		call modify_Some_TData_class2(A%TData,Adim,ia,ja,B,Bdim,ib,jb)
		return
	end subroutine
	subroutine modify_Some_Data2_real4(A,ia,ja,B,ib,jb)	
		class(Tensor),intent(inout) :: A
		real*4,intent(in)::B(:,:)
		integer,intent(in)::ia(2),ja(2),ib(2),jb(2)
		integer::Bdim(2),Adim(2)
		if(A%Dimension%getRank().ne.2)then
			call writemess('ERROR in SETTing value in Tensor.f90',-1)
			call writemess('The rank of the Tensor should be 2',-1)
			call error_stop
		end if
		Adim=A%dim()
		Bdim(1)=size(B,1)
		Bdim(2)=size(B,2)
		call modify_chech_inputParameter(ia,Adim(1))
		call modify_chech_inputParameter(ja,Adim(2))
		call modify_chech_inputParameter(ib,Bdim(1))
		call modify_chech_inputParameter(jb,Bdim(2))
		call modify_chech_inputindex(ia,ib)
		call modify_chech_inputindex(ja,jb)
		call modify_Some_TData_class2(A%TData,Adim,ia,ja,B,Bdim,ib,jb)
		return
	end subroutine
	subroutine modify_Some_Data2_real8(A,ia,ja,B,ib,jb)	
		class(Tensor),intent(inout) :: A
		real*8,intent(in)::B(:,:)
		integer,intent(in)::ia(2),ja(2),ib(2),jb(2)
		integer::Bdim(2),Adim(2)
		if(A%Dimension%getRank().ne.2)then
			call writemess('ERROR in SETTing value in Tensor.f90',-1)
			call writemess('The rank of the Tensor should be 2',-1)
			call error_stop
		end if
		Adim=A%dim()
		Bdim(1)=size(B,1)
		Bdim(2)=size(B,2)
		call modify_chech_inputParameter(ia,Adim(1))
		call modify_chech_inputParameter(ja,Adim(2))
		call modify_chech_inputParameter(ib,Bdim(1))
		call modify_chech_inputParameter(jb,Bdim(2))
		call modify_chech_inputindex(ia,ib)
		call modify_chech_inputindex(ja,jb)
		call modify_Some_TData_class2(A%TData,Adim,ia,ja,B,Bdim,ib,jb)
		return
	end subroutine
	subroutine modify_Some_Data2_com4(A,ia,ja,B,ib,jb)	
		class(Tensor),intent(inout) :: A
		complex*8,intent(in)::B(:,:)
		integer,intent(in)::ia(2),ja(2),ib(2),jb(2)
		integer::Bdim(2),Adim(2)
		if(A%Dimension%getRank().ne.2)then
			call writemess('ERROR in SETTing value in Tensor.f90',-1)
			call writemess('The rank of the Tensor should be 2',-1)
			call error_stop
		end if
		Adim=A%dim()
		Bdim(1)=size(B,1)
		Bdim(2)=size(B,2)
		call modify_chech_inputParameter(ia,Adim(1))
		call modify_chech_inputParameter(ja,Adim(2))
		call modify_chech_inputParameter(ib,Bdim(1))
		call modify_chech_inputParameter(jb,Bdim(2))
		call modify_chech_inputindex(ia,ib)
		call modify_chech_inputindex(ja,jb)
		call modify_Some_TData_class2(A%TData,Adim,ia,ja,B,Bdim,ib,jb)
		return
	end subroutine
	subroutine modify_Some_Data2_com8(A,ia,ja,B,ib,jb)	
		class(Tensor),intent(inout) :: A
		complex*16,intent(in)::B(:,:)
		integer,intent(in)::ia(2),ja(2),ib(2),jb(2)
		integer::Bdim(2),Adim(2)
		if(A%Dimension%getRank().ne.2)then
			call writemess('ERROR in SETTing value in Tensor.f90',-1)
			call writemess('The rank of the Tensor should be 2',-1)
			call error_stop
		end if
		Adim=A%dim()
		Bdim(1)=size(B,1)
		Bdim(2)=size(B,2)
		call modify_chech_inputParameter(ia,Adim(1))
		call modify_chech_inputParameter(ja,Adim(2))
		call modify_chech_inputParameter(ib,Bdim(1))
		call modify_chech_inputParameter(jb,Bdim(2))
		call modify_chech_inputindex(ia,ib)
		call modify_chech_inputindex(ja,jb)
		call modify_Some_TData_class2(A%TData,Adim,ia,ja,B,Bdim,ib,jb)
		return
	end subroutine
	subroutine modify_Some_Data2_logi(A,ia,ja,B,ib,jb)	
		class(Tensor),intent(inout) :: A
		logical,intent(in)::B(:,:)
		integer,intent(in)::ia(2),ja(2),ib(2),jb(2)
		integer::Bdim(2),Adim(2)
		if(A%Dimension%getRank().ne.2)then
			call writemess('ERROR in SETTing value in Tensor.f90',-1)
			call writemess('The rank of the Tensor should be 2',-1)
			call error_stop
		end if
		Adim=A%dim()
		Bdim(1)=size(B,1)
		Bdim(2)=size(B,2)
		call modify_chech_inputParameter(ia,Adim(1))
		call modify_chech_inputParameter(ja,Adim(2))
		call modify_chech_inputParameter(ib,Bdim(1))
		call modify_chech_inputParameter(jb,Bdim(2))
		call modify_chech_inputindex(ia,ib)
		call modify_chech_inputindex(ja,jb)
		call modify_Some_TData_class2(A%TData,Adim,ia,ja,B,Bdim,ib,jb)
		return
	end subroutine
	subroutine modify_Some_Data2_char(A,ia,ja,B,ib,jb)	
		class(Tensor),intent(inout) :: A
		character(len=*),intent(in)::B(:,:)
		integer,intent(in)::ia(2),ja(2),ib(2),jb(2)
		integer::Bdim(2),Adim(2)
		if(A%Dimension%getRank().ne.2)then
			call writemess('ERROR in SETTing value in Tensor.f90',-1)
			call writemess('The rank of the Tensor should be 2',-1)
			call error_stop
		end if
		Adim=A%dim()
		Bdim(1)=size(B,1)
		Bdim(2)=size(B,2)
		call modify_chech_inputParameter(ia,Adim(1))
		call modify_chech_inputParameter(ja,Adim(2))
		call modify_chech_inputParameter(ib,Bdim(1))
		call modify_chech_inputParameter(jb,Bdim(2))
		call modify_chech_inputindex(ia,ib)
		call modify_chech_inputindex(ja,jb)
		call modify_Some_TData_class2(A%TData,Adim,ia,ja,B,Bdim,ib,jb)
		return
	end subroutine
	subroutine modify_Some_Data2_Tensor(A,ia,ja,B,ib,jb)	
		class(Tensor),intent(inout) :: A
		class(Tensor),intent(in)::B
		integer,intent(in)::ia(2),ja(2),ib(2),jb(2)
		integer::Bdim(2),Adim(2)
		if(A%Dimension%getRank().ne.2)then
			call writemess('ERROR in SETTing value in Tensor.f90',-1)
			call writemess('The rank of the Tensor should be 2',-1)
			call error_stop
		end if
		if(B%Dimension%getRank().ne.2)then
			call writemess('ERROR in SETTing value in Tensor.f90',-1)
			call writemess('The rank of the Tensor should be 2',-1)
			call error_stop
		end if
		Adim=A%dim()
		Bdim=B%dim()
		call modify_chech_inputParameter(ia,Adim(1))
		call modify_chech_inputParameter(ja,Adim(2))
		call modify_chech_inputParameter(ib,Bdim(1))
		call modify_chech_inputParameter(jb,Bdim(2))
		call modify_chech_inputindex(ia,ib)
		call modify_chech_inputindex(ja,jb)
		call modify_Some_TData_TData2(A%TData,Adim,ia,ja,B%TData,Bdim,ib,jb)
		return
	end subroutine
	
	subroutine modify_Some_Data2_Tensor2(A,ia,ja,B)	
		class(Tensor),intent(inout) :: A
		class(Tensor),intent(in)::B
		integer,intent(in)::ia(2),ja(2)
		integer::Bdim(2),Adim(2),ib(2),jb(2),total
		if(A%Dimension%getRank().ne.2)then
			call writemess('ERROR in SETTing value in Tensor.f90',-1)
			call writemess('The rank of the Tensor should be 2',-1)
			call error_stop
		end if
		Adim=A%dim()
		total=B%TData%getTotalData()
		call modify_chech_inputParameter(ia,Adim(1))
		call modify_chech_inputParameter(ja,Adim(2))
		call modify_chech_inputindex2(ia,ja,total)
		Bdim(1)=ia(2)-ia(1)+1
		Bdim(2)=ja(2)-ja(1)+1
		ib(1)=1
		ib(2)=Bdim(1)
		jb(1)=1
		jb(2)=Bdim(2)
		call modify_Some_TData_TData2(A%TData,Adim,ia,ja,B%TData,Bdim,ib,jb)
		return
	end subroutine


	subroutine set_value_Tensor_class(T,value)
		class(Tensor),intent(inout)::T
		class(*),intent(in)::value
		if(.not.T%getflag()) then
			return
		end if
		select type(value)
			type is (integer)
				call set_all_data_int(T%TData,value)
			type is (real(kind=4))
				call set_all_data_real4(T%TData,value)
			type is (real(kind=8))
				call set_all_data_real8(T%TData,value)
			type is (complex(kind=4))
				call set_all_data_com4(T%TData,value)
			type is (complex(kind=8))
				call set_all_data_com8(T%TData,value)
			type is (logical)
				call set_all_data_logi(T%TData,value)
		end select
		return
	end subroutine
	subroutine set_zero_Tensor(T)
		class(Tensor),intent(inout)::T
		if(.not.T%getflag()) then
			return
		end if
		call set_all_data_int(T%TData,0)
		return
	end subroutine
	subroutine set_value_Tensor_int(T,value)
		class(Tensor),intent(inout)::T
		integer,intent(in)::value
		if(.not.T%getflag()) then
			return
		end if
		call set_all_data_int(T%TData,value)
		return
	end subroutine
	subroutine set_value_Tensor_real4(T,value)
		class(Tensor),intent(inout)::T
		real(kind=4),intent(in)::value
		if(.not.T%getflag()) then
			return
		end if
		call set_all_data_real4(T%TData,value)
		return
	end subroutine
	subroutine set_value_Tensor_real8(T,value)
		class(Tensor),intent(inout)::T
		real(kind=8),intent(in)::value
		if(.not.T%getflag()) then
			return
		end if
		call set_all_data_real8(T%TData,value)
		return
	end subroutine
	subroutine set_value_Tensor_com4(T,value)
		class(Tensor),intent(inout)::T
		complex(kind=4),intent(in)::value
		if(.not.T%getflag()) then
			return
		end if
		call set_all_data_com4(T%TData,value)
		return
	end subroutine
	subroutine set_value_Tensor_com8(T,value)
		class(Tensor),intent(inout)::T
		complex(kind=8),intent(in)::value
		if(.not.T%getflag()) then
			return
		end if
		call set_all_data_com8(T%TData,value)
		return
	end subroutine
	subroutine set_value_Tensor_logi(T,value)
		class(Tensor),intent(inout)::T
		logical,intent(in)::value
		if(.not.T%getflag()) then
			return
		end if
		call set_all_data_logi(T%TData,value)
		return
	end subroutine
	subroutine set_value_Tensor_char(T,value)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::value
		if(.not.T%getflag()) then
			return
		end if
		call set_all_data_char(T%TData,value)
		return
	end subroutine

	!******************************************************************************
	!******************************************************************************
	!
	!                    set random elements in Tensors    
	!
	!******************************************************************************
	!******************************************************************************

	subroutine set_random_Tensor(T,region)
		class(Tensor),intent(inout)::T
		real*8,optional,intent(in)::region(*)
		if(.not.T%getflag()) then
			return
		end if
		call generate_random_data(T%TData,region)
		return
	end subroutine
	subroutine set_random_Tensor_region4(T,region)
		class(Tensor),intent(inout)::T
		real*4,intent(in)::region(*)
		if(.not.T%getflag()) then
			return
		end if
		call generate_random_data_region4(T%TData,region)
		return
	end subroutine
	subroutine set_random_Tensor_regioni(T,region)
		class(Tensor),intent(inout)::T
		integer,intent(in)::region(*)
		if(.not.T%getflag()) then
			return
		end if
		call generate_random_data_regioni(T%TData,region)
		return
	end subroutine



	!******************************************************************************
	!******************************************************************************
	!
	!                                    print Tensor    
	!
	!******************************************************************************
	!******************************************************************************

	subroutine writemess_Tensor(mess,cpu_number)!overwrite writemess
		type(Tensor),intent(in)::mess
		integer,optional,intent(in)::cpu_number
		integer::i,j,k,totoal,rank
		integer,pointer::idata2(:,:),idata3(:,:,:),idata4(:,:,:,:)
		real*4,pointer::sdata2(:,:),sdata3(:,:,:),sdata4(:,:,:,:)
		real*8,pointer::ddata2(:,:),ddata3(:,:,:),ddata4(:,:,:,:)
		complex*8,pointer::cdata2(:,:),cdata3(:,:,:),cdata4(:,:,:,:)
		complex*16,pointer::zdata2(:,:),zdata3(:,:,:),zdata4(:,:,:,:)
		logical,pointer::ldata2(:,:),ldata3(:,:,:),ldata4(:,:,:,:)
		character(len=max_len_of_char)::w
		w=''
		if(.not.mess%getFlag())then
			call writemess('There is no data in the Tensor',cpu_number)
			return
		end if
		totoal=mess%TData%getTotalData()
		if(totoal.eq.0)then
			call writemess('There is no data in the Tensor',cpu_number)
			return
		end if
		rank=mess%Dimension%getRank()
		select case(mess%getType())
			case(1)
				select case(rank)
					case(1)
						call writemess(mess%ii(),cpu_number)
					case(2)
						call mess%pointer(idata2)
						do i=1,mess%dim(1)
							call writemess(idata2(i,:),cpu_number)
						end do
					case(3)
						call mess%pointer(idata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(idata3(i,:,j),cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(idata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(idata4(i,:,j,k),cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%ii(),cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(2)
				select case(rank)
					case(1)
						call writemess(mess%si(),cpu_number)
					case(2)
						call mess%pointer(sdata2)
						do i=1,mess%dim(1)
							call writemess(sdata2(i,:),cpu_number)
						end do
					case(3)
						call mess%pointer(sdata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(sdata3(i,:,j),cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(sdata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(sdata4(i,:,j,k),cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%si(),cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(3)
				select case(rank)
					case(1)
						call writemess(mess%di(),cpu_number)
					case(2)
						call mess%pointer(ddata2)
						do i=1,mess%dim(1)
							call writemess(ddata2(i,:),cpu_number)
						end do
					case(3)
						call mess%pointer(ddata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(ddata3(i,:,j),cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(ddata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(ddata4(i,:,j,k),cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%di(),cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(4)
				select case(rank)
					case(1)
						call writemess(mess%ci(),cpu_number)
					case(2)
						call mess%pointer(cdata2)
						do i=1,mess%dim(1)
							call writemess(cdata2(i,:),cpu_number)
						end do
					case(3)
						call mess%pointer(cdata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(cdata3(i,:,j),cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(cdata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(cdata4(i,:,j,k),cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%ci(),cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(5)
				select case(rank)
					case(1)
						call writemess(mess%zi(),cpu_number)
					case(2)
						call mess%pointer(zdata2)
						do i=1,mess%dim(1)
							call writemess(zdata2(i,:),cpu_number)
						end do
					case(3)
						call mess%pointer(zdata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(zdata3(i,:,j),cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(zdata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(zdata4(i,:,j,k),cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%zi(),cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(6)
				select case(rank)
					case(1)
						call writemess(mess%li(),cpu_number)
					case(2)
						call mess%pointer(ldata2)
						do i=1,mess%dim(1)
							call writemess(ldata2(i,:),cpu_number)
						end do
					case(3)
						call mess%pointer(ldata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(ldata3(i,:,j),cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(ldata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(ldata4(i,:,j,k),cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%li(),cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(7)
				do i=1,totoal-1
					w=w+mess%ai(i)+','
				end do
				w=w+mess%ai(totoal)
				call writemess(w,cpu_number)
		end select
		return
	end subroutine

	subroutine writemess_Tensor_form(mess,form,cpu_number)!overwrite writemess
		type(Tensor),intent(in)::mess
		character(len=*),intent(in)::form
		integer,optional,intent(in)::cpu_number
		integer::i,j,k,totoal,rank
		integer,pointer::idata2(:,:),idata3(:,:,:),idata4(:,:,:,:)
		real*4,pointer::sdata2(:,:),sdata3(:,:,:),sdata4(:,:,:,:)
		real*8,pointer::ddata2(:,:),ddata3(:,:,:),ddata4(:,:,:,:)
		complex*8,pointer::cdata2(:,:),cdata3(:,:,:),cdata4(:,:,:,:)
		complex*16,pointer::zdata2(:,:),zdata3(:,:,:),zdata4(:,:,:,:)
		logical,pointer::ldata2(:,:),ldata3(:,:,:),ldata4(:,:,:,:)
		character(len=max_len_of_char)::w
		w=''
		if(.not.mess%getFlag())then
			call writemess('There is no data in the Tensor',cpu_number)
			return
		end if
		totoal=mess%TData%getTotalData()
		if(totoal.eq.0)then
			call writemess('There is no data in the Tensor',cpu_number)
			return
		end if
		rank=mess%Dimension%getRank()
		select case(mess%getType())
			case(1)
				select case(rank)
					case(1)
						call writemess(mess%ii(),form,cpu_number)
					case(2)
						call mess%pointer(idata2)
						do i=1,mess%dim(1)
							call writemess(idata2(i,:),form,cpu_number)
						end do
					case(3)
						call mess%pointer(idata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(idata3(i,:,j),form,cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(idata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(idata4(i,:,j,k),form,cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%ii(),form,cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(2)
				select case(rank)
					case(1)
						call writemess(mess%si(),form,cpu_number)
					case(2)
						call mess%pointer(sdata2)
						do i=1,mess%dim(1)
							call writemess(sdata2(i,:),form,cpu_number)
						end do
					case(3)
						call mess%pointer(sdata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(sdata3(i,:,j),form,cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(sdata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(sdata4(i,:,j,k),form,cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%si(),form,cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(3)
				select case(rank)
					case(1)
						call writemess(mess%di(),form,cpu_number)
					case(2)
						call mess%pointer(ddata2)
						do i=1,mess%dim(1)
							call writemess(ddata2(i,:),form,cpu_number)
						end do
					case(3)
						call mess%pointer(ddata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(ddata3(i,:,j),form,cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(ddata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(ddata4(i,:,j,k),form,cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%di(),form,cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(4)
				select case(rank)
					case(1)
						call writemess(mess%ci(),form,cpu_number)
					case(2)
						call mess%pointer(cdata2)
						do i=1,mess%dim(1)
							call writemess(cdata2(i,:),form,cpu_number)
						end do
					case(3)
						call mess%pointer(cdata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(cdata3(i,:,j),form,cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(cdata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(cdata4(i,:,j,k),form,cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%ci(),form,cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(5)
				select case(rank)
					case(1)
						call writemess(mess%zi(),form,cpu_number)
					case(2)
						call mess%pointer(zdata2)
						do i=1,mess%dim(1)
							call writemess(zdata2(i,:),form,cpu_number)
						end do
					case(3)
						call mess%pointer(zdata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(zdata3(i,:,j),form,cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(zdata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(zdata4(i,:,j,k),form,cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%zi(),form,cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(6)
				select case(rank)
					case(1)
						call writemess(mess%li(),cpu_number)
					case(2)
						call mess%pointer(ldata2)
						do i=1,mess%dim(1)
							call writemess(ldata2(i,:),cpu_number)
						end do
					case(3)
						call mess%pointer(ldata3)
						do j=1,mess%dim(3)
							call writemess('(*,*,'+j+')') 
							do i=1,mess%dim(2)
								call writemess(ldata3(i,:,j),cpu_number)
							end do
						end do
					case(4)
						call mess%pointer(ldata4)
						do k=1,mess%dim(4)
							do j=1,mess%dim(3)
								call writemess('(*,*,'+j+','+k+')') 
								do i=1,mess%dim(2)
									call writemess(ldata4(i,:,j,k),cpu_number)
								end do
							end do
						end do
					case default
						call writemess('The Tensor Data are:',cpu_number)
						call writemess(mess%li(),cpu_number)
						call writemess(mess%Dimension,cpu_number)
					end select
			case(7)
				do i=1,totoal-1
					w=w+mess%ai(i)+','
				end do
				w=w+mess%ai(totoal)
				call writemess(w,cpu_number)
		end select
		return
	end subroutine
	
	subroutine Dprint2(Dimen)
		class(Tensor),intent(in) ::Dimen
		CHARACTER(len=20)::classTypeChar
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)""
		if(Dimen%getflag()) then
			classTypeChar=Dimen%getclassType()
			if(Dimen%ifDynamic())then
				write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(*,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(*,*) "The rank of the Tensor is"
			write(*,*) Dimen%Dimension%getRank()
			write(*,*) "The number of  data of the Tensor is"
			write(*,*) Dimen%TData%getTotalData()
			write(*,*) "The data of the Tensor is"
			call Tprintdata(Dimen%TData,0)
			call Dimen%Dimension%print()
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
	subroutine Tprint2(T,words,printType)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)trim(words)
		if(T%getflag()) then
			classTypeChar=T%getclassType()
			if(T%ifDynamic())then
				write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(*,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(*,*) "The rank of the Tensor is"
			write(*,*) T%Dimension%getRank()
			write(*,*) "The number of  data of the Tensor is"
			write(*,*) T%TData%getTotalData()
			write(*,*) "The data of the Tensor is"
			call Tprintdata(T%TData,0,printType)
			call T%Dimension%print()
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
	subroutine Tprint3(T,realpart,printType)
		class(Tensor),intent(in) :: T
		integer,intent(in)::realpart
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)""
		if(T%getflag()) then
			classTypeChar=T%getclassType()
			if(T%ifDynamic())then
				write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(*,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(*,*) "The rank of the Tensor is"
			write(*,*) T%Dimension%getRank()
			write(*,*) "The number of  data of the Tensor is"
			write(*,*) T%TData%getTotalData()
			write(*,*) "The data of the Tensor is"
			call Tprintdata(T%TData,realpart,printType)
			call T%Dimension%print()
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
	subroutine Tprint4(T,words,realpart,printType)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		CHARACTER(len=*),optional,intent(in)::printType
		integer,intent(in)::realpart
		CHARACTER(len=20)::classTypeChar
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)trim(words)
		if(T%getflag()) then
			classTypeChar=T%getclassType()
			if(T%ifDynamic())then
				write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(*,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(*,*) "The rank of the Tensor is"
			write(*,*) T%Dimension%getRank()
			write(*,*) "The number of  data of the Tensor is"
			write(*,*) T%TData%getTotalData()
			write(*,*) "The data of the Tensor is"
			call Tprintdata(T%TData,realpart,printType)
			call T%Dimension%print()
			write(*,*) "***end***"
			write(*,*) ""
		else
			write(*,*) "There is no data"
		end if
		return
	end subroutine
	subroutine Tprint_file1(T,words,uni,printType)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		integer,intent(in)::uni
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		write(uni,*)"=================="
		write(uni,*)"readable data"
		write(uni,*)trim(words)
		if(T%getflag()) then
			classTypeChar=T%getclassType()
			if(T%ifDynamic())then
				write(uni,*)'Dynamic class Tensor,data type is:'
				write(uni,*)classTypeChar
			else
				write(uni,*)'Static class Tensor,data type is:'
				write(uni,*)classTypeChar
			end if
			write(uni,*) "The rank of the Tensor is"
			write(uni,*) T%Dimension%getRank()
			write(uni,*) "The number of  data of the Tensor is"
			write(uni,*) T%TData%getTotalData()
			if((T%getType().eq.4).or.(T%getType().eq.5))then
				write(uni,*) "The data of the Tensor is(real)"
				call Tprintdata_file(T%TData,uni,1,printType)
				write(uni,*) "The data of the Tensor is(imag)"
				call Tprintdata_file(T%TData,uni,2,printType)
			else
				write(uni,*) "The data of the Tensor is"
				call Tprintdata_file(T%TData,uni,0,printType)
			end if
			call T%Dimension%info_file(uni)
			write(uni,*) "***END***"
			write(uni,*) ""
		else
			write(uni,*) "There is no data"
			write(uni,*)"END"
		end if
		return
	end subroutine
	subroutine Tprint_file2(T,uni,printType)
		class(Tensor),intent(in) :: T
		integer,intent(in)::uni
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		write(uni,*)"=================="
		write(uni,*)"readable data"
		write(uni,*)"-"
		if(T%getflag()) then
			classTypeChar=T%getclassType()
			if(T%ifDynamic())then
				write(uni,*)'Dynamic class Tensor,data type is:'
				write(uni,*)classTypeChar
			else
				write(uni,*)'Static class Tensor,data type is:'
				write(uni,*)classTypeChar
			end if
			write(uni,*) "The rank of the Tensor is"
			write(uni,*) T%Dimension%getRank()
			write(uni,*) "The number of  data of the Tensor is"
			write(uni,*) T%TData%getTotalData()
			if((T%getType().eq.4).or.(T%getType().eq.5))then
				write(uni,*) "The data of the Tensor is(real)"
				call Tprintdata_file(T%TData,uni,1,printType)
				write(uni,*) "The data of the Tensor is(imag)"
				call Tprintdata_file(T%TData,uni,2,printType)
			else
				write(uni,*) "The data of the Tensor is"
				call Tprintdata_file(T%TData,uni,0,printType)
			end if
			call T%Dimension%info_file(uni)
			write(uni,*) "***END***"
			write(uni,*) ""
		else
			write(uni,*) "There is no data"
			write(uni,*)"END"
		end if
		return
	end subroutine

	subroutine Tread_file(T,uni)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::uni
		CHARACTER(len=20)::classTypeChar
		CHARACTER(len=50)::notused
		integer::rank,TotalData,i
		integer,allocatable::idata(:)
		real*4,allocatable::sdata(:),sdata2(:)
		real*8,allocatable::ddata(:),ddata2(:)
		character(len=max_len_of_char_in_TData),allocatable::adata(:)
		logical,allocatable::ldata(:)
		read(uni,*)notused
		read(uni,*)notused
		if(notused.ne.'readable')then
			call writemess("error in reading",-1)
			call error_stop()
		end if
		read(uni,*)notused
		read(uni,*)notused
		read(uni,*)classTypeChar
		call T%empty()
		if(classTypeChar.equ.'END')return
		call T%setType(classTypeChar)
		read(uni,*) notused
		read(uni,*) rank
		read(uni,*) notused
		read(uni,*) TotalData
		select case(T%getType())
			case(1)
				allocate(idata(TotalData))
				read(uni,*) notused
				read(uni,*)(idata(i),i=1,TotalData)
				T=idata
			case(2)
				allocate(sdata(TotalData))
				read(uni,*) notused
				read(uni,*)(sdata(i),i=1,TotalData)
				T=sdata
			case(3)
				allocate(ddata(TotalData))
				read(uni,*) notused
				read(uni,*)(ddata(i),i=1,TotalData)
				T=ddata
			case(4)
				allocate(sdata(TotalData))
				allocate(sdata2(TotalData))
				read(uni,*) notused
				read(uni,*)(sdata(i),i=1,TotalData)
				read(uni,*) notused
				read(uni,*)(sdata2(i),i=1,TotalData)
				T=cmplx(sdata,sdata2,kind=4)
			case(5)
				allocate(ddata(TotalData))
				allocate(ddata2(TotalData))
				read(uni,*) notused
				read(uni,*)(ddata(i),i=1,TotalData)
				read(uni,*) notused
				read(uni,*)(ddata2(i),i=1,TotalData)
				T=dcmplx(ddata,ddata2)
			case(6)
				allocate(ldata(TotalData))
				read(uni,*) notused
				read(uni,*)(ldata(i),i=1,TotalData)
				T=ldata
			case(7)
				allocate(adata(TotalData))
				read(uni,*) notused
				read(uni,*)(adata(i),i=1,TotalData)
				T=adata
		end select
		call T%Dimension%read(uni)
		read(uni,*) notused
		return
	end subroutine

	subroutine Tread_data(T,uni)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::uni
		CHARACTER(len=50)::notused
		integer::rank,TotalData,i,dim1,dim2,j
		integer,pointer::idata(:,:),iidata(:)
		real*4,pointer::sdata(:,:),ssdata(:)
		real*8,pointer::ddata(:,:),dddata(:)
		character(len=max_len_of_char_in_TData),pointer::adata(:,:),aadata(:)
		logical,pointer::ldata(:,:),lldata(:)
		if(T%Dimension%getRank().gt.2)then
			call writemess('ERROR in reading data, only allow for rank<=2 Tensor',-1)
			call error_stop
		end if
		if(T%Dimension%getRank().eq.1)then
			TotalData=T%TData%getTotalData()
			select case(T%getType())
				case(1)
					call T%pointer(iidata)
					read(uni,*)(iidata(i),i=1,TotalData)
					nullify(iidata)
				case(2)
					call T%pointer(ssdata)
					read(uni,*)(ssdata(i),i=1,TotalData)
					nullify(ssdata)
				case(3)
					call T%pointer(dddata)
					read(uni,*)(dddata(i),i=1,TotalData)
					nullify(dddata)
				case(4)
					call writemess('ERROR in reading data, Tensor',-1)
					call writemess('Do not finished for complex data yet',-1)
					call writemess('You can real two real data and combine them into a complex one',-1)
					call writemess('for example:A and B are real*4 Tensor.',-1)
					call writemess(' call A%readData(unit1).',-1)
					call writemess(' call B%readData(unit2).',-1)
					call writemess(' C=cmplex(A,B).',-1)
					call error_stop
				case(5)
					call writemess('ERROR in reading data, Tensor',-1)
					call writemess('Do not finished for complex data yet',-1)
					call writemess('You can real two real data and combine them into a complex one',-1)
					call writemess('for example:A and B are real*8 Tensor.',-1)
					call writemess(' call A%readData(unit1).',-1)
					call writemess(' call B%readData(unit2).',-1)
					call writemess(' C=dcmplex(A,B).',-1)
					call error_stop
				case(6)
					call T%pointer(lldata)
					read(uni,*)(lldata(i),i=1,TotalData)
					nullify(lldata)
				case(7)
					call T%pointer(aadata)
					read(uni,*)(aadata(i),i=1,TotalData)
					nullify(aadata)
			end select
			return
		endif
		dim1=T%dim(1)
		dim2=T%dim(2)
		TotalData=T%TData%getTotalData()
		select case(T%getType())
			case(1)
				call T%pointer(idata)
				do i=1,dim1
					read(uni,*)(idata(i,j),j=1,dim2)
				end do
				nullify(idata)
			case(2)
				call T%pointer(sdata)
				do i=1,dim1
					read(uni,*)(sdata(i,j),j=1,dim2)
				end do
				nullify(sdata)
			case(3)
				call T%pointer(ddata)
				do i=1,dim1
					read(uni,*)(ddata(i,j),j=1,dim2)
				end do
				nullify(ddata)
			case(4)
				call writemess('ERROR in reading data, Tensor',-1)
				call writemess('Do not finished for complex data yet',-1)
				call writemess('You can real two real data and combine them into a complex one',-1)
				call writemess('for example:A and B are real*4 Tensor.',-1)
				call writemess(' call A%readData(unit1).',-1)
				call writemess(' call B%readData(unit2).',-1)
				call writemess(' C=cmplex(A,B).',-1)
				call error_stop
			case(5)
				call writemess('ERROR in reading data, Tensor',-1)
				call writemess('Do not finished for complex data yet',-1)
				call writemess('You can real two real data and combine them into a complex one',-1)
				call writemess('for example:A and B are real*8 Tensor.',-1)
				call writemess(' call A%readData(unit1).',-1)
				call writemess(' call B%readData(unit2).',-1)
				call writemess(' C=dcmplex(A,B).',-1)
				call error_stop
			case(6)
				call T%pointer(ldata)
					do i=1,dim1
						read(uni,*)(ldata(i,j),j=1,dim2)
					end do
					nullify(ldata)
			case(7)
				call T%pointer(adata)
				do i=1,dim1
					read(uni,*)(adata(i,j),j=1,dim2)
				end do
				nullify(adata)
		end select
		return
	end subroutine

	subroutine Dprint(Dimen)
		class(Tensor),intent(in) ::Dimen
		CHARACTER(len=20)::classTypeChar
		integer,allocatable::dimenen(:)
		write(*,*)"=================="
		write(*,*)"------------------"
		classTypeChar=Dimen%getclassType()
		if(Dimen%ifDynamic())then
			write(*,*)'Dynamic,',classTypeChar
		else
			write(*,*)'Static,',classTypeChar
		end if
		write(*,*) "*** START ***"
		if(.not.Dimen%getFlag())then
			write(*,*)"There is no data"
			write(*,*) "*** END ***"
			return
		end if
		
		select case(Dimen%Dimension%getRank())
			case(1)
				call Tprintdata(Dimen%TData,0)
				write(*,*) "*** END ***"
			case(2)
				allocate(dimenen(2))
				dimenen=Dimen%dim()
				call Tprint_as_matrix(Dimen%TData,0,Dimen%Dimension%getRank(),dimenen)
				write(*,*) "*** END ***"
			case(3)
				allocate(dimenen(3))
				dimenen=Dimen%dim()
				call Tprint_as_matrix(Dimen%TData,0,Dimen%Dimension%getRank(),dimenen)
				write(*,*) "*** END ***"
			case(4)
				allocate(dimenen(4))
				dimenen=Dimen%dim()
				call Tprint_as_matrix(Dimen%TData,0,Dimen%Dimension%getRank(),dimenen)
				write(*,*) "*** END ***"
			case default
				write(*,*) "rank of the Tensor is large than 4"
				classTypeChar=Dimen%getclassType()
				if(Dimen%ifDynamic())then
					write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
				else
					write(*,*)'static class Tensor,data type is,',classTypeChar
				end if
				write(*,*) "The data of the Tensor is"
				call Tprintdata(Dimen%TData,0)
				write(*,*) "***end***"
				write(*,*) "The dimension of the Tensor is"
				call Dimen%Dimension%print()
				write(*,*) "The rank,total data are"
				write(*,*) Dimen%Dimension%getRank(),Dimen%TData%getTotalData()
		end select
		return
	end subroutine
	subroutine TMprint2(T,words,printType)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		integer,allocatable::dimen(:)
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)trim(words)
		classTypeChar=T%getclassType()
		if(T%ifDynamic())then
			write(*,*)'Dynamic,',classTypeChar
		else
			write(*,*)'Static,',classTypeChar
		end if
		write(*,*) "*** START ***"
		if(.not.T%getFlag())then
			write(*,*)"There is no data"
			write(*,*) "*** END ***"
			return
		end if
		
		select case(T%Dimension%getRank())
			case(1)
				call Tprintdata(T%TData,0,printType)
				write(*,*) "*** END ***"
			case(2)
				allocate(dimen(2))
				dimen=T%dim()
				call Tprint_as_matrix(T%TData,0,T%Dimension%getRank(),dimen,printType)
				write(*,*) "*** END ***"
			case(3)
				allocate(dimen(3))
				dimen=T%dim()
				call Tprint_as_matrix(T%TData,0,T%Dimension%getRank(),dimen,printType)
				write(*,*) "*** END ***"
			case(4)
				allocate(dimen(4))
				dimen=T%dim()
				call Tprint_as_matrix(T%TData,0,T%Dimension%getRank(),dimen,printType)
				write(*,*) "*** END ***"
			case default
				write(*,*) "rank of the Tensor is large than 4"
				classTypeChar=T%getclassType()
				if(T%ifDynamic())then
					write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
				else
					write(*,*)'static class Tensor,data type is,',classTypeChar
				end if
				write(*,*) "The data of the Tensor is"
				call Tprintdata(T%TData,0,printType)
				write(*,*) "***end***"
				write(*,*) "The dimension of the Tensor is"
				call T%Dimension%print()
				write(*,*) "The rank,total data are"
				write(*,*) T%Dimension%getRank(),T%TData%getTotalData()
		end select
		return
	end subroutine
	subroutine TMprint3(T,realpart,printType)
		class(Tensor),intent(in) :: T
		integer,intent(in)::realpart
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		integer,allocatable::dimen(:)
		write(*,*)"=================="
		write(*,*)"------------------"
		classTypeChar=T%getclassType()
		if(T%ifDynamic())then
			write(*,*)'Dynamic,',classTypeChar
		else
			write(*,*)'Static,',classTypeChar
		end if
		write(*,*) "*** START ***"
		if(.not.T%getFlag())then
			write(*,*)"There is no data"
			write(*,*) "*** END ***"
			return
		end if
		select case(T%Dimension%getRank())
			case(1)
				call Tprintdata(T%TData,realpart,printType)
				write(*,*) "*** END ***"
			case(2)
				allocate(dimen(2))
				dimen=T%dim()
				call Tprint_as_matrix(T%TData,realpart,T%Dimension%getRank(),dimen,printType)
				write(*,*) "*** END ***"
			case(3)
				allocate(dimen(3))
				dimen=T%dim()
				call Tprint_as_matrix(T%TData,realpart,T%Dimension%getRank(),dimen,printType)
				write(*,*) "*** END ***"
			case(4)
				allocate(dimen(4))
				dimen=T%dim()
				call Tprint_as_matrix(T%TData,realpart,T%Dimension%getRank(),dimen,printType)
				write(*,*) "*** END ***"
			case default
				write(*,*) "rank of the Tensor is large than 4"
				classTypeChar=T%getclassType()
				if(T%ifDynamic())then
					write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
				else
					write(*,*)'static class Tensor,data type is,',classTypeChar
				end if
				write(*,*) "The data of the Tensor is"
				call Tprintdata(T%TData,realpart,printType)
				write(*,*) "***end***"
				write(*,*) "The dimension of the Tensor is"
				call T%Dimension%print()
				write(*,*) "The rank,total data are"
				write(*,*) T%Dimension%getRank(),T%TData%getTotalData()
		end select
		return
	end subroutine
	subroutine TMprint4(T,words,realpart,printType)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		integer,intent(in)::realpart
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		integer,allocatable::dimen(:)
		write(*,*)"=================="
		write(*,*)"------------------"
		write(*,*)trim(words)
		classTypeChar=T%getclassType()
		if(T%ifDynamic())then
			write(*,*)'Dynamic,',classTypeChar
		else
			write(*,*)'Static,',classTypeChar
		end if
		write(*,*) "*** START ***"
		if(.not.T%getFlag())then
			write(*,*)"There is no data"
			write(*,*) "*** END ***"
			return
		end if
		select case(T%Dimension%getRank())
			case(1)
				call Tprintdata(T%TData,realpart,printType)
				write(*,*) "*** END ***"
			case(2)
				allocate(dimen(2))
				dimen=T%dim()
				call Tprint_as_matrix(T%TData,realpart,T%Dimension%getRank(),dimen,printType)
				write(*,*) "*** END ***"
			case(3)
				allocate(dimen(3))
				dimen=T%dim()
				call Tprint_as_matrix(T%TData,realpart,T%Dimension%getRank(),dimen,printType)
				write(*,*) "*** END ***"
			case(4)
				allocate(dimen(4))
				dimen=T%dim()
				call Tprint_as_matrix(T%TData,realpart,T%Dimension%getRank(),dimen,printType)
				write(*,*) "*** END ***"
			case default
				write(*,*) "rank of the Tensor is large than 4"
				classTypeChar=T%getclassType()
				if(T%ifDynamic())then
					write(*,*)'Dynamic class Tensor,data type is,',classTypeChar
				else
					write(*,*)'static class Tensor,data type is,',classTypeChar
				end if
				write(*,*) "The data of the Tensor is"
				call Tprintdata(T%TData,realpart,printType)
				write(*,*) "***end***"
				write(*,*) "The dimension of the Tensor is"
				call T%Dimension%print()
				write(*,*) "The rank,total data are"
				write(*,*) T%Dimension%getRank(),T%TData%getTotalData()
		end select
		return
	end subroutine

	subroutine TMprint_file1(T,words,uni,printType)
		class(Tensor),intent(in) :: T
		integer,intent(in)::uni
		CHARACTER(len=*),intent(in)::words
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		integer,allocatable::dimen(:)
		write(uni,*)"=================="
		write(uni,*)"------------------"
		write(uni,*)trim(words)
		classTypeChar=T%getclassType()
		if(T%ifDynamic())then
				write(uni,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(uni,*)'static class Tensor,data type is,',classTypeChar
			end if
		write(uni,*) "*** START ***"
		if(.not.T%getFlag())then
			write(uni,*)"There is no data"
			write(uni,*) "*** END ***"
			return
		end if
		select case(T%Dimension%getRank())
			case(1)
				call Tprintdata_file(T%TData,uni,0,printType)
				write(uni,*) "*** END ***"
			case(2)
				allocate(dimen(2))
				dimen=T%dim()
				call Tprint_as_matrix_file(T%TData,uni,0,T%Dimension%getRank(),dimen,printType)
				write(uni,*) "*** END ***"
			case(3)
				allocate(dimen(3))
				dimen=T%dim()
				call Tprint_as_matrix_file(T%TData,uni,0,T%Dimension%getRank(),dimen,printType)
				write(uni,*) "*** END ***"
			case(4)
				allocate(dimen(4))
				dimen=T%dim()
				call Tprint_as_matrix_file(T%TData,uni,0,T%Dimension%getRank(),dimen,printType)
				write(uni,*) "*** END ***"
			case default
				write(uni,*) "rank of the Tensor is large than 4"
				classTypeChar=T%getclassType()
				if(T%ifDynamic())then
					write(uni,*)'Dynamic class Tensor,data type is,',classTypeChar
				else
					write(uni,*)'static class Tensor,data type is,',classTypeChar
				end if
				write(uni,*) "The data of the Tensor is"
				call Tprintdata_file(T%TData,uni,0,printType)
				write(uni,*) "***end***"
				write(uni,*) "The dimension of the Tensor is"
				call T%Dimension%print_file(uni)
				write(uni,*) "The rank,total data are"
				write(uni,*) T%Dimension%getRank(),T%TData%getTotalData()
		end select
		return
	end subroutine
	subroutine TMprint_file2(T,uni,printType)
		class(Tensor),intent(in) :: T
		integer,intent(in)::uni
		CHARACTER(len=*),optional,intent(in)::printType
		CHARACTER(len=20)::classTypeChar
		integer,allocatable::dimen(:)
		write(uni,*)"=================="
		write(uni,*)"------------------"
		classTypeChar=T%getclassType()
		if(T%ifDynamic())then
				write(uni,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(uni,*)'static class Tensor,data type is,',classTypeChar
			end if
		write(uni,*) "*** START ***"
		if(.not.T%getFlag())then
			write(uni,*)"There is no data"
			write(uni,*) "*** END ***"
			return
		end if
		select case(T%Dimension%getRank())
			case(1)
				call Tprintdata_file(T%TData,uni,0,printType)
				write(uni,*) "*** END ***"
			case(2)
				allocate(dimen(2))
				dimen=T%dim()
				call Tprint_as_matrix_file(T%TData,uni,0,T%Dimension%getRank(),dimen,printType)
				write(uni,*) "*** END ***"
			case(3)
				allocate(dimen(3))
				dimen=T%dim()
				call Tprint_as_matrix_file(T%TData,uni,0,T%Dimension%getRank(),dimen,printType)
				write(uni,*) "*** END ***"
			case(4)
				allocate(dimen(4))
				dimen=T%dim()
				call Tprint_as_matrix_file(T%TData,uni,0,T%Dimension%getRank(),dimen,printType)
				write(uni,*) "*** END ***"
			case default
				write(uni,*) "rank of the Tensor is large than 4"
				classTypeChar=T%getclassType()
				if(T%ifDynamic())then
					write(uni,*)'Dynamic class Tensor,data type is,',classTypeChar
				else
					write(uni,*)'static class Tensor,data type is,',classTypeChar
				end if
				write(uni,*) "The data of the Tensor is"
				call Tprintdata_file(T%TData,uni,0,printType)
				write(uni,*) "***end***"
				write(uni,*) "The dimension of the Tensor is"
				call T%dimension%print_file(uni)
				write(uni,*) "The rank,total data are"
				write(uni,*) T%Dimension%getRank(),T%TData%getTotalData()
		end select
		return
	end subroutine

	subroutine readdimension(dimen,uni)
		class(Tensor),target,intent(inout)::dimen
		integer,intent(in)::uni
		class(Tensor),pointer::T
		CHARACTER(len=20)::classTypeChar
		CHARACTER(len=50)::notused
		integer::rank,TotalData,i
		integer,allocatable::idata(:)
		real*4,allocatable::sdata(:),sdata2(:)
		real*8,allocatable::ddata(:),ddata2(:)
		character(len=max_len_of_char_in_TData),allocatable::adata(:)
		logical,allocatable::ldata(:)
		T=>Dimen
		read(uni,*)notused
		read(uni,*)notused
		if(notused.ne.'readable')then
			call writemess("error in reading",-1)
			call error_stop()
		end if
		read(uni,*)notused
		read(uni,*)notused
		read(uni,*)classTypeChar
		call T%empty()
		if(classTypeChar.equ.'END')return
		call T%setType(classTypeChar)
		read(uni,*) notused
		read(uni,*) rank
		read(uni,*) notused
		read(uni,*) TotalData
		select case(T%getType())
			case(1)
				allocate(idata(TotalData))
				read(uni,*) notused
				read(uni,*)(idata(i),i=1,TotalData)
				T=idata
			case(2)
				allocate(sdata(TotalData))
				read(uni,*) notused
				read(uni,*)(sdata(i),i=1,TotalData)
				T=sdata
			case(3)
				allocate(ddata(TotalData))
				read(uni,*) notused
				read(uni,*)(ddata(i),i=1,TotalData)
				T=ddata
			case(4)
				allocate(sdata(TotalData))
				allocate(sdata2(TotalData))
				read(uni,*) notused
				read(uni,*)(sdata(i),i=1,TotalData)
				read(uni,*) notused
				read(uni,*)(sdata2(i),i=1,TotalData)
				T=cmplx(sdata,sdata2,kind=4)
			case(5)
				allocate(ddata(TotalData))
				allocate(ddata2(TotalData))
				read(uni,*) notused
				read(uni,*)(ddata(i),i=1,TotalData)
				read(uni,*) notused
				read(uni,*)(ddata2(i),i=1,TotalData)
				T=dcmplx(ddata,ddata2)
			case(6)
				allocate(ldata(TotalData))
				read(uni,*) notused
				read(uni,*)(ldata(i),i=1,TotalData)
				T=ldata
			case(7)
				allocate(adata(TotalData))
				read(uni,*) notused
				read(uni,*)(adata(i),i=1,TotalData)
				T=adata
		end select
		call T%Dimension%read(uni)
		read(uni,*) notused
		return
	end subroutine

	!********************* print dimension *********************    

	subroutine TDprint1(T,words,uni)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::words
		integer,intent(in)::uni
		CHARACTER(len=20)::classTypeChar
		write(uni,*)"=================="
		write(uni,*)"------------------"
		write(uni,*)trim(words)
		if(T%getflag()) then!if1
			write(uni,*) "*** START ***"
			classTypeChar=T%getclassType()
			if(T%ifDynamic())then
				write(uni,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(uni,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(uni,*) "The rank of the Tensor is"
			write(uni,*) T%Dimension%getRank()
			write(uni,*) "The number of  data of the Tensor is"
			write(uni,*) T%TData%getTotalData()
			call T%Dimension%print_file(uni)
			write(uni,*) "***end***"
			write(uni,*) ""
		else!if1
			write(uni,*) "There is no data"
		end if!if1
		return
	end subroutine
	subroutine TDprint2(T,words)
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),optional,intent(in)::words
		CHARACTER(len=20)::classTypeChar
		call writemess("==================",-1)
		call writemess("------------------",-1)
		call writemess(words)
		if(T%getflag()) then!if1
			call writemess("*** START ***",-1)
			classTypeChar=T%getclassType()
			if(T%ifDynamic())then
				call writemess('Dynamic class Tensor,data type is '+(' '+classTypeChar),-1)
			else
				call writemess('static class Tensor,data type is '+(' '+classTypeChar),-1)
			end if
			call writemess("The rank of the Tensor is",-1)
			call writemess(''+T%Dimension%getRank())
			call writemess("The number of  data of the Tensor is",-1)
			call writemess(''+T%TData%getTotalData(),-1)
			call T%Dimension%print()
			call writemess( "***end***",-1)
			call writemess("",-1)
		else!if1
			call writemess("There is no data",-1)
		end if!if1
		return
	end subroutine
	subroutine TDprint3(T,uni)
		class(Tensor),intent(in) :: T
		integer,intent(in)::uni
		CHARACTER(len=20)::classTypeChar
		write(uni,*)"=================="
		write(uni,*)"------------------"
		if(T%getflag()) then!if1
			write(uni,*) "*** START ***"
			classTypeChar=T%getclassType()
			if(T%ifDynamic())then
				write(uni,*)'Dynamic class Tensor,data type is,',classTypeChar
			else
				write(uni,*)'static class Tensor,data type is,',classTypeChar
			end if
			write(uni,*) "The rank of the Tensor is"
			write(uni,*) T%Dimension%getRank()
			write(uni,*) "The number of  data of the Tensor is"
			write(uni,*) T%TData%getTotalData()
			call T%Dimension%print_file(uni)
			write(uni,*) "***end***"
			write(uni,*) ""
		else!if1
			write(uni,*) "There is no data"
		end if!if1
		return
	end subroutine
	
	!******************************************************************************
	!******************************************************************************
	!
	!               output elements    
	!
	!******************************************************************************
	!******************************************************************************

	integer function addressToIndes(T,Adim)
		class(Tensor),intent(in) :: T
		integer,intent(in) :: Adim(:)
		integer,allocatable::Tdim(:)
		integer :: i,Dimlen
		call copydimension(Tdim,T%Dimension)
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
	integer function iElement(T,TDim)
		class(Tensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		integer::inde
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(T%Dimension%getRank().eq.1)then
			if(Tdim(1).gt.T%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_int(iElement,T%TData,Tdim,(/T%dim(1)/),1)
			return
		end if
		if(T%Dimension%getRank().eq.2)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_int(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2)/),2)
			return
		end if
		if(T%Dimension%getRank().eq.3)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt. (T%dim(2))).or. (Tdim(3).gt.T%dim(3)) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_int(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3)/),3)
			return
		end if
		if(T%Dimension%getRank().eq.4)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))).or.(Tdim(3).gt.(T%dim(3))).or.(Tdim(4).gt.(T%dim(4))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_int(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3),T%dim(4)/),4)
			return
		end if
		inde=addressToIndes(T,Tdim)
		if(inde.gt.T%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
			call error_stop()
		end if
		call element_routine_int(iElement,T%TData,(/inde/),(/T%TData%getTotalData()/),1)
		return
	end function		  
	function sElement(T,Tdim)result(iElement)
		real(kind=4) ::iElement
		class(Tensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		integer::inde
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(T%Dimension%getRank().eq.1)then
			if(Tdim(1).gt.T%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_real4(iElement,T%TData,Tdim,(/T%dim(1)/),1)
			return
		end if
		if(T%Dimension%getRank().eq.2)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_real4(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2)/),2)
			return
		end if
		if(T%Dimension%getRank().eq.3)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt. (T%dim(2))).or. (Tdim(3).gt.(T%dim(3))))Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_real4(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3)/),3)
			return
		end if
		if(T%Dimension%getRank().eq.4)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))).or.(Tdim(3).gt.(T%dim(3))).or.(Tdim(4).gt.(T%dim(4))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_real4(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3),T%dim(4)/),4)
			return
		end if
		inde=addressToIndes(T,Tdim)
		if(inde.gt.T%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
			call error_stop()
		end if
		call element_routine_real4(iElement,T%TData,(/inde/),(/T%TData%getTotalData()/),1)
		return
	end function
	function dElement(T,Tdim)result(iElement)
		real(kind=8) ::iElement
		class(Tensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		integer::inde
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(T%Dimension%getRank().eq.1)then
			if(Tdim(1).gt.T%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_real8(iElement,T%TData,Tdim,(/T%dim(1)/),1)
			return
		end if
		if(T%Dimension%getRank().eq.2)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_real8(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2)/),2)
			return
		end if
		if(T%Dimension%getRank().eq.3)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt. (T%dim(2))).or. (Tdim(3).gt.(T%dim(3))))Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_real8(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3)/),3)
			return
		end if
		if(T%Dimension%getRank().eq.4)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))).or.(Tdim(3).gt.(T%dim(3))).or.(Tdim(4).gt.(T%dim(4))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_real8(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3),T%dim(4)/),4)
			return
		end if
		inde=addressToIndes(T,Tdim)
		if(inde.gt.T%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
			call error_stop()
		end if
		call element_routine_real8(iElement,T%TData,(/inde/),(/T%TData%getTotalData()/),1)
		return
	end function		
	function cElement(T,Tdim)result(iElement)
		complex(kind=4) ::iElement
		class(Tensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		integer::inde
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(T%Dimension%getRank().eq.1)then
			if(Tdim(1).gt.T%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_com4(iElement,T%TData,Tdim,(/T%dim(1)/),1)
			return
		end if
		if(T%Dimension%getRank().eq.2)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_com4(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2)/),2)
			return
		end if
		if(T%Dimension%getRank().eq.3)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt. (T%dim(2))).or. (Tdim(3).gt.(T%dim(3))))Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_com4(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3)/),3)
			return
		end if
		if(T%Dimension%getRank().eq.4)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))).or.(Tdim(3).gt.(T%dim(3))).or.(Tdim(4).gt.(T%dim(4))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_com4(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3),T%dim(4)/),4)
			return
		end if
		inde=addressToIndes(T,Tdim)
		if(inde.gt.T%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
			call error_stop()
		end if
		call element_routine_com4(iElement,T%TData,(/inde/),(/T%TData%getTotalData()/),1)
		return
	end function  
	function zElement(T,Tdim)result(iElement)
		complex(kind=8) ::iElement
		class(Tensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		integer::inde
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(T%Dimension%getRank().eq.1)then
			if(Tdim(1).gt.T%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call writemess("you have input:"+Tdim(1))
				call error_stop()
			end if
			call element_routine_com8(iElement,T%TData,Tdim,(/T%dim(1)/),1)
			return
		end if
		if(T%Dimension%getRank().eq.2)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call writemess("you have input:("+Tdim(1)+","+Tdim(2)+")")
				call writemess("dimension of T is input:("+(T%dim(1))+","+(T%dim(2))+")")
				call error_stop()
			end if
			call element_routine_com8(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2)/),2)
			return
		end if
		if(T%Dimension%getRank().eq.3)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt. (T%dim(2))).or. (Tdim(3).gt.(T%dim(3))))Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_com8(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3)/),3)
			return
		end if
		if(T%Dimension%getRank().eq.4)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))).or.(Tdim(3).gt.(T%dim(3))).or.(Tdim(4).gt.(T%dim(4))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_com8(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3),T%dim(4)/),4)
			return
		end if
		inde=addressToIndes(T,Tdim)
		if(inde.gt.T%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
			call error_stop()
		end if
		call element_routine_com8(iElement,T%TData,(/inde/),(/T%TData%getTotalData()/),1)
		return
	end function  
	logical function lElement(T,Tdim)result(iElement)
		class(Tensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		integer::inde
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(T%Dimension%getRank().eq.1)then
			if(Tdim(1).gt.T%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_logi(iElement,T%TData,Tdim,(/T%dim(1)/),1)
			return
		end if
		if(T%Dimension%getRank().eq.2)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_logi(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2)/),2)
			return
		end if
		if(T%Dimension%getRank().eq.3)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt. (T%dim(2))).or. (Tdim(3).gt.(T%dim(3))))Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_logi(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3)/),3)
			return
		end if
		if(T%Dimension%getRank().eq.4)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))).or.(Tdim(3).gt.(T%dim(3))).or.(Tdim(4).gt.(T%dim(4))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_logi(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3),T%dim(4)/),4)
			return
		end if
		inde=addressToIndes(T,Tdim)
		if(inde.gt.T%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
			call error_stop()
		end if
		call element_routine_logi(iElement,T%TData,(/inde/),(/T%TData%getTotalData()/),1)
		return
	end function  

	function aElement(T,Tdim)result(iElement)
		character(len=max_len_of_char_in_TData) ::iElement
		class(Tensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		integer::inde
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(T%Dimension%getRank().eq.1)then
			if(Tdim(1).gt.T%TData%getTotalData())Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_char(iElement,T%TData,Tdim,(/T%dim(1)/),1)
			return
		end if
		if(T%Dimension%getRank().eq.2)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_char(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2)/),2)
			return
		end if
		if(T%Dimension%getRank().eq.3)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt. (T%dim(2))).or. (Tdim(3).gt.(T%dim(3))))Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_char(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3)/),3)
			return
		end if
		if(T%Dimension%getRank().eq.4)then
			if( (Tdim(1).gt.(T%dim(1))) .or. (Tdim(2).gt.(T%dim(2))).or.(Tdim(3).gt.(T%dim(3))).or.(Tdim(4).gt.(T%dim(4))) )Then
				call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
				call error_stop()
			end if
			call element_routine_char(iElement,T%TData,Tdim,(/T%dim(1),T%dim(2),T%dim(3),T%dim(4)/),4)
			return
		end if
		inde=addressToIndes(T,Tdim)
		if(inde.gt.T%TData%getTotalData())Then
			call writemess("Index is larger than the len of Tensor,totalData="+T%TData%getTotalData(),-1)
			call error_stop()
		end if
		call element_routine_char(iElement,T%TData,(/inde/),(/T%TData%getTotalData()/),1)
		return
	end function  
	function TElement(T,Tdim)result(Element)
		class(Tensor),allocatable::Element
		class(Tensor),intent(in) ::T
		integer,intent(in)::Tdim(:)
		allocate(Element,mold=T)
		select case(T%getType())
			case (1)
				Element=T%ii(Tdim)
			case (2)
				Element=T%si(Tdim)
			case (3)
				Element=T%di(Tdim)	
			case (4)
				Element=T%ci(Tdim)
			case (5)
				Element=T%zi(Tdim)
			case (6)
				Element=T%li(Tdim)
			case (7)
				Element=T%ai(Tdim)
			case default
				call writemess("ERROR in TElement",-1)
				call error_stop()
		end select
		return
	end function
	function TElement2(T,Tdim)result(Element)
		class(Tensor),allocatable::Element
		class(Tensor),intent(in) ::T
		integer,intent(in)::Tdim
		allocate(Element,mold=T)
		select case(T%getType())
			case (1)
				Element=T%ii(Tdim)
			case (2)
				Element=T%si(Tdim)
			case (3)
				Element=T%di(Tdim)	
			case (4)
				Element=T%ci(Tdim)
			case (5)
				Element=T%zi(Tdim)
			case (6)
				Element=T%li(Tdim)
			case (7)
				Element=T%ai(Tdim)
			case default
				call writemess("ERROR in TElement",-1)
				call error_stop()
		end select
		return
	end function
	integer function iElement2 (T,inde)result(Element)
		class(Tensor),intent(in) ::T
		integer,intent(in)::inde
		integer::length
		length=T%TData%getTotalData()
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(inde.gt.length)then
			write(*,*)"index is larger then the length of the Tensor,(.i.)"
			call error_stop()
		end if
		call element_routine_int(Element,T%TData,(/inde/),(/length/),1)
		return
	end function
	function sElement2 (T,inde)result(Element)
		real(kind=4) ::Element
		class(Tensor),intent(in) ::T
		integer,intent(in)::inde
		integer::length
		length=T%TData%getTotalData()
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(inde.gt.length)then
			write(*,*)"index is larger then the length of the Tensor,(.i.)"
			call error_stop()
		end if
		call element_routine_real4(Element,T%TData,(/inde/),(/length/),1)
		return
	end function
	function dElement2 (T,inde)result(Element)
		real(kind=8) ::Element
		class(Tensor),intent(in) ::T
		integer,intent(in)::inde
		integer::length
		length=T%TData%getTotalData()
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(inde.gt.length)then
			write(*,*)"index is larger then the length of the Tensor,(.i.)"
			call error_stop()
		end if
		call element_routine_real8(Element,T%TData,(/inde/),(/length/),1)
		return
	end function
	function cElement2 (T,inde)result(Element)
		complex(kind=4) ::Element
		class(Tensor),intent(in) ::T
		integer,intent(in)::inde
		integer::length
		length=T%TData%getTotalData()
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(inde.gt.length)then
			write(*,*)"index is larger then the length of the Tensor,(.i.)"
			call error_stop()
		end if
		call element_routine_com4(Element,T%TData,(/inde/),(/length/),1)
		return
	end function
	function zElement2 (T,inde)result(Element)
		complex(kind=8) ::Element
		class(Tensor),intent(in) ::T
		integer,intent(in)::inde
		integer::length
		length=T%TData%getTotalData()
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(inde.gt.length)then
			write(*,*)"index is larger then the length of the Tensor,(.i.)"
			call error_stop()
		end if
		call element_routine_com8(Element,T%TData,(/inde/),(/length/),1)
		return
	end function
	logical function lElement2 (T,inde)result(Element)
		class(Tensor),intent(in) ::T
		integer,intent(in)::inde
		integer::length
		length=T%TData%getTotalData()
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(inde.gt.length)then
			write(*,*)"index is larger then the length of the Tensor,(.i.)"
			call error_stop()
		end if
		call element_routine_logi(Element,T%TData,(/inde/),(/length/),1)
		return
	end function
	function aElement2 (T,inde)result(Element)
		character(len=max_len_of_char_in_TData) ::Element
		class(Tensor),intent(in) ::T
		integer,intent(in)::inde
		integer::length
		length=T%TData%getTotalData()
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(.i.)"
			call error_stop()
		end if
		if(inde.gt.length)then
			write(*,*)"index is larger then the length of the Tensor,(.i.)"
			call error_stop()
		end if
		call element_routine_char(Element,T%TData,(/inde/),(/length/),1)
		return
	end function
	function ielementAll(T)result(Res)
		integer,allocatable::Res(:)
		class(Tensor),intent(in) ::T
		allocate(Res(T%TData%getTotalData()))
		Res=T
		return
	end function
	function selementAll(T)result(Res)
		real*4,allocatable::Res(:)
		class(Tensor),intent(in) ::T
		allocate(Res(T%TData%getTotalData()))
		Res=T
		return
	end function
	function delementAll(T)result(Res)
		real*8,allocatable::Res(:)
		class(Tensor),intent(in) ::T
		allocate(Res(T%TData%getTotalData()))
		Res=T
		return
	end function
	function celementAll(T)result(Res)
		complex*8,allocatable::Res(:)
		class(Tensor),intent(in) ::T
		allocate(Res(T%TData%getTotalData()))
		Res=T
		return
	end function
	function zelementAll(T)result(Res)
		complex*16,allocatable::Res(:)
		class(Tensor),intent(in) ::T
		allocate(Res(T%TData%getTotalData()))
		Res=T
		return
	end function
	function lelementAll(T)result(Res)
		logical,allocatable::Res(:)
		class(Tensor),intent(in) ::T
		allocate(Res(T%TData%getTotalData()))
		Res=T
		return
	end function
	function aelementAll(T)result(Res)
		character(len=max_len_of_char_in_TData),allocatable::Res(:)
		class(Tensor),intent(in) ::T
		allocate(Res(T%TData%getTotalData()))
		Res=T
		return
	end function

	!**************************************************************************
	!**************************************************************************
	!
	!                                  + - * /
	!
	!**************************************************************************
	!**************************************************************************

	!
		! int + int --->int
		! int + real4 --->rea4
		! int + real8 --->rea8
		! int + complex(kin=4) --->complex(kin=4)
		! int + complex(kin=8) --->complex(kin=8)

		! real4 + real4 --->rea4
		! real4 + real8 --->rea8
		! real4 + complex(kin=4) --->complex(kin=4)
		! real4 + complex(kin=8) --->complex(kin=8)

		! real8 + real8 --->rea8
		! real8 + complex(kin=4) --->complex(kin=4)
		! real8 + complex(kin=8) --->complex(kin=8)

		! complex(kin=4) + complex(kin=4) --->complex(kin=4)
		! complex(kin=4) + complex(kin=8) --->complex(kin=8)

	
	
	function  AddFunc(Dimen,Dimen2)
		class(Dimension),allocatable::AddFunc
		class(Tensor),intent(in) :: Dimen
		class(Dimension),intent(in) :: Dimen2
		integer::classtype
		allocate(AddFunc,mold=Dimen)
		select type(Dimen2)
		class is (Tensor)
			select type(AddFunc)
			class is (Tensor)
				if((.not.Dimen%getflag()).or.(.not.Dimen2%getflag()))then
					call writemess("There is no data in the Tensor,(+)",-1)
					call error_stop()
				end if
				if(Dimen%TData%getTotalData().ne.Dimen2%TData%getTotalData()) then
					call writemess("The totalData of T1 and T2 are not the same,(+)",-1)
					call writemess(Dimen%TData%getTotalData()+','+Dimen2%TData%getTotalData(),-1)
					call writemess("The program will stop",-1)
					call error_stop()
				end if
				classtype=select_type_in_add_minu(Dimen%TData,Dimen2%TData)
				call AddFunc%allocate(Dimen%Dimension,classtype)
				call add_minu_TData(AddFunc%TData,Dimen%Tdata,Dimen2%TData,1)
			class default
				call writemess(' ERROR in (+), input class error',-1)
				call error_stop
			end select
		class default
			call writemess(' ERROR in (+), input class error',-1)
			call error_stop
		end select
		return
	end function
	function  AddFunc2(Dimen,Dimenvec)
		class(Dimension),allocatable::AddFunc2
		class(Tensor),intent(in) :: Dimen
		integer,intent(in) :: Dimenvec(:)
		call writemess('ERROR in (+), DO not have the type of array + Tensor',-1)
		call error_stop
	end function
	function  AddFunc3(Dimenvec,Dimen)
		class(Dimension),allocatable::AddFunc3
		class(Tensor),intent(in) :: Dimen
		integer,intent(in) :: Dimenvec(:)
		call writemess('ERROR in (+), DO not have the type of array + Tensor',-1)
		call error_stop
	end function
	function add_int(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		integer,intent(in)::num
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,1)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_int(add%TData,T1%Tdata,num,1)
		return
	end function
	function add_int_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		integer,intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,1)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_int_TData(add%TData,num,T1%Tdata,1)
		return
	end function
	function add_real4(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		real(kind=4),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,2)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_real4(add%TData,T1%Tdata,num,1)
		return
	end function
	function add_real4_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		real(kind=4),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,2)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_real4_TData(add%TData,num,T1%Tdata,1)
		return
	end function
	function add_real8(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		real(kind=8),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,3)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_real8(add%TData,T1%Tdata,num,1)
		return
	end function
	function add_real8_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		real(kind=8),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,3)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_real8_TData(add%TData,num,T1%Tdata,1)
		return
	end function
	function add_com4(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		complex(kind=4),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,4)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_com4(add%TData,T1%Tdata,num,1)
		return
	end function
	function add_com4_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		complex(kind=4),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,4)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_com4_TData(add%TData,num,T1%Tdata,1)
		return
	end function
	function add_com8(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		complex(kind=8),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,5)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_com8(add%TData,T1%Tdata,num,1)
		return
	end function
	function add_com8_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		complex(kind=8),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,5)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_com8_TData(add%TData,num,T1%Tdata,1)
		return
	end function
	function add_char(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		character(len=*),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,7)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_char(add%TData,T1%Tdata,num,1)
		return
	end function
	function add_char_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		character(len=*),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(+)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,7)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_char_TData(add%TData,num,T1%Tdata,1)
		return
	end function
	function minus(T1,T2)
		class(Tensor),intent(in) :: T1
		class(Tensor),intent(in) :: T2
		type(Dimension)::dim1,dim2
		integer::classtype
		class(Tensor),allocatable::minus
		allocate(minus,mold=T1)
		if((.not.T1%getflag()).or.(.not.T2%getflag()))then
			call writemess("There is no data in the Tensor,(-)",-1)
			call error_stop()
		end if
		dim1=T1%Dimension
		dim2=T2%Dimension
		if(.not.(dim1.equ.dim2)) then
			call writemess("The dimension of T1 and T2 are not the same,in (-)")
			call writemess("The program will stop")
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,T2%TData)
		call minus%allocate(T1%Dimension,classtype)
		call add_minu_TData(minus%TData,T1%Tdata,T2%TData,-1)
		return
	end function
	function minus_int(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		integer,intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(-)")
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,1)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_int(add%TData,T1%Tdata,num,-1)
		return
	end function
	function minus_int_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		integer,intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(-)")
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,1)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_int_TData(add%TData,num,T1%Tdata,-1)
		return
	end function
	function minus_real4(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		real(kind=4),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(-)")
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,2)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_real4(add%TData,T1%Tdata,num,-1)
		return
	end function
	function minus_real4_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		real(kind=4),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(-)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,2)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_real4_TData(add%TData,num,T1%Tdata,-1)
		return
	end function
	function minus_real8(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		real(kind=8),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(-)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,3)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_real8(add%TData,T1%Tdata,num,-1)
		return
	end function
	function minus_real8_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		real(kind=8),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(-)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,3)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_real8_TData(add%TData,num,T1%Tdata,-1)
		return
	end function
	function minus_com4(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		complex(kind=4),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(-)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,4)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_com4(add%TData,T1%Tdata,num,-1)
		return
	end function
	function minus_com4_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		complex(kind=4),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(-)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,4)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_com4_TData(add%TData,num,T1%Tdata,-1)
		return
	end function
	function minus_com8(T1,num)result(add)
		class(Tensor),intent(in) :: T1
		complex(kind=8),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(-)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,5)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_TData_com8(add%TData,T1%Tdata,num,-1)
		return
	end function
	function minus_com8_(num,T1)result(add)
		class(Tensor),intent(in) :: T1
		complex(kind=8),intent(in)::num
		type(Dimension)::dim1
		integer::classtype
		class(Tensor),allocatable::add
		allocate(add,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(-)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,5)
		call add%allocate(T1%Dimension,classtype)
		call add_minu_com8_TData(add%TData,num,T1%Tdata,-1)
		return
	end function
	
	
	function multiply_number_int(T1,num)result(Res)
		class(Tensor),intent(in) :: T1
		integer,intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(*)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,1)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_int(Res%TData,T1%TData,num)
		return
	end function
	function multiply_number_real4(T1,num)result(Res)
		class(Tensor),intent(in) :: T1
		real(kind=4),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(*)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,2)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_real4(Res%TData,T1%TData,num)
		return
	end function
	function multiply_number_real8(T1,num)result(Res)
		class(Tensor),intent(in) :: T1
		real(kind=8),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(*)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,3)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_real8(Res%TData,T1%TData,num)
		return
	end function
	function multiply_number_com4(T1,num)result(Res)
		class(Tensor),intent(in) :: T1
		complex(kind=4),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(*)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,4)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_com4(Res%TData,T1%TData,num)
		return
	end function
	function multiply_number_com8(T1,num)result(Res)
		class(Tensor),intent(in) :: T1
		complex(kind=8),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(*)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,5)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_com8(Res%TData,T1%TData,num)
		return
	end function
	
	function multiply_number_int_(num,T1)result(Res)
		class(Tensor),intent(in) :: T1
		integer,intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(*)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,1)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_int(Res%TData,T1%TData,num)
		return
	end function
	function multiply_number_real4_(num,T1)result(Res)
		class(Tensor),intent(in) :: T1
		real(kind=4),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(*)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,2)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_real4(Res%TData,T1%TData,num)
		return
	end function
	function multiply_number_real8_(num,T1)result(Res)
		class(Tensor),intent(in) :: T1
		real(kind=8),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(*)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,3)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_real8(Res%TData,T1%TData,num)
		return
	end function
	function multiply_number_com4_(num,T1)result(Res)
		class(Tensor),intent(in) :: T1
		complex(kind=4),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(*)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,4)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_com4(Res%TData,T1%TData,num)
		return
	end function
	function multiply_number_com8_(num,T1)result(Res)
		class(Tensor),intent(in) :: T1
		complex(kind=8),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(*)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,5)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_com8(Res%TData,T1%TData,num)
		return
	end function



	function divide_Tensor(T1,T) result(Res)
		class(Tensor),intent(in) :: T1
		class(Tensor),intent(in) :: T
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T)
		if(T%TData%getTotalData().ne.1)then
			call writemess("ERROR in T1/T2, T2 should be lengh=1",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,T%TData)
		call Res%allocate(T1%Dimension,max(classtype,2))
		call TDatadivideTensor(Res%TData,T1%TData,T%TData)
		return
	end function
	function int_divide_Tensor(num,T) result(Res)
		integer,intent(in) :: num
		class(Tensor),intent(in) :: T
		integer::classtype,Aclasstype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T)
		if(T%TData%getTotalData().ne.1)then
			call writemess("ERROR in number/T, T should be lengh=1",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T%TData,2)
		Aclasstype=1
		call Res%allocate(T%Dimension,classtype)
		call TDatadivideTensor2(Res%TData,num,T%TData,Aclasstype)
		return
	end function
	function real4_divide_Tensor(num,T) result(Res)
		real(kind=4),intent(in) :: num
		class(Tensor),intent(in) :: T
		integer::classtype,Aclasstype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T)
		if(T%TData%getTotalData().ne.1)then
			call writemess("ERROR in number/T, T should be lengh=1",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T%TData,2)
		Aclasstype=2
		call Res%allocate(T%Dimension,classtype)
		call TDatadivideTensor2(Res%TData,num,T%TData,Aclasstype)
		return
	end function
	function real8_divide_Tensor(num,T) result(Res)
		real(kind=8),intent(in) :: num
		class(Tensor),intent(in) :: T
		integer::classtype,Aclasstype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T)
		if(T%TData%getTotalData().ne.1)then
			call writemess("ERROR in number/T, T should be lengh=1",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T%TData,3)
		Aclasstype=3
		call Res%allocate(T%Dimension,classtype)
		call TDatadivideTensor2(Res%TData,num,T%TData,Aclasstype)
		return
	end function
	function com4_divide_Tensor(num,T) result(Res)
		complex(kind=4),intent(in) :: num
		class(Tensor),intent(in) :: T
		integer::classtype,Aclasstype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T)
		if(T%TData%getTotalData().ne.1)then
			call writemess("ERROR in number/T, T should be lengh=1",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T%TData,4)
		Aclasstype=4
		call Res%allocate(T%Dimension,classtype)
		call TDatadivideTensor2(Res%TData,num,T%TData,Aclasstype)
		return
	end function
	function com8_divide_Tensor(num,T) result(Res)
		complex(kind=8),intent(in) :: num
		class(Tensor),intent(in) :: T
		integer::classtype,Aclasstype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T)
		if(T%TData%getTotalData().ne.1)then
			call writemess("ERROR in number/T, T should be lengh=1",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T%TData,5)
		Aclasstype=5
		call Res%allocate(T%Dimension,classtype)
		call TDatadivideTensor2(Res%TData,num,T%TData,Aclasstype)
		return
	end function
	function divide_num_int(T1,num) result(Res)
		class(Tensor),intent(in) :: T1
		integer,intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(/)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,3)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_real8(Res%TData,T1%TData,1d0/dble(num))
		return
	end function
	function divide_num_real4(T1,num) result(Res)
		class(Tensor),intent(in) :: T1
		real(kind=4),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(/)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,3)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_real8(Res%TData,T1%TData,1d0/dble(num))
		return
	end function
	function divide_num_real8(T1,num) result(Res)
		class(Tensor),intent(in) :: T1
		real(kind=8),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(/)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,3)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_real8(Res%TData,T1%TData,1d0/num)
		return
	end function
	function divide_num_com4(T1,num) result(Res)
		class(Tensor),intent(in) :: T1
		complex(kind=4),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(/)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,5)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_com8(Res%TData,T1%TData,dcmplx(1d0)/dcmplx(num))
		return
	end function
	function divide_num_com8(T1,num) result(Res)
		class(Tensor),intent(in) :: T1
		complex(kind=8),intent(in) ::   num
		integer::classtype
		class(Tensor),allocatable::Res
		allocate(Res,mold=T1)
		if(.not.T1%getflag())then
			call writemess("There is no data in the Tensor,(/)",-1)
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,5)
		call Res%allocate(T1%Dimension,classtype)
		call TDatamultiply_number_com8(Res%TData,T1%TData,dcmplx(1.)/num)
		return
	end function


	!**************** ProductTensor  ***************************
		!		ProductTensor regard the last index of T1 and the first index
		!	  of T2 as the dimenion for matrix-product,other index will be see
		!	  as another dimenison.T1 and T2 can be any rank,but the last dimenion
		!	  of T1 and the first diemsion of T2 should be equal.

	function ProductTensor (T1,T2) 
		class(Tensor),allocatable::ProductTensor
		class(Tensor),intent(in) :: T1
		class(Tensor),intent(in) :: T2
		integer::rank1,rank2,flag,T1m,T1n,T2m,T2n,T1l,T2l
		type(Dimension),pointer::D1,D2,newD
		integer::i,classtype,total1,total2
		allocate(ProductTensor,mold=T1)
		if(.not.T1%getflag()) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the first Tensor"
			write(*,*)"stop"
			call error_stop()
		end if
		if(.not.T2%getflag()) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the second Tensor"
			write(*,*)"stop"
			call error_stop()
		end if	
		D1=>workingDimension1
		D2=>workingDimension2
		newD=>WorkingDimension3
		rank1=T1%Dimension%getRank()
		rank2=T2%Dimension%getRank()
		D1=T1%dimension
		D2=T2%dimension
		total1=T1%TData%getTotalData()
		total2=T2%TData%getTotalData()
		if((total1.eq.1).or.(total2.eq.1)) then
			flag=0
		else if((rank1.eq.1).and.(rank2.eq.1)) then
			flag=1
		else if((rank1.eq.1).and.(rank2.ge.2)) then
			flag=2
		else if((rank1.ge.2).and.(rank2.eq.1)) then
			flag=3	
		else if((rank1.ge.2).and.(rank2.ge.2)) then
			flag=4
		else
			write(*,*)"ERROR in ProductTensor",rank1,rank2
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,T2%TData)
		select case (flag)
			case (0)
				if((total1.eq.1).and.(total2.eq.1)) then!there is only one element in both T1 and T2
						if((rank1.eq.1).and.(rank2.eq.1))then!Number*number,(1) *(1)
							call ProductTensor%allocate(T1%Dimension,classtype)
							call product_NumNum(ProductTensor%TData,T1%TData,T2%TData)	
							return
						else if((rank1.eq.1).and.(rank2.ne.1))then!Number*Tensor,(1) * (1,1,1)
							call ProductTensor%allocate(T2%Dimension,classtype)
							call product_NumNum(ProductTensor%TData,T1%TData,T2%TData)	
							return
						else if((rank1.ne.1).and.(rank2.eq.1))then!Tensor*Number, (1,1,1) * (1)
							call ProductTensor%allocate(T1%Dimension,classtype)
							call product_NumNum(ProductTensor%TData,T1%TData,T2%TData)	
							return
						else if((rank1.ne.1).and.(rank2.ne.1))then!Tensor*Tensor, (1,1,1) * (1,1)
							if(rank1.ge.2) then
								newD=D1.subdim.[1,rank1-1]
							end if
							if(rank2.ge.2) then
								newD=newD+(D2.subdim.[2,rank2])
							end if
							call ProductTensor%allocate(newD,classtype)
							call product_NumNum(ProductTensor%TData,T1%TData,T2%TData)	
							return
						end if
				else if((total1.eq.1).and.(total2.ne.1)) then!there is only one element in both T1, but not T2
						if(rank1.eq.1)then!Number*Tensor,(1) *(3) or (1) *(2,3) 
							call ProductTensor%allocate(T2%Dimension,classtype)
							call product_Mnum_dim1(ProductTensor%TData,T2%TData,T1%TData)	
							return
						else 
							if(rank2.eq.1)then!Tensor*number,(1,1) *(3)
								call writemess("ERROR in ProductTensor,case -1,stop",-1)
								call T1%diminfo('dimension of T1')
								call T2%diminfo('dimension of T2')
								call error_stop()
							else!Tensor*Tensor,(1,1) *(1,2,1,2)
								if(T2%dim(1).ne.1)then
									call writemess("ERROR in ProductTensor,case -2,stop",-1)
									call T1%diminfo('dimension of T1')
									call T2%diminfo('dimension of T2')
									call error_stop()
								end if
								if(rank1.ge.2) then
									newD=D1.subdim.[1,rank1-1]
								end if
								if(rank2.ge.2) then
									newD=newD+(D2.subdim.[2,rank2])
								end if
								call ProductTensor%allocate(newD,classtype)
								call product_Mnum_dim1(ProductTensor%TData,T2%TData,T1%TData)	
								return
							end if
						end if
					else if((total1.ne.1).and.(total2.eq.1)) then!there is only one element in both T2, but not T1
							if(rank2.eq.1)then!Tensor*number,(3) *(1) or (2,3) *(1) 
								call ProductTensor%allocate(T1%Dimension,classtype)
								call product_Mnum_dim1(ProductTensor%TData,T1%TData,T2%TData)	
								return
							else
								if(rank1.eq.1)then!Tensor*number,(3)*(1,1)
									call writemess("ERROR in ProductTensor,case -3,stop",-1)
									call T1%diminfo('dimension of T1')
									call T2%diminfo('dimension of T2')
									call error_stop()
								else!Tensor*Tensor,(1,2,2,1) *(1,1)
									if(T1%dim(rank1).ne.1)then
										call writemess("ERROR in ProductTensor,case -4,stop",-1)
										call T1%diminfo('dimension of T1')
										call T2%diminfo('dimension of T2')
										call error_stop()
									end if
									if(rank1.ge.2) then
										newD=D1.subdim.[1,rank1-1]
									end if
									if(rank2.ge.2) then
										newD=newD+(D2.subdim.[2,rank2])
									end if
									call ProductTensor%allocate(newD,classtype)
									call product_Mnum_dim1(ProductTensor%TData,T1%TData,T2%TData)	
									return
								end if
							end if
					else
							call writemess("ERROR in ProductTensor,case -5,stop",-1)
							call T1%diminfo('dimension of T1')
							call T2%diminfo('dimension of T2')
							call error_stop()
					end if
			case (1)
				T1m=T1.dim.1
				T2n=T2.dim.1
				if(T1m.ne.T2n) then
					call writemess("ERROR in ProductTensor,case 1,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				call ProductTensor%allocate((/1/),classtype)
				call product_VV_dim1(ProductTensor%TData,T1%TData,T2%TData)	
				return
			case (2)
				if(rank2.ge.2) then
					newD=D2.subdim.[2,rank2]
					D2=D2%fuseIndex(2,rank2)
				end if
				T2m=D2%dim(1)
				T2n=D2%dim(2)
				if((D1%dim(1)) .ne. T2m) then
					call writemess("ERROR in ProductTensor,case 2,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				call ProductTensor%allocate(newD,classtype)
				if(T2n.eq.1)then![m]*[m,1]case
					call product_VV_dim1(ProductTensor%TData,T1%TData,T2%TData)	
				else
					call product_VM_dim1(ProductTensor%TData,T1%TData,T2%TData,T2m,T2n)	
				end if
				return
		case (3)
				if(rank1.ge.2) then
					newD=D1.subdim.[1,rank1-1]
					D1=D1%fuseIndex(1,rank1-2)
				end if
				T1m=D1%dim(1)
				T1n=D1%dim(2)
				if((D2%dim(1)) .ne. T1n) then
					call writemess("ERROR in ProductTensor,case 3,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				call ProductTensor%allocate(newD,classtype)
				if(T1m.eq.1)then![1,m]*[m]case
					call product_VV_dim1(ProductTensor%TData,T1%TData,T2%TData)	
				else
					call product_MV_dim1(ProductTensor%TData,T1%TData,T2%TData,T1m,T1n)	
				end if
				return
		case (4)
				if(rank1.ge.2) then
					newD=D1.subdim.[1,rank1-1]
					D1=D1%fuseIndex(1,rank1-2)
				end if
				if(rank2.ge.2) then
					newD=newD+(D2.subdim.[2,rank2])
					D2=D2%fuseIndex(2,rank2)
				end if
				if((D1%dim(2)).ne.(D2%dim(1))) then
					call writemess("ERROR in ProductTensor,case 4,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				T1m=D1%dim(1)
				T1n=D1%dim(2)
				T2m=D2%dim(1)
				T2n=D2%dim(2)
				call ProductTensor%allocate(newD,classtype)
				if((T1m.eq.1).and.(T2n.eq.1))then![1,m] [m,1]
					call product_VV_dim1(ProductTensor%TData,T1%TData,T2%TData)	
				else if((T1m.eq.1).and.(T2n.ne.1))then![1,m] [m,n]
					call product_VM_dim1(ProductTensor%TData,T1%TData,T2%TData,T2m,T2n)	
				else if((T1m.ne.1).and.(T2n.eq.1))then![m,n] [n,1]
					call product_MV_dim1(ProductTensor%TData,T1%TData,T2%TData,T1m,T1n)	
				else
					call product_MM_dim1(ProductTensor%TData,T1%TData,T2%Tdata,T1m,T1n,T2n)	
				end if
				return
		case default 
			write(*,*) "ERROR in ProductTensor,no such data"
			call error_stop()
		end 	select
		return
	end function
	
	
	subroutine ProductTensorRoutine1(Res,T1,T2,alpha,beta)! Res = alpha* T1 *  T2  + beta*Res
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: T1,T2
		class(*),intent(in)::alpha,beta
		integer::rank1,rank2,flag,T1m,T1n,T2m,T2n,T1l,T2l
		type(Dimension),pointer::D1,D2,newD
		integer::i,classtype,total1,total2
		class(Tensor),pointer::Resp,T1p,T2p
		if(.not.T1%getflag()) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the first Tensor"
			write(*,*)"stop"
			call error_stop()
		end if
		if(.not.T2%getflag()) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the second Tensor"
			write(*,*)"stop"
			call error_stop()
		end if	
		
		Resp=>Res
		T1p=>T1
		T2p=>T2
		if(associated(Resp,T1p).or.associated(Resp,T2p))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%ProductTensorRoutine(Res,T1,T2,alpha,beta)')
			call writemess('Res and T1, or Res and T2, can not be a same variable')
			call error_stop
		end if
		Resp=>null()
		T1p=>null()
		T2p=>null()
		
		
		D1=>workingDimension1
		D2=>workingDimension2
		newD=>WorkingDimension3
		rank1=T1%Dimension%getRank()
		rank2=T2%Dimension%getRank()
		D1=T1%dimension
		D2=T2%dimension
		total1=T1%TData%getTotalData()
		total2=T2%TData%getTotalData()
		if((total1.eq.1).or.(total2.eq.1)) then
			flag=0
		else if((rank1.eq.1).and.(rank2.eq.1)) then
			flag=1
		else if((rank1.eq.1).and.(rank2.ge.2)) then
			flag=2
		else if((rank1.ge.2).and.(rank2.eq.1)) then
			flag=3	
		else if((rank1.ge.2).and.(rank2.ge.2)) then
			flag=4
		else
			write(*,*)"ERROR in ProductTensor",rank1,rank2
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,T2%TData)
		select case (flag)
			case (0)
				if((total1.eq.1).and.(total2.eq.1)) then!there is only one element in both T1 and T2
						if((rank1.eq.1).and.(rank2.eq.1))then!Number*number,(1) *(1)
							if(.not.Res%getFlag())then
								call Res%allocate(T1%Dimension,classtype)
								call Res%zero()
							end if
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,beta)	
							return
						else if((rank1.eq.1).and.(rank2.ne.1))then!Number*Tensor,(1) * (1,1,1)
							if(.not.Res%getFlag())then
								call Res%allocate(T2%Dimension,classtype)
								call Res%zero()
							end if
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,beta)		
							return
						else if((rank1.ne.1).and.(rank2.eq.1))then!Tensor*Number, (1,1,1) * (1)
							if(.not.Res%getFlag())then
								call Res%allocate(T1%Dimension,classtype)
								call Res%zero()
							end if
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,beta)			
							return
						else if((rank1.ne.1).and.(rank2.ne.1))then!Tensor*Tensor, (1,1,1) * (1,1)
							if(rank1.ge.2) then
								newD=D1.subdim.[1,rank1-1]
							end if
							if(rank2.ge.2) then
								newD=newD+(D2.subdim.[2,rank2])
							end if
							if(.not.Res%getFlag())then
								call Res%allocate(newD,classtype)
								call Res%zero()
							end if
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,beta)		
							return
						end if
				else if((total1.eq.1).and.(total2.ne.1)) then!there is only one element in both T1, but not T2
						if(rank1.eq.1)then!Number*Tensor,(1) *(3) or (1) *(2,3) 
							if(.not.Res%getFlag())then
								call Res%allocate(T2%Dimension,classtype)
								call Res%zero()
							end if
							call product_Mnum_dim1_par(Res%TData,T2%TData,T1%TData,alpha,beta)			
							return
						else 
							if(rank2.eq.1)then!Tensor*number,(1,1) *(3)
								call writemess("ERROR in ProductTensor,case -1,stop",-1)
								call T1%diminfo('dimension of T1')
								call T2%diminfo('dimension of T2')
								call error_stop()
							else!Tensor*Tensor,(1,1) *(1,2,1,2)
								if(T2%dim(1).ne.1)then
									call writemess("ERROR in ProductTensor,case -2,stop",-1)
									call T1%diminfo('dimension of T1')
									call T2%diminfo('dimension of T2')
									call error_stop()
								end if
								if(rank1.ge.2) then
									newD=D1.subdim.[1,rank1-1]
								end if
								if(rank2.ge.2) then
									newD=newD+(D2.subdim.[2,rank2])
								end if
								if(.not.Res%getFlag())then
									call Res%allocate(newD,classtype)
									call Res%zero()
								end if
								call product_Mnum_dim1_par(Res%TData,T2%TData,T1%TData,alpha,beta)	
								return
							end if
						end if
					else if((total1.ne.1).and.(total2.eq.1)) then!there is only one element in both T2, but not T1
							if(rank2.eq.1)then!Tensor*number,(3) *(1) or (2,3) *(1) 
								if(.not.Res%getFlag())then
									call Res%allocate(T1%Dimension,classtype)
									call Res%zero()
								end if
								call product_Mnum_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)	
								return
							else
								if(rank1.eq.1)then!Tensor*number,(3)*(1,1)
									call writemess("ERROR in ProductTensor,case -3,stop",-1)
									call T1%diminfo('dimension of T1')
									call T2%diminfo('dimension of T2')
									call error_stop()
								else!Tensor*Tensor,(1,2,2,1) *(1,1)
									if(T1%dim(rank1).ne.1)then
										call writemess("ERROR in ProductTensor,case -4,stop",-1)
										call T1%diminfo('dimension of T1')
										call T2%diminfo('dimension of T2')
										call error_stop()
									end if
									if(rank1.ge.2) then
										newD=D1.subdim.[1,rank1-1]
									end if
									if(rank2.ge.2) then
										newD=newD+(D2.subdim.[2,rank2])
									end if
									if(.not.Res%getFlag())then
										call Res%allocate(newD,classtype)
										call Res%zero()
									end if
									call product_Mnum_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)	
									return
								end if
							end if
					else
							call writemess("ERROR in ProductTensor,case -5,stop",-1)
							call T1%diminfo('dimension of T1')
							call T2%diminfo('dimension of T2')
							call error_stop()
					end if
			case (1)
				T1m=T1.dim.1
				T2n=T2.dim.1
				if(T1m.ne.T2n) then
					call writemess("ERROR in ProductTensor,case 1,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				if(.not.Res%getFlag())then
					call Res%allocate((/1/),classtype)
					call Res%zero()
				end if
				call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)	
				return
			case (2)
				if(rank2.ge.2) then
					newD=D2.subdim.[2,rank2]
					D2=D2%fuseIndex(2,rank2)
				end if
				T2m=D2%dim(1)
				T2n=D2%dim(2)
				if((D1%dim(1)) .ne. T2m) then
					call writemess("ERROR in ProductTensor,case 2,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				if(.not.Res%getFlag()) then
					call Res%allocate(newD,classtype)
					call Res%zero()
				end if
				if(T2n.eq.1)then![m]*[m,1]case
					call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)
				else
					call product_VM_dim1_par(Res%TData,T1%TData,T2%TData,T2m,T2n,alpha,beta)	
				end if
				return
		case (3)
				if(rank1.ge.2) then
					newD=D1.subdim.[1,rank1-1]
					D1=D1%fuseIndex(1,rank1-2)
				end if
				T1m=D1%dim(1)
				T1n=D1%dim(2)
				if((D2%dim(1)) .ne. T1n) then
					call writemess("ERROR in ProductTensor,case 3,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				if(.not.Res%getFlag())then
					call Res%allocate(newD,classtype)
					call Res%zero()
				end if
				if(T1m.eq.1)then![1,m]*[m]case
					call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)
				else
					call product_MV_dim1_par(Res%TData,T1%TData,T2%TData,T1m,T1n,alpha,beta)
				end if
				return
		case (4)
				if(rank1.ge.2) then
					newD=D1.subdim.[1,rank1-1]
					D1=D1%fuseIndex(1,rank1-2)
				end if
				if(rank2.ge.2) then
					newD=newD+(D2.subdim.[2,rank2])
					D2=D2%fuseIndex(2,rank2)
				end if
				if((D1%dim(2)).ne.(D2%dim(1))) then
					call writemess("ERROR in ProductTensor,case 4,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				T1m=D1%dim(1)
				T1n=D1%dim(2)
				T2m=D2%dim(1)
				T2n=D2%dim(2)
				if(.not.Res%getFlag())then
					call Res%allocate(newD,classtype)
					call Res%zero()
				end if
				if((T1m.eq.1).and.(T2n.eq.1))then![1,m] [m,1]
					call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)
				else if((T1m.eq.1).and.(T2n.ne.1))then![1,m] [m,n]
					call product_VM_dim1_par(Res%TData,T1%TData,T2%TData,T2m,T2n,alpha,beta)	
				else if((T1m.ne.1).and.(T2n.eq.1))then![m,n] [n,1]
					call product_MV_dim1_par(Res%TData,T1%TData,T2%TData,T1m,T1n,alpha,beta)	
				else
					call product_MM_dim1_par(Res%TData,T1%TData,T2%Tdata,T1m,T1n,T2n,alpha,beta)
				end if
				return
		case default 
			write(*,*) "ERROR in ProductTensor,no such data"
			call error_stop()
		end 	select
	end subroutine 
	
	
	subroutine ProductTensorRoutine2(Res,T1,T2,alpha)! Res = alpha* T1 *  T2
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: T1,T2
		class(*),intent(in)::alpha
		integer::rank1,rank2,flag,T1m,T1n,T2m,T2n,T1l,T2l
		type(Dimension),pointer::D1,D2,newD
		integer::i,classtype,total1,total2
		class(Tensor),pointer::Resp,T1p,T2p
		if(.not.T1%getflag()) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the first Tensor"
			write(*,*)"stop"
			call error_stop()
		end if
		if(.not.T2%getflag()) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the second Tensor"
			write(*,*)"stop"
			call error_stop()
		end if	
		Resp=>Res
		T1p=>T1
		T2p=>T2
		if(associated(Resp,T1p).or.associated(Resp,T2p))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%ProductTensorRoutine(Res,T1,T2,alpha,beta)')
			call writemess('Res and T1, or Res and T2, can not be a same variable')
			call error_stop
		end if
		Resp=>null()
		T1p=>null()
		T2p=>null()
		D1=>workingDimension1
		D2=>workingDimension2
		newD=>WorkingDimension3
		call Res%empty()
		rank1=T1%Dimension%getRank()
		rank2=T2%Dimension%getRank()
		D1=T1%dimension
		D2=T2%dimension
		total1=T1%TData%getTotalData()
		total2=T2%TData%getTotalData()
		if((total1.eq.1).or.(total2.eq.1)) then
			flag=0
		else if((rank1.eq.1).and.(rank2.eq.1)) then
			flag=1
		else if((rank1.eq.1).and.(rank2.ge.2)) then
			flag=2
		else if((rank1.ge.2).and.(rank2.eq.1)) then
			flag=3	
		else if((rank1.ge.2).and.(rank2.ge.2)) then
			flag=4
		else
			write(*,*)"ERROR in ProductTensor",rank1,rank2
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,T2%TData)
		select case (flag)
			case (0)
				if((total1.eq.1).and.(total2.eq.1)) then!there is only one element in both T1 and T2
						if((rank1.eq.1).and.(rank2.eq.1))then!Number*number,(1) *(1)
							call Res%allocate(T1%Dimension,classtype)
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,0)	
							return
						else if((rank1.eq.1).and.(rank2.ne.1))then!Number*Tensor,(1) * (1,1,1)
							call Res%allocate(T2%Dimension,classtype)
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,0)		
							return
						else if((rank1.ne.1).and.(rank2.eq.1))then!Tensor*Number, (1,1,1) * (1)
								call Res%allocate(T1%Dimension,classtype)
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,0)			
							return
						else if((rank1.ne.1).and.(rank2.ne.1))then!Tensor*Tensor, (1,1,1) * (1,1)
							if(rank1.ge.2) then
								newD=D1.subdim.[1,rank1-1]
							end if
							if(rank2.ge.2) then
								newD=newD+(D2.subdim.[2,rank2])
							end if
							call Res%allocate(newD,classtype)
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,0)		
							return
						end if
				else if((total1.eq.1).and.(total2.ne.1)) then!there is only one element in both T1, but not T2
						if(rank1.eq.1)then!Number*Tensor,(1) *(3) or (1) *(2,3) 
							call Res%allocate(T2%Dimension,classtype)
							call product_Mnum_dim1_par(Res%TData,T2%TData,T1%TData,alpha,0)			
							return
						else 
							if(rank2.eq.1)then!Tensor*number,(1,1) *(3)
								call writemess("ERROR in ProductTensor,case -1,stop",-1)
								call T1%diminfo('dimension of T1')
								call T2%diminfo('dimension of T2')
								call error_stop()
							else!Tensor*Tensor,(1,1) *(1,2,1,2)
								if(T2%dim(1).ne.1)then
									call writemess("ERROR in ProductTensor,case -2,stop",-1)
									call T1%diminfo('dimension of T1')
									call T2%diminfo('dimension of T2')
									call error_stop()
								end if
								if(rank1.ge.2) then
									newD=D1.subdim.[1,rank1-1]
								end if
								if(rank2.ge.2) then
									newD=newD+(D2.subdim.[2,rank2])
								end if
								call Res%allocate(newD,classtype)
								call product_Mnum_dim1_par(Res%TData,T2%TData,T1%TData,alpha,0)	
								return
							end if
						end if
					else if((total1.ne.1).and.(total2.eq.1)) then!there is only one element in both T2, but not T1
							if(rank2.eq.1)then!Tensor*number,(3) *(1) or (2,3) *(1) 
								call Res%allocate(T1%Dimension,classtype)
								call product_Mnum_dim1_par(Res%TData,T1%TData,T2%TData,alpha,0)	
								return
							else
								if(rank1.eq.1)then!Tensor*number,(3)*(1,1)
									call writemess("ERROR in ProductTensor,case -3,stop",-1)
									call T1%diminfo('dimension of T1')
									call T2%diminfo('dimension of T2')
									call error_stop()
								else!Tensor*Tensor,(1,2,2,1) *(1,1)
									if(T1%dim(rank1).ne.1)then
										call writemess("ERROR in ProductTensor,case -4,stop",-1)
										call T1%diminfo('dimension of T1')
										call T2%diminfo('dimension of T2')
										call error_stop()
									end if
									if(rank1.ge.2) then
										newD=D1.subdim.[1,rank1-1]
									end if
									if(rank2.ge.2) then
										newD=newD+(D2.subdim.[2,rank2])
									end if
									call Res%allocate(newD,classtype)
									call product_Mnum_dim1_par(Res%TData,T1%TData,T2%TData,alpha,0)	
									return
								end if
							end if
					else
							call writemess("ERROR in ProductTensor,case -5,stop",-1)
							call T1%diminfo('dimension of T1')
							call T2%diminfo('dimension of T2')
							call error_stop()
					end if
			case (1)
				T1m=T1.dim.1
				T2n=T2.dim.1
				if(T1m.ne.T2n) then
					call writemess("ERROR in ProductTensor,case 1,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				call Res%allocate((/1/),classtype)
				call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,0)	
				return
			case (2)
				if(rank2.ge.2) then
					newD=D2.subdim.[2,rank2]
					D2=D2%fuseIndex(2,rank2)
				end if
				T2m=D2%dim(1)
				T2n=D2%dim(2)
				if((D1%dim(1)) .ne. T2m) then
					call writemess("ERROR in ProductTensor,case 2,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				call Res%allocate(newD,classtype)
				if(T2n.eq.1)then![m]*[m,1]case
					call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,0)
				else
					call product_VM_dim1_par(Res%TData,T1%TData,T2%TData,T2m,T2n,alpha,0)	
				end if
				return
		case (3)
				if(rank1.ge.2) then
					newD=D1.subdim.[1,rank1-1]
					D1=D1%fuseIndex(1,rank1-2)
				end if
				T1m=D1%dim(1)
				T1n=D1%dim(2)
				if((D2%dim(1)) .ne. T1n) then
					call writemess("ERROR in ProductTensor,case 3,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				call Res%allocate(newD,classtype)
				if(T1m.eq.1)then![1,m]*[m]case
					call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,0)
				else
					call product_MV_dim1_par(Res%TData,T1%TData,T2%TData,T1m,T1n,alpha,0)
				end if
				return
		case (4)
				if(rank1.ge.2) then
					newD=D1.subdim.[1,rank1-1]
					D1=D1%fuseIndex(1,rank1-2)
				end if
				if(rank2.ge.2) then
					newD=newD+(D2.subdim.[2,rank2])
					D2=D2%fuseIndex(2,rank2)
				end if
				if((D1%dim(2)).ne.(D2%dim(1))) then
					call writemess("ERROR in ProductTensor,case 4,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				T1m=D1%dim(1)
				T1n=D1%dim(2)
				T2m=D2%dim(1)
				T2n=D2%dim(2)
				call Res%allocate(newD,classtype)
				if((T1m.eq.1).and.(T2n.eq.1))then![1,m] [m,1]
					call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,0)
				else if((T1m.eq.1).and.(T2n.ne.1))then![1,m] [m,n]
					call product_VM_dim1_par(Res%TData,T1%TData,T2%TData,T2m,T2n,alpha,0)	
				else if((T1m.ne.1).and.(T2n.eq.1))then![m,n] [n,1]
					call product_MV_dim1_par(Res%TData,T1%TData,T2%TData,T1m,T1n,alpha,0)	
				else
					call product_MM_dim1_par(Res%TData,T1%TData,T2%Tdata,T1m,T1n,T2n,alpha,0)
				end if
				return
		case default 
			write(*,*) "ERROR in ProductTensor,no such data"
			call error_stop()
		end 	select
	end subroutine 
	
	
	
	!ProductTensorRoutine3 Not OK yet

	subroutine ProductTensorRoutine3(Res,T1,T2,alpha,beta,TRANSA,TRANSB)! Res = alpha*TRANSA( T1) * TRANSB( T2)  + beta*Res
		class(Tensor),target,intent(inout)::Res
		class(Tensor),target,intent(in) :: T1,T2
		character*1,intent(in)::TRANSA,TRANSB
		class(*),intent(in)::alpha,beta
		integer::rank1,rank2,flag,T1m,T1n,T2m,T2n,T1l,T2l
		type(Dimension),pointer::D1,D2,newD
		integer::i,classtype,total1,total2
		class(Tensor),pointer::Resp,T1p,T2p
		if(.not.T1%getflag()) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the first Tensor"
			write(*,*)"stop"
			call error_stop()
		end if
		if(.not.T2%getflag()) then
			write(*,*)"ERROR in (*)"
			write(*,*)"Tensor is no data is the second Tensor"
			write(*,*)"stop"
			call error_stop()
		end if	
		
		Resp=>Res
		T1p=>T1
		T2p=>T2
		if(associated(Resp,T1p).or.associated(Resp,T2p))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%ProductTensorRoutine(Res,T1,T2,alpha,beta)')
			call writemess('Res and T1, or Res and T2, can not be a same variable')
			call error_stop
		end if
		Resp=>null()
		T1p=>null()
		T2p=>null()
		
		
		D1=>workingDimension1
		D2=>workingDimension2
		newD=>WorkingDimension3
		rank1=T1%Dimension%getRank()
		rank2=T2%Dimension%getRank()
		D1=T1%dimension
		D2=T2%dimension
		total1=T1%TData%getTotalData()
		total2=T2%TData%getTotalData()
		if((total1.eq.1).or.(total2.eq.1)) then
			flag=0
		else if((rank1.eq.1).and.(rank2.eq.1)) then
			flag=1
		else if((rank1.eq.1).and.(rank2.ge.2)) then
			flag=2
		else if((rank1.ge.2).and.(rank2.eq.1)) then
			flag=3	
		else if((rank1.ge.2).and.(rank2.ge.2)) then
			flag=4
		else
			write(*,*)"ERROR in ProductTensor",rank1,rank2
			call error_stop()
		end if
		classtype=select_type_in_add_minu(T1%TData,T2%TData)
		select case (flag)
			case (0)
				if((total1.eq.1).and.(total2.eq.1)) then!there is only one element in both T1 and T2
						if((rank1.eq.1).and.(rank2.eq.1))then!Number*number,(1) *(1)
							if(.not.Res%getFlag())then
								call Res%allocate(T1%Dimension,classtype)
								call Res%zero()
							end if
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,beta)	
							return
						else if((rank1.eq.1).and.(rank2.ne.1))then!Number*Tensor,(1) * (1,1,1)
							if(.not.Res%getFlag())then
								call Res%allocate(T2%Dimension,classtype)
								call Res%zero()
							end if
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,beta)		
							return
						else if((rank1.ne.1).and.(rank2.eq.1))then!Tensor*Number, (1,1,1) * (1)
							if(.not.Res%getFlag())then
								call Res%allocate(T1%Dimension,classtype)
								call Res%zero()
							end if
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,beta)			
							return
						else if((rank1.ne.1).and.(rank2.ne.1))then!Tensor*Tensor, (1,1,1) * (1,1)
							if(rank1.ge.2) then
								newD=D1.subdim.[1,rank1-1]
							end if
							if(rank2.ge.2) then
								newD=newD+(D2.subdim.[2,rank2])
							end if
							if(.not.Res%getFlag())then
								call Res%allocate(newD,classtype)
								call Res%zero()
							end if
							call product_NumNum_par(Res%TData,T1%TData,T2%TData,alpha,beta)		
							return
						end if
				else if((total1.eq.1).and.(total2.ne.1)) then!there is only one element in both T1, but not T2
						if(rank1.eq.1)then!Number*Tensor,(1) *(3) or (1) *(2,3) 
							if(.not.Res%getFlag())then
								call Res%allocate(T2%Dimension,classtype)
								call Res%zero()
							end if
							call product_Mnum_dim1_par(Res%TData,T2%TData,T1%TData,alpha,beta)			
							return
						else 
							if(rank2.eq.1)then!Tensor*number,(1,1) *(3)
								call writemess("ERROR in ProductTensor,case -1,stop",-1)
								call T1%diminfo('dimension of T1')
								call T2%diminfo('dimension of T2')
								call error_stop()
							else!Tensor*Tensor,(1,1) *(1,2,1,2)
								if(T2%dim(1).ne.1)then
									call writemess("ERROR in ProductTensor,case -2,stop",-1)
									call T1%diminfo('dimension of T1')
									call T2%diminfo('dimension of T2')
									call error_stop()
								end if
								if(rank1.ge.2) then
									newD=D1.subdim.[1,rank1-1]
								end if
								if(rank2.ge.2) then
									newD=newD+(D2.subdim.[2,rank2])
								end if
								if(.not.Res%getFlag())then
									call Res%allocate(newD,classtype)
									call Res%zero()
								end if
								call product_Mnum_dim1_par(Res%TData,T2%TData,T1%TData,alpha,beta)	
								return
							end if
						end if
					else if((total1.ne.1).and.(total2.eq.1)) then!there is only one element in both T2, but not T1
							if(rank2.eq.1)then!Tensor*number,(3) *(1) or (2,3) *(1) 
								if(.not.Res%getFlag())then
									call Res%allocate(T1%Dimension,classtype)
									call Res%zero()
								end if
								call product_Mnum_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)	
								return
							else
								if(rank1.eq.1)then!Tensor*number,(3)*(1,1)
									call writemess("ERROR in ProductTensor,case -3,stop",-1)
									call T1%diminfo('dimension of T1')
									call T2%diminfo('dimension of T2')
									call error_stop()
								else!Tensor*Tensor,(1,2,2,1) *(1,1)
									if(T1%dim(rank1).ne.1)then
										call writemess("ERROR in ProductTensor,case -4,stop",-1)
										call T1%diminfo('dimension of T1')
										call T2%diminfo('dimension of T2')
										call error_stop()
									end if
									if(rank1.ge.2) then
										newD=D1.subdim.[1,rank1-1]
									end if
									if(rank2.ge.2) then
										newD=newD+(D2.subdim.[2,rank2])
									end if
									if(.not.Res%getFlag())then
										call Res%allocate(newD,classtype)
										call Res%zero()
									end if
									call product_Mnum_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)	
									return
								end if
							end if
					else
							call writemess("ERROR in ProductTensor,case -5,stop",-1)
							call T1%diminfo('dimension of T1')
							call T2%diminfo('dimension of T2')
							call error_stop()
					end if
			case (1)
				T1m=T1.dim.1
				T2n=T2.dim.1
				if(T1m.ne.T2n) then
					call writemess("ERROR in ProductTensor,case 1,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				if(.not.Res%getFlag())then
					call Res%allocate((/1/),classtype)
					call Res%zero()
				end if
				call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)	
				return
			case (2)
				if(rank2.ge.2) then
					newD=D2.subdim.[2,rank2]
					D2=D2%fuseIndex(2,rank2)
				end if
				T2m=D2%dim(1)
				T2n=D2%dim(2)
				if((D1%dim(1)) .ne. T2m) then
					call writemess("ERROR in ProductTensor,case 2,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				if(.not.Res%getFlag()) then
					call Res%allocate(newD,classtype)
					call Res%zero()
				end if
				if(T2n.eq.1)then![m]*[m,1]case
					call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)
				else
					call product_VM_dim1_par(Res%TData,T1%TData,T2%TData,T2m,T2n,alpha,beta)	
				end if
				return
		case (3)
				if(rank1.ge.2) then
					newD=D1.subdim.[1,rank1-1]
					D1=D1%fuseIndex(1,rank1-2)
				end if
				T1m=D1%dim(1)
				T1n=D1%dim(2)
				if((D2%dim(1)) .ne. T1n) then
					call writemess("ERROR in ProductTensor,case 3,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				if(.not.Res%getFlag())then
					call Res%allocate(newD,classtype)
					call Res%zero()
				end if
				if(T1m.eq.1)then![1,m]*[m]case
					call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)
				else
					call product_MV_dim1_par(Res%TData,T1%TData,T2%TData,T1m,T1n,alpha,beta)
				end if
				return
		case (4)
				if(rank1.ge.2) then
					newD=D1.subdim.[1,rank1-1]
					D1=D1%fuseIndex(1,rank1-2)
				end if
				if(rank2.ge.2) then
					newD=newD+(D2.subdim.[2,rank2])
					D2=D2%fuseIndex(2,rank2)
				end if
				if((D1%dim(2)).ne.(D2%dim(1))) then
					call writemess("ERROR in ProductTensor,case 4,stop",-1)
					call T1%diminfo('dimension of T1')
					call T2%diminfo('dimension of T2')
					call error_stop()
				end if
				T1m=D1%dim(1)
				T1n=D1%dim(2)
				T2m=D2%dim(1)
				T2n=D2%dim(2)
				if(.not.Res%getFlag())then
					call Res%allocate(newD,classtype)
					call Res%zero()
				end if
				if((T1m.eq.1).and.(T2n.eq.1))then![1,m] [m,1]
					call product_VV_dim1_par(Res%TData,T1%TData,T2%TData,alpha,beta)
				else if((T1m.eq.1).and.(T2n.ne.1))then![1,m] [m,n]
					call product_VM_dim1_par(Res%TData,T1%TData,T2%TData,T2m,T2n,alpha,beta)	
				else if((T1m.ne.1).and.(T2n.eq.1))then![m,n] [n,1]
					call product_MV_dim1_parameter(Res%TData,T1%TData,T2%TData,T1m,T1n,alpha,beta,TRANSA)	
				else
					call product_MM_dim1_parameter(Res%TData,T1%TData,T2%Tdata,T1m,T1n,T2n,alpha,beta,TRANSA,TRANSB)
				end if
				return
		case default 
			write(*,*) "ERROR in ProductTensor,no such data"
			call error_stop()
		end 	select
	end subroutine 


	!**************************************************************************
	!**************************************************************************
	!
	!                                  reshape tensors
	!
	!**************************************************************************
	!**************************************************************************

	subroutine reset_dim_no_check1(Ten,dimen)
		class(Tensor),intent(inout)::Ten
		type(dimension),intent(in)::dimen
		Ten%Dimension=dimen
		call Ten%reset_total_Data_no_check(dimen%size())
		return
	end subroutine	
	subroutine reset_dim_no_check2(Ten,dimen)
		class(Tensor),intent(inout)::Ten
		integer,intent(in)::dimen(:)
		Ten%Dimension=dimen
		call Ten%reset_total_Data_no_check(product(dimen))
		return
	end subroutine	

	!**************** fuse   ****************
	!combine two index of the Tensor,which is con_index and con_index+1

	function fuseDimension_val(Dimen,con_index)result(fuseDim)
		class(Dimension),allocatable::fuseDim
		class(Tensor),intent(in) :: Dimen
		integer,intent(in) :: con_index
		allocate(fuseDim,mold=Dimen)
		select type(fuseDim)
		class is (Tensor)
			if(.not.Dimen%getFlag())then
				call fuseDim%empty()
				return
			end if
			call fuseDim%allocate(Dimen%Dimension.fuse.con_index,Dimen%getType())
			call assignmentTData(fuseDim%TData,Dimen%TData)
		class default
			call writemess('ERROR in operator(.fuse.):input class error',-1)
			call error_stop
		end select
		return
	end function
	function fuseDimension_vec(dimen,vector)result(fuseDim)
		class(Dimension),allocatable::fuseDim
		integer,intent(in) ::vector(2)
		class(Tensor),intent(in) :: dimen
		allocate(fuseDim,mold=Dimen)
		select type(fuseDim)
		class is (Tensor)
			if(.not.dimen%getFlag())then
				call fuseDim%empty()
				return
			end if
			call fuseDim%allocate(dimen%Dimension.fuse.vector,dimen%getType())
			call assignmentTData(fuseDim%TData,dimen%TData)
		class default
			call writemess('ERROR in operator(.fuse.):input class error',-1)
			call error_stop
		end select
		return
	end function		

	function splitDimension2(dimen,vector)result(splitDimension)
		class(Dimension),allocatable::splitDimension
		class(Tensor),intent(in) :: dimen
		integer,intent(in) ::vector(2)
		allocate(splitDimension,mold=Dimen)
		select type(splitDimension)
		class is (Tensor)
			if(.not.dimen%getFlag())then
				call splitDimension%empty()
				return
			end if
			call splitDimension%allocate(dimen%Dimension.split.vector,dimen%getType())
			call assignmentTData(splitDimension%TData,dimen%TData)
		class default
			call writemess('ERROR in operator(.split.):input class error',-1)
			call error_stop
		end select
		return
	end function
			
	function splitDimension3(dimen,de_index)result(splitDimension)
		class(Dimension),allocatable::splitDimension
		class(Tensor),intent(in) :: dimen
		integer,intent(in) :: de_index
		allocate(splitDimension,mold=Dimen)
		select type(splitDimension)
		class is (Tensor)
			if(.not.dimen%getFlag())then
				call splitDimension%empty()
				return
			end if
			call splitDimension%allocate(dimen%Dimension.split.de_index,dimen%getType())
			call assignmentTData(splitDimension%TData,dimen%TData)
		class default
			call writemess('ERROR in operator(.split.):input class error',-1)
			call error_stop
		end select
		return
	end function
			
	function splitDimensionAll(Dimen)result(splitDimension)
		class(Dimension),allocatable::splitDimension
		class(Tensor),intent(in) :: Dimen
		allocate(splitDimension,mold=Dimen)
		select type(splitDimension)
			class is (Tensor)
			if(.not.Dimen%getFlag())then
				call splitDimension%empty()
				return
			end if
			call splitDimension%allocate(.split.Dimen%Dimension,Dimen%getType())
			call assignmentTData(splitDimension%TData,Dimen%TData)
		class default
			call writemess('ERROR in operator(.split.):input class error',-1)
			call error_stop
		end select
		return
	end function	


	!**************************************************************************
	!**************************************************************************
	!
	!                      logical function
	!
	!**************************************************************************
	!**************************************************************************


	logical function equal_of_dim(dim1,dim2)
		class(Tensor),intent(in) :: dim1
		class(Dimension),intent(in) :: dim2
		type(Tensor) :: T
		integer :: l,flag
		select type(dim2)
		class is (Tensor)
			equal_of_dim=.true.
			if(.not.dim1%getFlag())then
				call writemess('There is no data in T1,(.eq.)',-1)
				call error_stop()
			end if
			if(.not.dim2%getFlag())then
				call writemess('There is no data in T2,(.eq.)',-1)
				call error_stop()
			end if
			if(dim1%Dimension%getRank().ne.dim2%Dimension%getRank()) then
				equal_of_dim=.false.
				return
			end if
			if(.not.(dim1%Dimension.equ.dim2%Dimension)) then
				equal_of_dim=.false.
				return
			end if
			flag=dim1%getType()*10+dim2%getType()
			select case(dim1%getType())
				case (11)
					do l=1,dim1%TData%gettotalData()
						if(dim1%ii(l).ne.dim2%ii(l))then
							equal_of_dim=.false.
						return
						end if
					end do
				case default
					T=dim1-dim2
					if(T%dmax('maxa').gt.default_zero_double_number)then
						equal_of_dim=.false.
						return
					end if
					return
			end select
		class default
			call writemess('ERROR in operator(.equ.), class error',-1)
			call error_stop
		end select
	end function
	logical function T_eq_int(T,num)
		class(Tensor),intent(in) :: T
		integer,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.eq.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).eq.num
			case (2)
				T_eq_int=T%si(1).eq.num
			case (3)
				T_eq_int=T%di(1).eq.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.eq.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function int_eq_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		integer,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.eq.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).eq.num
			case (2)
				T_eq_int=T%si(1).eq.num
			case (3)
				T_eq_int=T%di(1).eq.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.eq.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_eq_real4(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=4),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.eq.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).eq.num
			case (2)
				T_eq_int=T%si(1).eq.num
			case (3)
				T_eq_int=T%di(1).eq.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.eq.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function real4_eq_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=4),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.eq.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).eq.num
			case (2)
				T_eq_int=T%si(1).eq.num
			case (3)
				T_eq_int=T%di(1).eq.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.eq.)!",-1)
				call error_stop()
		end select
		return
	end function			
	logical function T_eq_real8(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=8),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.eq.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).eq.num
			case (2)
				T_eq_int=T%si(1).eq.num
			case (3)
				T_eq_int=T%di(1).eq.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.eq.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function real8_eq_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=8),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.eq.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			T_eq_int=.false.
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).eq.num
			case (2)
				T_eq_int=T%si(1).eq.num
			case (3)
				T_eq_int=T%di(1).eq.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.eq.)!",-1)
				call error_stop()
		end select
		return
	end function		

	!use only there is one element and Tensor are real or integer

	logical function le_of_Tensor(T1,T2)
		class(Tensor),intent(in) :: T1
		class(Tensor),intent(in) :: T2
		integer :: l
		if(.not.T1%getFlag())then
			call writemess('There is no data in T1,(.le.)',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in T2,(.le.)',-1)
			call error_stop()
		end if
		if(T1%TData%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.le.)!!",-1)
			call error_stop()
		end if
		if(T1%TData%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.le.)!!",-1)
			call error_stop()
		end if
		if(T1%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.le.)!",-1)
			call error_stop()
		end if
		if(T2%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.le.)!",-1)
			call error_stop()
		end if
		le_of_Tensor=T1%di(1).le.T2%di(1)
		return
	end function
	
	logical function T_le_int(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		integer,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.le.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.le.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).le.num
			case (2)
				T_eq_int=T%si(1).le.num
			case (3)
				T_eq_int=T%di(1).le.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.le.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function int_le_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		integer,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.le.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.le.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.le.T%ii(1)
			case (2)
				T_eq_int=num.le.T%si(1)
			case (3)
				T_eq_int=num.le.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.le.)!",-1)
				call error_stop()
		end select
		return
	end function
	
	logical function T_le_real4(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=4),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.le.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.le.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).le.num
			case (2)
				T_eq_int=T%si(1).le.num
			case (3)
				T_eq_int=T%di(1).le.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.le.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function real4_le_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=4),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.le.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.le.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.le.T%ii(1)
			case (2)
				T_eq_int=num.le.T%si(1)
			case (3)
				T_eq_int=num.le.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.le.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_le_real8(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=8),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.le.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.le.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).le.num
			case (2)
				T_eq_int=T%si(1).le.num
			case (3)
				T_eq_int=T%di(1).le.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.le.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function real8_le_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=8),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.le.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.le.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.le.T%ii(1)
			case (2)
				T_eq_int=num.le.T%si(1)
			case (3)
				T_eq_int=num.le.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.le.)!",-1)
				call error_stop()
		end select
		return
	end function
	
	logical function lt_of_Tensor(T1,T2)
		class(Tensor),intent(in) :: T1
		class(Tensor),intent(in) :: T2
		integer :: l
		if(.not.T1%getFlag())then
			call writemess('There is no data in T1,(.lt.)',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in T2,(.lt.)',-1)
			call error_stop()
		end if
		if(T1%TData%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.lt.)!!",-1)
			call error_stop()
		end if
		if(T1%TData%gettotalData().gt.2)then
			call writemess("ONLY FOR Tensor with one element,(.lt.)!!",-1)
			call error_stop()
		end if
		if(T1%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.lt.)!",-1)
			call error_stop()
		end if
		if(T2%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.lt.)!",-1)
			call error_stop()
		end if
		lt_of_Tensor=T1%di(1).lt.T2%di(1)
		return
	end function
	
	logical function T_lt_int(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		integer,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.lt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.lt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).lt.num
			case (2)
				T_eq_int=T%si(1).lt.num
			case (3)
				T_eq_int=T%di(1).lt.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.lt.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function int_lt_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		integer,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.lt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.lt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.lt.T%ii(1)
			case (2)
				T_eq_int=num.lt.T%si(1)
			case (3)
				T_eq_int=num.lt.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.lt.)!",-1)
				call error_stop()
		end select
		return
	end function
	
	logical function T_lt_real4(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=4),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.lt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.lt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).lt.num
			case (2)
				T_eq_int=T%si(1).lt.num
			case (3)
				T_eq_int=T%di(1).lt.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.lt.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function real4_lt_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=4),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.lt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.lt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.lt.T%ii(1)
			case (2)
				T_eq_int=num.lt.T%si(1)
			case (3)
				T_eq_int=num.lt.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.lt.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_lt_real8(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=8),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.lt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.lt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).lt.num
			case (2)
				T_eq_int=T%si(1).lt.num
			case (3)
				T_eq_int=T%di(1).lt.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.lt.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function real8_lt_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=8),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.lt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.lt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.lt.T%ii(1)
			case (2)
				T_eq_int=num.lt.T%si(1)
			case (3)
				T_eq_int=num.lt.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.lt.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function ge_of_Tensor(T1,T2)
		class(Tensor),intent(in) :: T1
		class(Tensor),intent(in) :: T2
		integer :: l
		if(.not.T1%getFlag())then
			call writemess('There is no data in T1,(.ge.)',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in T2,(.ge.)',-1)
			call error_stop()
		end if
		if(T1%TData%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.ge.)!!",-1)
			call error_stop()
		end if
		if(T1%TData%gettotalData().gt.2)then
			call writemess("ONLY FOR Tensor with one element,(.ge.)!!",-1)
			call error_stop()
		end if
		if(T1%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.ge.)!",-1)
			call error_stop()
		end if
		if(T2%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.ge.)!",-1)
			call error_stop()
		end if
		ge_of_Tensor=T1%di(1).ge.T2%di(1)
		return
	end function
	
	logical function T_ge_int(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		integer,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.ge.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.ge.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).ge.num
			case (2)
				T_eq_int=T%si(1).ge.num
			case (3)
				T_eq_int=T%di(1).ge.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.ge.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function int_ge_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		integer,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.ge.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.ge.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.ge.T%ii(1)
			case (2)
				T_eq_int=num.ge.T%si(1)
			case (3)
				T_eq_int=num.ge.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.ge.)!",-1)
				call error_stop()
		end select
		return
	end function
	
	logical function T_ge_real4(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=4),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.ge.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.ge.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).ge.num
			case (2)
				T_eq_int=T%si(1).ge.num
			case (3)
				T_eq_int=T%di(1).ge.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.ge.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function real4_ge_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=4),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.ge.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.ge.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.ge.T%ii(1)
			case (2)
				T_eq_int=num.ge.T%si(1)
			case (3)
				T_eq_int=num.ge.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.ge.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_ge_real8(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=8),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.ge.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.ge.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).ge.num
			case (2)
				T_eq_int=T%si(1).ge.num
			case (3)
				T_eq_int=T%di(1).ge.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.ge.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function real8_ge_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=8),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.ge.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.ge.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.ge.T%ii(1)
			case (2)
				T_eq_int=num.ge.T%si(1)
			case (3)
				T_eq_int=num.ge.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.ge.)!",-1)
				call error_stop()
		end select
		return
	end function
	
	logical function gt_of_Tensor(T1,T2)
		class(Tensor),intent(in) :: T1,T2
		integer :: l
		if(.not.T1%getFlag())then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		if(T1%TData%gettotalData().gt.1)then
			call writemess("ONLY FOR Tensor with one element,(.gt.)!!",-1)
			call error_stop()
		end if
		if(T1%TData%gettotalData().gt.2)then
			call writemess("ONLY FOR Tensor with one element,(.gt.)!!",-1)
			call error_stop()
		end if
		if(T1%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
			call error_stop()
		end if
		if(T2%getType().gt.3)then
			call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
			call error_stop()
		end if
		gt_of_Tensor=T1%di(1).gt.T2%di(1)
		return
	end function
	
	logical function T_gt_int(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		integer,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).gt.num
			case (2)
				T_eq_int=T%si(1).gt.num
			case (3)
				T_eq_int=T%di(1).gt.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function int_gt_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		integer,intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.gt.T%ii(1)
			case (2)
				T_eq_int=num.gt.T%si(1)
			case (3)
				T_eq_int=num.gt.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function
	
	logical function T_gt_real4(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=4),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).gt.num
			case (2)
				T_eq_int=T%si(1).gt.num
			case (3)
				T_eq_int=T%di(1).gt.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function real4_gt_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=4),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.gt.T%ii(1)
			case (2)
				T_eq_int=num.gt.T%si(1)
			case (3)
				T_eq_int=num.gt.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function
	logical function T_gt_real8(T,num)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=8),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T1,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=T%ii(1).gt.num
			case (2)
				T_eq_int=T%si(1).gt.num
			case (3)
				T_eq_int=T%di(1).gt.num
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function	
	logical function real8_gt_T(num,T)result(T_eq_int)
		class(Tensor),intent(in) :: T
		real(kind=8),intent(in)::num
		if(.not.T%getFlag())then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		if(T%TData%gettotalData().gt.1)then
			call writemess('There is no data in T2,(.gt.)',-1)
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				T_eq_int=num.gt.T%ii(1)
			case (2)
				T_eq_int=num.gt.T%si(1)
			case (3)
				T_eq_int=num.gt.T%di(1)
			case default
				call writemess("ONLY FOR Tensor of integer or real,(.gt.)!",-1)
				call error_stop()
		end select
		return
	end function


	!****************************************************************
	!********************************************************************
	!
	!                 max or min element in Tensor
	!
	!*******************************************************************
	!*******************************************************************

	integer function intmaxElement(T)
		class(Tensor),intent(in) :: T
		call intmaxTData(intmaxElement,T%TData)
		return
	end function
	function realmaxElement(T)
		real*4 ::realmaxElement
		class(Tensor),intent(in) :: T
		call real4_maxTData(realmaxElement,T%TData)
		return
	end function
	function dblemaxElement(T)
		real*8 ::dblemaxElement
		class(Tensor),intent(in) :: T
		call real8_maxTData(dblemaxElement,T%TData)
		return
	end function
	function cmplxmaxElement(T)
		complex(kind=4) ::cmplxmaxElement
		class(Tensor),intent(in) :: T
		call com4_maxTData(cmplxmaxElement,T%TData)
		return
	end function
	function dcmplxmaxElement(T)
		complex(kind=8) ::dcmplxmaxElement
		class(Tensor),intent(in) :: T
		call com8_maxTData(dcmplxmaxElement,T%TData)
		return
	end function
	function TmaxElement(T)
		class(Tensor),allocatable :: TmaxElement
		class(Tensor),intent(in) :: T
		allocate(TmaxElement,mold=T)
		select case(T%getType())
			case(1)
				TmaxElement=intmaxElement(T)
			case(2)
				TmaxElement=realmaxElement(T)
			case(3)
				TmaxElement=dblemaxElement(T)
			case(4)
				TmaxElement=cmplxmaxElement(T)
			case(5)
				TmaxElement=dcmplxmaxElement(T)
			case default
				call writemess("ERROR in type of input Tensor,(max)",-1)
				call error_stop()
		end select
		return
	end function
	integer function intminElement(T)
		class(Tensor),intent(in) :: T
		call intminTData(intminElement,T%TData)
		return
	end function
	function realminElement(T)
		real*4 ::realminElement
		class(Tensor),intent(in) :: T
		call real4_minTData(realminElement,T%TData)
		return
	end function
	function dbleminElement(T)
		real*8 ::dbleminElement
		class(Tensor),intent(in) :: T
		call real8_minTData(dbleminElement,T%TData)
		return
	end function
	function cmplxminElement(T)
		complex(kind=4) ::cmplxminElement
		class(Tensor),intent(in) :: T
		call com4_minTData(cmplxminElement,T%TData)
		return
	end function
	function dcmplxminElement(T)
		complex(kind=8) ::dcmplxminElement
		class(Tensor),intent(in) :: T
		call com8_minTData(dcmplxminElement,T%TData)
		return
	end function
	function TminElement(T)
		class(Tensor),allocatable ::TminElement
		class(Tensor),intent(in) :: T
		allocate(TminElement,mold=T)
		select case(T%getType())
			case(1)
				TminElement=intminElement(T)
			case(2)
				TminElement=realminElement(T)
			case(3)
				TminElement=dbleminElement(T)
			case(4)
				TminElement=cmplxminElement(T)
			case(5)
				TminElement=dcmplxminElement(T)
			case default
				call writemess("ERROR in type of input Tensor,(min)",-1)
				call error_stop()
		end select
		return
	end function

	!maxminflag=
		!'maxa': max abs 
		!'mina': min abs 
		!'maxr': max real
		!'minr': min real
		!'maxi': 0(not com) or max imag
		!'mini': 0(not com) or min imag

	integer function intmaxminElement(T,maxminflag)result(Res)
		class(Tensor),intent(in) :: T
		character(len=4),intent(in)::maxminflag
		select case(maxminflag)
			case('maxa')
				call intmaxabsTData(Res,T%TData)
			case('mina')
				call intminabsTData(Res,T%TData)
			case('maxr')
				call intmaxrealTData(Res,T%TData)
			case('minr')
				call intminrealTData(Res,T%TData)
			case('maxi')
				call intmaximagTData(Res,T%TData)
			case('mini')
				call intminimagTData(Res,T%TData)
		end select
		return
	end function
	function realmaxminElement(T,maxminflag)result(Res)
		real(kind=4) ::Res
		class(Tensor),intent(in) :: T
		character(len=4),intent(in)::maxminflag
		select case(maxminflag)
			case('maxa')
				call real4_maxabsTData(Res,T%TData)
			case('mina')
				call real4_minabsTData(Res,T%TData)
			case('maxr')
				call real4_maxrealTData(Res,T%TData)
			case('minr')
				call real4_minrealTData(Res,T%TData)
			case('maxi')
				call real4_maximagTData(Res,T%TData)
			case('mini')
				call real4_minimagTData(Res,T%TData)
		end select
		return
	end function
	function dblemaxminElement(T,maxminflag)result(Res)
		real(kind=8) ::Res
		class(Tensor),intent(in) :: T
		character(len=4),intent(in)::maxminflag
		select case(maxminflag)
			case('maxa')
				call real8_maxabsTData(Res,T%TData)
			case('mina')
				call real8_minabsTData(Res,T%TData)
			case('maxr')
				call real8_maxrealTData(Res,T%TData)
			case('minr')
				call real8_minrealTData(Res,T%TData)
			case('maxi')
				call real8_maximagTData(Res,T%TData)
			case('mini')
				call real8_minimagTData(Res,T%TData)
		end select
		return
	end function
	function cmplxmaxminElement(T,maxminflag)result(Res)
		complex(kind=4) ::Res
		class(Tensor),intent(in) :: T
		character(len=4),intent(in)::maxminflag
		select case(maxminflag)
			case('maxa')
				call com4_maxabsTData(Res,T%TData)
			case('mina')
				call com4_minabsTData(Res,T%TData)
			case('maxr')
				call com4_maxrealTData(Res,T%TData)
			case('minr')
				call com4_minrealTData(Res,T%TData)
			case('maxi')
				call com4_maximagTData(Res,T%TData)
			case('mini')
				call com4_minimagTData(Res,T%TData)
		end select
		return
	end function
	function dcmplxmaxminElement(T,maxminflag)result(Res)
		complex(kind=8) ::Res
		class(Tensor),intent(in) :: T
		character(len=4),intent(in)::maxminflag
		select case(maxminflag)
			case('maxa')
				call com8_maxabsTData(Res,T%TData)
			case('mina')
				call com8_minabsTData(Res,T%TData)
			case('maxr')
				call com8_maxrealTData(Res,T%TData)
			case('minr')
				call com8_minrealTData(Res,T%TData)
			case('maxi')
				call com8_maximagTData(Res,T%TData)
			case('mini')
				call com8_minimagTData(Res,T%TData)
		end select
		return
	end function
	function TmaxminElement(T,maxminflag)
		class(Tensor),allocatable ::TmaxminElement
		class(Tensor),intent(in) :: T
		character(len=4),intent(in)::maxminflag
		real*4::s
		real*8::D
		allocate(TmaxminElement,mold=T)
		select case(T%getType())
			case(1)
				TmaxminElement=intmaxminElement(T,maxminflag)
			case(2)
				TmaxminElement=realmaxminElement(T,maxminflag)
			case(3)
				TmaxminElement=dblemaxminElement(T,maxminflag)
			case(4)
				s=cmplxmaxminElement(T,maxminflag)
				TmaxminElement=s
			case(5)
				d=dcmplxmaxminElement(T,maxminflag)
				TmaxminElement=d
			case default
				call writemess("ERROR in type of input Tensor,(maxmin)",-1)
				call error_stop()
		end select
		return
	end function
	
	
	integer function sumTensori(T)result(Res)
		class(Tensor),intent(in) :: T
		call sum_TDatai(Res,T%TData)
		return
	end function
	function sumTensors(T)result(Res)
		real*4 ::Res
		class(Tensor),intent(in) :: T
		call sum_TDatas(Res,T%TData)
		return
	end function
	function sumTensord(T)result(Res)
		real*8 ::Res
		class(Tensor),intent(in) :: T
		call sum_TDatad(Res,T%TData)
		return
	end function
	function sumTensorc(T)result(Res)
		complex*8 ::Res
		class(Tensor),intent(in) :: T
		call sum_TDatac(Res,T%TData)
		return
	end function
	function sumTensorz(T)result(Res)
		complex*16 ::Res
		class(Tensor),intent(in) :: T
		call sum_TDataz(Res,T%TData)
		return
	end function
	function sumTensorT(T)result(Res)
		class(Tensor),allocatable ::Res
		class(Tensor),intent(in) :: T
		integer::tempi
		real*4::temps
		real*8::tempd
		complex*8::tempc
		complex*16::tempz
		allocate(Res,mold=T)
		select case (T%getType())
			case (1)
				call sum_TDatai(tempi,T%TData)
				Res=tempi
			case (2)
				call sum_TDatas(temps,T%TData)
				Res=temps
			case (3)
				call sum_TDatad(tempd,T%TData)
				Res=tempd
			case (4)
				call sum_TDatac(tempc,T%TData)
				Res=tempc
			case (5)
				call sum_TDataz(tempz,T%TData)
				Res=tempz
		end select
		return
	end function

	!***************************************************************************
	!***************************************************************************
	!
	!                                useful function
	!
	!***************************************************************************
	!**************************************************************************


	!***************   dot   ******************
	!		return <phi1|phi2>
	!	

	integer function idotTensor(phi1,phi2)result(dotTensor)
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		if(.not.phi1%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		if(.not.phi2%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		N1=phi1%TData%getTotalData()
		N2=phi2%TData%getTotalData()
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			call error_stop()
		end if
		call  product_dotc_int(dotTensor,phi1%TData,phi2%TData)	
		RETURN
	end function
	function sdotTensor(phi1,phi2)result(dotTensor)
		real(kind=4) ::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		if(.not.phi1%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		if(.not.phi2%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		N1=phi1%TData%getTotalData()
		N2=phi2%TData%getTotalData()
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			call error_stop()
		end if
		call  product_dotc_real4(dotTensor,phi1%TData,phi2%TData)	
		RETURN
	end function
	function ddotTensor(phi1,phi2)result(dotTensor)
		real(kind=8)::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		if(.not.phi1%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		if(.not.phi2%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		N1=phi1%TData%getTotalData()
		N2=phi2%TData%getTotalData()
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			call error_stop()
		end if
		call  product_dotc_real8(dotTensor,phi1%TData,phi2%TData)	
		RETURN
	end function
	function cdotTensor(phi1,phi2)result(dotTensor)
		complex(kind=4) ::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		if(.not.phi1%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		if(.not.phi2%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		N1=phi1%TData%getTotalData()
		N2=phi2%TData%getTotalData()
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			call error_stop()
		end if
		call  product_dotc_com4(dotTensor,phi1%TData,phi2%TData)	
		RETURN
	end function
	function zdotTensor(phi1,phi2)result(dotTensor)
		complex(kind=8) ::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		if(.not.phi1%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		if(.not.phi2%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		N1=phi1%TData%getTotalData()
		N2=phi2%TData%getTotalData()
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			call error_stop()
		end if
		call  product_dotc_com8(dotTensor,phi1%TData,phi2%TData)	
		RETURN
	end function
	function TdotTensor(phi1,phi2)result(dotTensor)
		class(Tensor),allocatable::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		integer::i
		allocate(dotTensor,mold=phi1)
		classtype=select_type_in_add_minu(phi1%TData,phi2%TData)
		select case(classtype)
			case(1)
				dotTensor=idotTensor(phi1,phi2)
			case(2)
				dotTensor=sdotTensor(phi1,phi2)
			case(3)
				dotTensor=ddotTensor(phi1,phi2)
			case(4)
				dotTensor=cdotTensor(phi1,phi2)
			case(5)
				dotTensor=zdotTensor(phi1,phi2)
			case default
				call writemess("ERROR in .x.",-1)
				call error_stop()
		end select
		RETURN
	end function


	!***************   dot   ******************
	!		return <phi1|phi2>, dot product DO NOT conjugating the first vector,The Tensor will be regard as a vector

	function idotUTensor(phi1,phi2)result(dotTensor)
		integer::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		if(.not.phi1%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		if(.not.phi2%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		N1=phi1%TData%getTotalData()
		N2=phi2%TData%getTotalData()
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			call error_stop()
		end if
		call  product_dot_int(dotTensor,phi1%TData,phi2%TData)	
		RETURN
	end function
	function sdotUTensor(phi1,phi2)result(dotTensor)
		real(kind=4) ::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		if(.not.phi1%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		if(.not.phi2%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		N1=phi1%TData%getTotalData()
		N2=phi2%TData%getTotalData()
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			call error_stop()
		end if
		call  product_dot_real4(dotTensor,phi1%TData,phi2%TData)	
		RETURN
	end function
	function ddotUTensor(phi1,phi2)result(dotTensor)
		real(kind=8) ::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		if(.not.phi1%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		if(.not.phi2%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		N1=phi1%TData%getTotalData()
		N2=phi2%TData%getTotalData()
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			call error_stop()
		end if
		call  product_dot_real8(dotTensor,phi1%TData,phi2%TData)	
		RETURN
	end function
	function cdotUTensor(phi1,phi2)result(dotTensor)
		complex(kind=4) ::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		if(.not.phi1%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		if(.not.phi2%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		N1=phi1%TData%getTotalData()
		N2=phi2%TData%getTotalData()
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			call error_stop()
		end if
		call  product_dot_com4(dotTensor,phi1%TData,phi2%TData)	
		RETURN
	end function
	function zdotUTensor(phi1,phi2)result(dotTensor)
		complex(kind=8) ::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		if(.not.phi1%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		if(.not.phi2%getflag())then
			write(*,*)"There is no data in the Tensor,(dot)"
			call error_stop()
		end if
		N1=phi1%TData%getTotalData()
		N2=phi2%TData%getTotalData()
		if(N1.ne.N2) then
			write(*,*)"ERROR in dot"
			write(*,*)"Total number of input is"
			write(*,*)N1,N2
			call error_stop()
		end if
		call  product_dot_com8(dotTensor,phi1%TData,phi2%TData)	
		RETURN
	end function
	
	function TdotUTensor(phi1,phi2)result(dotTensor)
		class(Tensor),allocatable::dotTensor
		class(Tensor),intent(in)::phi1
		class(Tensor),intent(in)::phi2
		integer::N1,N2,classtype
		integer::i
		allocate(dotTensor,mold=phi1)
		classtype=select_type_in_add_minu(phi1%TData,phi2%TData)
		select case(classtype)
			case(1)
				dotTensor=idotUTensor(phi1,phi2)
			case(2)
				dotTensor=sdotUTensor(phi1,phi2)
			case(3)
				dotTensor=ddotUTensor(phi1,phi2)
			case(4)
				dotTensor=cdotUTensor(phi1,phi2)
			case(5)
				dotTensor=zdotUTensor(phi1,phi2)
			case default
				call writemess("ERROR in .dot.",-1)
				call error_stop()
		end select
		RETURN
	end function
	
	
	!*****************  norm   *****************
	!		return  sqrt(<phi|phi>	)

	function inormTensor(T)result(normTensor)
		integer::normTensor
		class(Tensor),intent(in) :: T
		real*4::s
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		call norm_routine_real4(s,T%TData)
		normTensor=s
		return
	end function	
	function snormTensor(T)result(normTensor)
		real(kind=4) ::normTensor
		class(Tensor),intent(in) :: T
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		call norm_routine_real4(normTensor,T%TData)
		return
	end function	
	function dnormTensor(T)result(normTensor)
		real(kind=8) ::normTensor
		class(Tensor),intent(in) :: T
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		call norm_routine_real8(normTensor,T%TData)
		return
	end function	
	function cnormTensor(T)result(normTensor)
		complex(kind=4) ::normTensor
		class(Tensor),intent(in) :: T
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		call norm_routine_com4(normTensor,T%TData)
		return
	end function	
	function znormTensor(T)result(normTensor)
		complex(kind=8) ::normTensor
		class(Tensor),intent(in) :: T
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		call norm_routine_com8(normTensor,T%TData)
		return
	end function	
	function TnormTensor(T)result(normTensor)
		class(Tensor),allocatable::normTensor
		class(Tensor),intent(in) :: T
		integer::i
		real*4::s
		real*8::d
		complex*8::c
		complex*16::z
		allocate(normTensor,mold=T)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		select case(T%getType())
			case (1) 
				call norm_routine_real4(s,T%TData)
				normTensor=s
			case (2) 
				call norm_routine_real4(s,T%TData)
				normTensor=s
			case (3) 
				call norm_routine_real8(d,T%TData)
				normTensor=d
			case (4) 
				call norm_routine_com4(c,T%TData)
				normTensor=real(c,kind=4)
			case (5) 
				call norm_routine_com8(z,T%TData)
				normTensor=dble(z)
			case default
				call writemess("ERROR in norm",-1)
				call error_stop()
		end select
		return
	end function	

	!*****************  norm2   *****************
	!		return  <phi|phi>	

	function inorm2Tensor(T)result(normTensor)
		integer::normTensor
		class(Tensor),intent(in) :: T
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		call norm2_routine_int(normTensor,T%TData)
		return
	end function	
	function snorm2Tensor(T)result(normTensor)
		real(kind=4) ::normTensor
		class(Tensor),intent(in) :: T
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		normTensor=snormTensor(T)
		normTensor=normTensor*normTensor
		return
	end function	
	function dnorm2Tensor(T)result(normTensor)
		real(kind=8) ::normTensor
		class(Tensor),intent(in) :: T
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		normTensor=dnormTensor(T)
		normTensor=normTensor*normTensor
		return
	end function	
	function cnorm2Tensor(T)result(normTensor)
		complex(kind=4) ::normTensor
		class(Tensor),intent(in) :: T
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		normTensor=cnormTensor(T)
		normTensor=normTensor*normTensor
		return
	end function	
	function znorm2Tensor(T)result(normTensor)
		complex(kind=8) ::normTensor
		class(Tensor),intent(in) :: T
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		normTensor=znormTensor(T)
		normTensor=normTensor*normTensor
		return
	end function
	function Tnorm2Tensor(T)result(normTensor)
		class(Tensor),allocatable::normTensor
		class(Tensor),intent(in) :: T
		allocate(normTensor,mold=T)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(norm)"
			call error_stop()
		end if
		select case(T%getType())
			case (1)
				normTensor=inorm2Tensor(T)
			case default
				normTensor=TnormTensor(T)
				normTensor=normTensor*normTensor
		end select
		return
	end function

	!***************  directProduct  *********************

	function directProduct(T1,T2)
		class(Tensor),allocatable::directProduct
		class(Tensor),intent(in) :: T1
		class(Tensor),intent(in) :: T2
		integer :: m1,n1,m2,n2,classtype
		type(Dimension):: D1,Dtemp
		allocate(directProduct,mold=T1)
		if((T1%Dimension%getRank().eq.1).and.(T2%Dimension%getRank().eq.1)) then
			D1=T1%Dimension
			D1=D1+T2%Dimension
			classtype=select_type_in_add_minu(T1%TData,T2%TData)
			call directProduct%allocate(D1,classtype)
			n1=T1.dim.1
			n2=T2.dim.1
			call directProduct1_routine(directProduct%TDAta,T1%TData,T2%TData,n1,n2)
			return
		end if
		if((T1%Dimension%getRank().eq.2).and.(T2%Dimension%getRank().eq.2)) then
			D1=T1%Dimension.subdim.1
			Dtemp=T2%Dimension.subdim.1
			D1=D1+Dtemp
			Dtemp=T1%Dimension.subdim.2
			D1=D1+Dtemp
			Dtemp=T2%Dimension.subdim.2
			D1=D1+Dtemp
			m1=T1.dim.1
			n1=T1.dim.2
			m2=T2.dim.1
			n2=T2.dim.2
			classtype=select_type_in_add_minu(T1%TData,T2%TData)
			call directProduct%allocate(D1,classtype)
			call directProduct2_routine(directProduct%TData,T1%TData,T2%TData,m1,n1,m2,n2)
		else
			write(*,*)"The directProduct is just for matrix"
			call error_stop()
		end if
		return
	end function
	function directProductM(T1,T2)result(directProduct)!return a matrix
		class(Tensor),allocatable::directProduct
		class(Tensor),intent(in) :: T1
		class(Tensor),intent(in) :: T2
		integer :: m1,n1,m2,n2,classtype
		type(Dimension):: D1,Dtemp
		allocate(directProduct,mold=T1)
		if((T1%Dimension%getRank().eq.1).and.(T2%Dimension%getRank().eq.1)) then
			D1=T1%Dimension
			D1=D1+T2%Dimension
			call D1%fuse(1)
			classtype=select_type_in_add_minu(T1%TData,T2%TData)
			call directProduct%allocate(D1,classtype)
			n1=T1.dim.1
			n2=T2.dim.1
			call directProduct1_routine(directProduct%TDAta,T1%TData,T2%TData,n1,n2)
			return
		end if
		if((T1%Dimension%getRank().eq.2).and.(T2%Dimension%getRank().eq.2)) then
			D1=T1%Dimension.subdim.1
			Dtemp=T2%Dimension.subdim.1
			D1=D1+Dtemp
			Dtemp=T1%Dimension.subdim.2
			D1=D1+Dtemp
			Dtemp=T2%Dimension.subdim.2
			D1=D1+Dtemp
			call D1%fuse(1)
			call D1%fuse(2)
			m1=T1.dim.1
			n1=T1.dim.2
			m2=T2.dim.1
			n2=T2.dim.2
			classtype=select_type_in_add_minu(T1%TData,T2%TData)
			call directProduct%allocate(D1,classtype)
			call directProduct2_routine(directProduct%TData,T1%TData,T2%TData,m1,n1,m2,n2)
		else
			write(*,*)"The directProduct is just for matrix"
			call error_stop()
		end if
		return
	end function
	
	function directProductTensor(T1,T2)result(Res)
		class(Tensor),allocatable::Res
		class(Tensor),intent(in) :: T1
		class(Tensor),intent(in) :: T2
		type(Dimension)::newDim
		integer::classtype
		allocate(Res,mold=T1)
		newDim=T1%Dimension+T2%Dimension
		classtype=select_type_in_add_minu(T1%TData,T2%TData)
		call Res%allocate(newDim,classtype)
		call directProductTensor_routine(Res%TData,T1%TData,T2%TData)
		return
	end function


	!trace 

	function TtraceTensor(T)result(trace)
		class(Tensor),allocatable::trace
		class(Tensor),intent(in) :: T
		integer::rank,i
		allocate(trace,mold=T)
		if(T%TData%getTotalData().eq.1)then
			trace=T%i(1)
			return
		end if
		rank=T%Dimension%getRank()
		if(rank.ne.2) then
			call writemess("error in trace",-1)
			call writemess("input Tensor should be a matrix",-1)
			call error_stop()
		end if
		if((T.dim.1).ne.(T.dim.2)) then
			call writemess("error in trace",-1)
			call writemess("input Tensor should be a matrix",-1)
			call writemess((T.dim.1)+','+(T.dim.2),-1)
			call error_stop()
		end if
		select case(T%getType())
			case(1)
				trace=itraceTensor(T)
			case(2)
				trace=straceTensor(T)
			case(3)
				trace=dtraceTensor(T)
			case(4)
				trace=ctraceTensor(T)
			case(5)
				trace=ztraceTensor(T)
			case default
				call writemess("ERROR in trace",-1)
				call error_stop()
		end 	select
		return
	end function	
	function ztraceTensor(T)result(trace)
		complex*16 ::trace
		class(Tensor),intent(in) :: T
		integer::rank,i
		if(T%TData%getTotalData().eq.1)then
			trace=T.zi.1
			return
		end if
		rank=T%Dimension%getRank()
		if(rank.ne.2) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		end if
		if((T.dim.1).ne.(T.dim.2)) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			write(*,*)(T.dim.1),(T.dim.2)
			call error_stop()
		end if
		trace=dcmplx(0d0,0d0)
		do i=1,(T.dim.1)
			trace=trace+(T.zi.(/i,i/))
		end do
		return
	end function	
	function ctraceTensor(T)result(trace)
		complex*8 ::trace
		class(Tensor),intent(in) :: T
		integer::rank,i
		if(T%TData%getTotalData().eq.1)then
			trace=T.ci.1
			return
		end if
		rank=T%Dimension%getRank()
		if(rank.ne.2) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		end if
		if((T.dim.1).ne.(T.dim.2)) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			write(*,*)(T.dim.1),(T.dim.2)
			call error_stop()
		end if
		trace=cmplx(0d0,0d0)
		do i=1,(T.dim.1)
			trace=trace+(T.ci.(/i,i/))
		end do
		return
	end function	
	function dtraceTensor(T)result(trace)
		real(kind=8) ::trace
		class(Tensor),intent(in) :: T
		integer::rank,i
		if(T%TData%getTotalData().eq.1)then
			trace=T.di.1
			return
		end if
		rank=T%Dimension%getRank()
		if(rank.ne.2) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		end if
		if((T.dim.1).ne.(T.dim.2)) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			write(*,*)(T.dim.1),(T.dim.2)
			call error_stop()
		end if
		trace=0d0
		do i=1,(T.dim.1)
			trace=trace+(T.di.(/i,i/))
		end do
		return
	end function	
	function straceTensor(T)result(trace)
		real(kind=4) ::trace
		class(Tensor),intent(in) :: T
		integer::rank,i
		if(T%TData%getTotalData().eq.1)then
			trace=T.si.1
			return
		end if
		rank=T%Dimension%getRank()
		if(rank.ne.2) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		end if
		if((T.dim.1).ne.(T.dim.2)) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			write(*,*)(T.dim.1),(T.dim.2)
			call error_stop()
		end if
		trace=0.0
		do i=1,(T.dim.1)
			trace=trace+(T.si.(/i,i/))
		end do
		return
	end function	
	function itraceTensor(T)result(trace)
		integer::trace
		class(Tensor),intent(in) :: T
		integer::rank,i
		if(T%TData%getTotalData().eq.1)then
			trace=T.ii.1
			return
		end if
		rank=T%Dimension%getRank()
		if(rank.ne.2) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		end if
		if((T.dim.1).ne.(T.dim.2)) then
			write(*,*)"error in trace"
			write(*,*)"input Tensor should be a matrix"
			write(*,*)(T.dim.1),(T.dim.2)
			call error_stop()
		end if
		trace=0
		do i=1,(T.dim.1)
			trace=trace+(T.ii.(/i,i/))
		end do
		return
	end function	

	!*****************  Htranspose  *****************

	function Htranspose(T)
		class(Tensor),allocatable::Htranspose
		class(Tensor),intent(in) :: T
		integer::rank,m,n,i
		integer,allocatable::indices(:)
		type(dimension)::dimen
		allocate(Htranspose,mold=T)
		rank=T%Dimension%getRank()
		if(rank.eq.1) then
			m=T.dim. 1
			call Htranspose%allocate(T,T%getType())
			call conjugate_subroutine(Htranspose%TData,T%TData)
		else if(rank.eq.2) then
			m=T.dim.1
			n=T.dim.2
			dimen=T%Dimension%Dimpermute((/2,1/))
			call Htranspose%allocate( dimen,T%getType())
			call Htranspose_subroutine(Htranspose%TData,T%TData,m,n)
		else
			allocate(indices(rank))
			do i=1,rank
				indices(rank-i+1)=i
			end do
			Htranspose=(.con.T).p.indices
		end if
		return
	end function
	function Htranspose2(T)result(Htranspose)
		class(Tensor),allocatable::Htranspose
		class(Tensor),intent(in) :: T
		integer :: rank,i,charlen
		integer,allocatable::indices(:)
		character(len=max_len_of_char_in_TData)::Tname
		allocate(Htranspose,mold=T)
		rank=T%Dimension%getRank()
		allocate(indices(rank))
		do i=1,rank
			indices(i)=rank-i+1
		end do
		Htranspose=.con.(T.p.indices)
		do i=1,rank
			Tname=Htranspose%Dimension%getName(i)
			charlen=len(trim(Tname))
			if(Tname(charlen:charlen).eq.dag_mark) then
				call Htranspose%setName(i,Tname(1:charlen-1))
			else
				call Htranspose%setName(i,Tname+dag_mark)
			end if
		end do
		return
	end function

	!*****************  conjugate  *****************

	function conjugate(T)
		class(Tensor),allocatable::conjugate
		class(Tensor),intent(in) :: T
		complex*16,allocatable :: newdata(:)
		integer :: m
		allocate(conjugate,mold=T)
		call conjugate%allocate(T,T%getType())
		call conjugate_subroutine(conjugate%TData,T%TData)
		return
	end function

	!transpose

	function transposeTensor(T)
		class(Tensor),allocatable::transposeTensor
		class(Tensor),intent(in) :: T
		integer::rank,m,n
		type(dimension)::dimen
		allocate(transposeTensor,mold=T)
		rank=T%Dimension%getRank()
		if(rank.eq.1) then
			transposeTensor=T
		else if(rank.eq.2) then
			m=T.dim.1
			n=T.dim.2
			call transposeTensor%allocate( T%Dimension%Dimpermute((/2,1/)) ,T%getType())
			call transpose_subroutine(transposeTensor%TData,T%TData,m,n)
		else
			write(*,*) "Tensor should be 1 or 2 dimension"
			write(*,*) "ERROR in Htranspose,stop"
			call error_stop()
		end if
		return
	end function

	!cccccccccccccccc	pauli_matrix  ccccccccccccccccccccccc
	!	the output is sigma_i*num,where i=x,y,z

	subroutine pauli_matrix(Tx,Ty,Tz,num)
		class(Tensor) :: Tx,Ty,Tz
		class(*),optional,intent(in)::num
		complex*16::numr
		complex*16 :: x(2,2)
		complex*16 :: y(2,2)
		complex*16 :: z(2,2)
		complex*16::EE=(1d0,0d0)
		complex*16::II=(0d0,1d0)
		if(present(num)) then
			numr=dselect(num)
		else
			numr=dcmplx(1d0,0d0)
		end if
		x(1,1)=dcmplx(0,0)
		x(1,2)=EE*numr
		x(2,1)=EE*numr
		x(2,2)=dcmplx(0,0)

		y(1,1)=dcmplx(0,0)
		y(1,2)=-1d0*II*numr
		y(2,1)=II*numr
		y(2,2)=dcmplx(0,0)

		z(1,1)=EE*numr
		z(1,2)=dcmplx(0,0)
		z(2,1)=dcmplx(0,0)
		z(2,2)=-1d0*EE*numr
		Tx=x
		Ty=y
		Tz=z
		return
	end subroutine	

	! return exp(H)

	 function expmTensor(H)
	 	class(Tensor),allocatable::expmTensor
		class(Tensor),intent(in) ::H
		type(Tensor)::temp
		real*8::a
		integer::i
		allocate(expmTensor,mold=H)
		if(H%Dimension%getRank().ne.2) then
			call writemess("ERROR in expm",-1)
			call writemess("input Tensor should be a matrix",-1)
			call error_stop()
		end if
		if(H%dim(1).ne.H%dim(2)) then
			call writemess("ERROR in expm",-1)
			call error_stop()
		end if
		call expmTensor%setType(H%getType())
		temp=H
		expmTensor=H
		a=1
		do i=2,99999
			temp=temp*H
			if(temp%isZero())exit
			temp=temp/i
			expmTensor=expmTensor+temp
		end do
		expmTensor=expmTensor+eye(H%dim(1),H%dim(2))
		return

	end function
	!H = vec*eye(val)*(vec**H)

	subroutine eigvalue(H,val,vec)
		class(Tensor),intent(in) ::H
		class(Tensor),intent(inout) ::val
		class(Tensor),optional,intent(inout)::vec
		integer::hdim(2)
		hdim(1)=H.dim.1
		hdim(2)=H.dim.2
		if(H%Dimension%getRank().ne.2) then
			write(*,*)"ERROR in eng"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		end if
		if(hdim(1).ne.hdim(2)) then
			write(*,*)"ERROR in eng"
			call error_stop()
		end if
		call val%empty()
		select case(H%getType())
			case (5)
				call val%allocate((/hdim(1)/),5)
				if(present(vec))then
					call vec%empty()
					call vec%allocate(H%Dimension,5)
				end if
			case (4)
				call val%allocate((/hdim(1)/),4)
				if(present(vec))then
					call vec%empty()
					call vec%allocate(H%Dimension,4)
				end if
			case (3)
				call val%allocate((/hdim(1)/),5)
				if(present(vec))then
					call vec%empty()
					call vec%allocate(H%Dimension,3)
				end if
			case (2)
				call val%allocate((/hdim(1)/),4)
				if(present(vec))then
					call vec%empty()
					call vec%allocate(H%Dimension,2)
				end if
			case (1)
				call val%allocate((/hdim(1)/),4)
				if(present(vec))then
					call vec%empty()
					call vec%allocate(H%Dimension,2)
				end if
		end select
		call eigenvalue_TData_routine(H%TData,hdim(1),val%TData,vec%TData)
		return
	end subroutine	
	
	function eigTensor(T,outvex)	result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),target,intent(in)::T
		logical,optional,intent(in)::outvex
		if(present(outvex).and.outvex)then
			allocate(res(2),mold=T)
			call eigvalue(T,res(1),res(2))
		else
			allocate(res(1),mold=T)
			call eigvalue(T,res(1))
		end if
		return
	end function

	!***************************************************************************
	!***************************************************************************
	!
	!                             eye tensor
	!
	!***************************************************************************
	!**************************************************************************

	function eye_int(s,dim1,dim2)
		class(Tensor),allocatable::eye_int
		integer,intent(in) :: s(:)
		integer,optional,intent(in) :: dim1,dim2
		integer :: i,m,n,lens
		allocate(Tensor::eye_int)
		lens=size(s)
		if(present(dim1))then
			m=dim1
			n=dim2
		else
			m=lens
			n=lens
		end if
		call eye_int%allocate((/m,n/),1)
		call eye_int%zero()
		do i=1,min(m,n,lens)
			call eye_int%setValue((/i,i/),s(i))
		end do
		return
	end function
	function eye_real4(s,dim1,dim2)result(eye)
		class(Tensor),allocatable::eye
		real(kind=4),intent(in) :: s(:)
		integer,optional,intent(in) :: dim1,dim2
		integer :: i,m,n,lens
		allocate(Tensor::eye)
		lens=size(s)
		if(present(dim1))then
			m=dim1
			n=dim2
		else
			m=lens
			n=lens
		end if
		call eye%allocate((/m,n/),2)
		call eye%zero()
		do i=1,min(m,n,lens)
			call eye%setValue((/i,i/),s(i))
		end do
		return
	end function
	function eye_real8(s,dim1,dim2)result(eye)
		class(Tensor),allocatable::eye
		real(kind=8),intent(in) :: s(:)
		integer,optional,intent(in) :: dim1,dim2
		integer :: i,m,n,lens
		allocate(Tensor::eye)
		lens=size(s)
		if(present(dim1))then
			m=dim1
			n=dim2
		else
			m=lens
			n=lens
		end if
		call eye%allocate((/m,n/),3)
		call eye%zero()
		do i=1,min(m,n,lens)
			call eye%setValue((/i,i/),s(i))
		end do
		return
	end function
	function eye_com4(s,dim1,dim2)result(eye)
		class(Tensor),allocatable::eye
		complex(kind=4),intent(in) :: s(:)
		integer,optional,intent(in) :: dim1,dim2
		integer :: i,m,n,lens
		allocate(Tensor::eye)
		lens=size(s)
		if(present(dim1))then
			m=dim1
			n=dim2
		else
			m=lens
			n=lens
		end if
		call eye%allocate((/m,n/),4)
		call eye%zero()
		do i=1,min(m,n,lens)
			call eye%setValue((/i,i/),s(i))
		end do
		return
	end function
	function eye_com8(s,dim1,dim2)result(eye)
		class(Tensor),allocatable::eye
		complex(kind=8),intent(in) :: s(:)
		integer,optional,intent(in) :: dim1,dim2
		integer :: i,m,n,lens
		allocate(Tensor::eye)
		lens=size(s)
		if(present(dim1))then
			m=dim1
			n=dim2
		else
			m=lens
			n=lens
		end if
		call eye%allocate((/m,n/),5)
		call eye%zero()
		do i=1,min(m,n,lens)
			call eye%setValue((/i,i/),s(i))
		end do
		return
	end function
	function eye_logi(s,dim1,dim2)result(eye)
		class(Tensor),allocatable::eye
		logical,intent(in) :: s(:)
		integer,optional,intent(in) :: dim1,dim2
		integer :: i,m,n,lens
		allocate(Tensor::eye)
		lens=size(s)
		if(present(dim1))then
			m=dim1
			n=dim2
		else
			m=lens
			n=lens
		end if
		call eye%allocate((/m,n/),6)
		call eye%zero()
		do i=1,min(m,n,lens)
			call eye%setValue((/i,i/),s(i))
		end do
		return
	end function
	function eye_char(s,dim1,dim2)result(eye)
		class(Tensor),allocatable::eye
		character(len=*),intent(in) :: s(:)
		integer,optional,intent(in) :: dim1,dim2
		integer :: i,m,n,lens
		allocate(Tensor::eye)
		lens=size(s)
		if(present(dim1))then
			m=dim1
			n=dim2
		else
			m=lens
			n=lens
		end if
		call eye%allocate((/m,n/),7)
		call eye%zero()
		do i=1,min(m,n,lens)
			call eye%setValue((/i,i/),s(i))
		end do
		return
	end function
	function eye_Tensor(s,dim1,dim2)result(eye)
		class(Tensor),allocatable::eye
		class(Tensor),intent(in) :: s
		integer,optional,intent(in) :: dim1,dim2
		integer :: i,m,n,lens
		allocate(Tensor::eye)
		lens=s%TData%getTotalData()
		if(s%Dimension%getRank().ne.1)then
			call writemess('ERROR in eye(T), input T should be a vector(rank=1)')
			call s%diminfo()
			call error_stop
		end if
		if(present(dim1))then
			m=dim1
			n=dim2
		else
			m=lens
			n=lens
		end if
		call eye%allocate((/m,n/),s%getType())
		call eye%zero()
		do i=1,min(m,n,lens)
			call eye%setValue((/i,i/),s%zi((/i/)))
		end do
		return
	end function
	function diag_int(s,m,n)result(diag)
		class(Tensor),allocatable::diag
		integer,intent(in) :: s
		integer,intent(in) :: m,n
		integer :: i
		allocate(Tensor::diag)
		call diag%allocate((/m,n/),1)
		call diag%zero()
		do i=1,min(m,n)
			call diag%setValue((/i,i/),s)
		end do
		return
	end function
	function identity_matrix_Tensor(m,n)result(diag)
		class(Tensor),allocatable::diag
		integer,intent(in) :: m,n
		integer :: i
		allocate(Tensor::diag)
		call diag%allocate((/m,n/),1)
		call diag%zero()
		do i=1,min(m,n)
			call diag%setValue((/i,i/),1)
		end do
		return
	end function
	function diag_real4(s,m,n)result(diag)
		class(Tensor),allocatable::diag
		real(kind=4),intent(in) :: s
		integer,intent(in) :: m,n
		integer :: i
		allocate(Tensor::diag)
		call diag%allocate((/m,n/),2)
		call diag%zero()
		do i=1,min(m,n)
			call diag%setValue((/i,i/),s)
		end do
		return
	end function
	function diag_real8(s,m,n)result(diag)
		class(Tensor),allocatable::diag
		real(kind=8),intent(in) :: s
		integer,intent(in) :: m,n
		integer :: i
		allocate(Tensor::diag)
		call diag%allocate((/m,n/),3)
		call diag%zero()
		do i=1,min(m,n)
			call diag%setValue((/i,i/),s)
		end do
		return
	end function
	function diag_com4(s,m,n)result(diag)
		class(Tensor),allocatable::diag
		complex(kind=4),intent(in) :: s
		integer,intent(in) :: m,n
		integer :: i
		allocate(Tensor::diag)
		call diag%allocate((/m,n/),4)
		call diag%zero()
		do i=1,min(m,n)
			call diag%setValue((/i,i/),s)
		end do
		return
	end function
	function diag_com8(s,m,n)result(diag)
		class(Tensor),allocatable::diag
		complex(kind=8),intent(in) :: s
		integer,intent(in) :: m,n
		integer :: i
		allocate(Tensor::diag)
		call diag%allocate((/m,n/),5)
		call diag%zero()
		do i=1,min(m,n)
			call diag%setValue((/i,i/),s)
		end do
		return
	end function
	function diag_logi(s,m,n)result(diag)
		class(Tensor),allocatable::diag
		logical,intent(in) :: s
		integer,intent(in) :: m,n
		integer :: i
		allocate(Tensor::diag)
		call diag%allocate((/m,n/),6)
		call diag%zero()
		do i=1,min(m,n)
			call diag%setValue((/i,i/),s)
		end do
		return
	end function
	function diag_char(s,m,n)result(diag)
		class(Tensor),allocatable::diag
		character(len=*),intent(in) :: s
		integer,intent(in) :: m,n
		integer :: i
		allocate(Tensor::diag)
		call diag%allocate((/m,n/),7)
		call diag%zero()
		do i=1,min(m,n)
			call diag%setValue((/i,i/),s)
		end do
		return
	end function
	function diag_classtype(m,n,classtype)result(diag)
		class(Tensor),allocatable::diag
		integer,intent(in) :: m,n,classtype
		integer :: i
		allocate(Tensor::diag)
		call diag%allocate((/m,n/),classtype)
		call diag%zero()
		do i=1,min(m,n)
			call diag%setValue((/i,i/),1)
		end do
		return
	end function
	function diag_type(m,n,classtype_)result(diag)
		class(Tensor),allocatable::diag
		integer,intent(in) :: m,n
		character(len=*),intent(in)::classtype_
		integer :: i,classtype
		allocate(Tensor::diag)
		classtype=select_data_type_char(classtype_)
		call diag%allocate((/m,n/),classtype)
		call diag%zero()
		do i=1,min(m,n)
			call diag%setValue((/i,i/),1)
		end do
		return
	end function

	subroutine eye_Tensor2(T)
		class(Tensor),intent(inout)::T
		integer::i,rank
		rank=T%Dimension%getRank()
		if(rank.ne.2)then
			call writemess("ERROR in set eye in Tensor,input should be a matrix",-1)
			call error_stop()
		end if
		call T%zero()
		do i=1,min(T%dim(1),T%dim(2))
			call T%setvalue((/i,i/),1)
		end do
		return
	end subroutine
	subroutine eye_Tensor3(T,m,n,classtype)
		class(Tensor),intent(inout)::T
		integer,intent(in)::m,n
		character(len=*),intent(in)::classtype
		integer::i
		call T%empty()
		call T%allocate((/m,n/),classtype)
		call T%zero()
		do i=1,min(T%dim(1),T%dim(2))
			call T%setvalue((/i,i/),1)
		end do
		return
	end subroutine
	subroutine eye_Tensor4(T,m,n,classtype)
		class(Tensor),intent(inout)::T
		integer,intent(in)::m,n,classtype
		integer::i
		call T%empty()
		call T%allocate((/m,n/),classtype)
		call T%zero()
		do i=1,min(T%dim(1),T%dim(2))
			call T%setvalue((/i,i/),1)
		end do
		return
	end subroutine
	subroutine eye_Tensor5(T,dimen,classtype)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::classtype
		integer,intent(in)::dimen(2)
		call T%empty()
		call T%allocate(dimen,classtype)
		call eye_Tensor2(T)
		return
	end subroutine
	subroutine eye_Tensor6(T,dimen,classtype)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::classtype
		type(dimension),intent(in)::dimen
		call T%empty()
		call T%allocate(dimen,classtype)
		call eye_Tensor2(T)
		return
	end subroutine

	subroutine sortTensor1(T)
		class(Tensor),intent(inout)::T
		logical::increase,realpart
		increase=.true.
		realpart=.true.
		call sortTData2(T%TData,increase,realpart)
		return
	end subroutine
	subroutine sortTensor2(T,indices)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(inout)::indices
		logical::increase,realpart
		increase=.true.
		realpart=.true.
		call indices%empty()
		call indices%allocate([T%TData%getTotalData()],'integer')
		call sortTData2(T%TData,increase,realpart)
		return
	end subroutine
	
	subroutine sortTensor3(T,increase)
		class(Tensor),intent(inout)::T
		logical,intent(in)::increase
		logical::realpart
		realpart=.true.
		call sortTData2(T%TData,increase,realpart)
		return
	end subroutine
	subroutine sortTensor4(T,indices,increase)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(inout)::indices
		logical,intent(in)::increase
		logical::realpart
		realpart=.true.
		call indices%empty()
		call indices%allocate([T%TData%getTotalData()],'integer')
		call sortTData2(T%TData,increase,realpart)
		return
	end subroutine
	
	subroutine sortTensor5(T,increase,realpart)
		class(Tensor),intent(inout)::T
		logical,intent(in)::increase,realpart
		call sortTData2(T%TData,increase,realpart)
		return
	end subroutine
	subroutine sortTensor6(T,indices,increase,realpart)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(inout)::indices
		logical,intent(in)::increase,realpart
		call indices%empty()
		call indices%allocate([T%TData%getTotalData()],'integer')
		call sortTData2(T%TData,increase,realpart)
		return
	end subroutine
	
	

	integer function iwhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		integer,intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (1)
				do whichindex=1,T%TData%getTotalData()
					if(element.eq.T%ii(whichindex))return
				end do
			case default
				call writemess('The Tensor is integer, input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	integer function swhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		real*4,intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (2)
				do whichindex=1,T%TData%getTotalData()
					if(abs(element-T%si(whichindex)).le. default_zero_real_number)return
				end do
			case default
				call writemess('The Tensor is real, input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	
	integer function dwhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		real*8,intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (3)
				do whichindex=1,T%TData%getTotalData()
					if(dabs(element-T%di(whichindex)).le. default_zero_double_number)return
				end do
			case default
				call writemess('The Tensor is real*8, input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	
	integer function cwhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		complex*8,intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (4)
				do whichindex=1,T%TData%getTotalData()
					if(cabs(element-T%ci(whichindex)).le. default_zero_real_number)return
				end do
			case default
				call writemess('The Tensor is complex(kind=4), input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	
	integer function zwhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		complex*16,intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (5)
				do whichindex=1,T%TData%getTotalData()
					if(cdabs(element-T%zi(whichindex)).le. default_zero_double_number)return
				end do
			case default
				call writemess('The Tensor is complex(kind=4), input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	integer function awhichindex(T,element)result(whichindex)
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (7)
				do whichindex=1,T%TData%getTotalData()
					if(element.equ.T%ai(whichindex))return
				end do
			case default
				call writemess('The Tensor is character, input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function
	
	integer function awhichindex2(T,element,maxlen)result(whichindex)
		class(Tensor),intent(in)::T
		integer,intent(in)::maxlen
		character(len=*),intent(in)::element
		if(.not.T%getflag())then
			whichindex=0
			return
		end if
		select case(T%getType())
			case (7)
				do whichindex=1,maxlen
					if(element.equ.T%ai(whichindex))return
				end do
			case default
				call writemess('The Tensor is character, input type error, (which)',-1)
				call error_stop
		end select
		whichindex=0
		return
	end function


	function enlargeTensorReal8(T,ithleg,newD,randomNumberScal)result(enlargeTensor)
		class(Tensor),allocatable::enlargeTensor
		class(Tensor),intent(in)::T
		integer,intent(in)::newD,ithleg
		real*8,intent(in)::randomNumberScal
		integer::D,rank
		type(Dimension)::dimen
		allocate(enlargeTensor,mold=T)
		rank=T%Dimension%getRank()
		D=T%dim(ithleg)
		if(D.eq.newD)then
			enlargeTensor=T
			return
		end if
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.T%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,T%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if	
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call enlargeTensor%allocate(dimen,T%getType())
		if(T%ifDynamic())call enlargeTensor%Dynamic()
		call enlargeTensor%random([-randomNumberScal,randomNumberScal])
		call enlargeTensor%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		enlargeTensor=enlargeTensor.pbi.ithleg
		return
	end function

	function enlargeTensorReal8_Name(T,nameleg,newD,randomNumberScal)result(enlargeTensor)
		class(Tensor),allocatable::enlargeTensor
		class(Tensor),intent(in)::T
		integer,intent(in)::newD
		character(len=*),intent(in)::nameleg
		real*8,intent(in)::randomNumberScal
		integer::D,rank,ithleg
		type(Dimension)::dimen
		allocate(enlargeTensor,mold=T)
		rank=T%Dimension%getRank()
		ithleg=T%FindOrder(nameleg)
		D=T%dim(ithleg)
		if(D.eq.newD)then
			enlargeTensor=T
			return
		end if
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.T%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,T%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if	
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call enlargeTensor%allocate(dimen,T%getType())
		if(T%ifDynamic())call enlargeTensor%Dynamic()
		call enlargeTensor%random([-randomNumberScal,randomNumberScal])
		call enlargeTensor%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		enlargeTensor=enlargeTensor.pbi.ithleg
		return
	end function

	function enlargeTensorAllReal8(T,newD,randomNumberScal)result(enlargeTensor)
		class(Tensor),allocatable::enlargeTensor
		class(Tensor),intent(in)::T
		integer,intent(in)::newD(:)
		real*8,intent(in)::randomNumberScal
		integer::rank,i
		allocate(enlargeTensor,mold=T)
		rank=T%Dimension%getRank()
		enlargeTensor=enlargeTensorReal8(T,1,newD(1),randomNumberScal)
		do i=2,rank
			enlargeTensor=enlargeTensorReal8(enlargeTensor,i,newD(i),randomNumberScal)
		end do
		return
	end function
	subroutine enlargeTensorReal8Routine(inoutT,ithleg,newD,randomNumberScal)
		class(Tensor),intent(inout)::inoutT
		integer,intent(in)::newD,ithleg
		real*8,intent(in)::randomNumberScal
		integer::D,rank
		type(Dimension)::dimen
		class(Tensor),allocatable::T
		logical::ifDynamic
		allocate(T,mold=inoutT)
		D=inoutT%dim(ithleg)
		if(D.eq.newD)return
		T=inoutT
		ifDynamic=inoutT%ifDynamic()
		rank=inoutT%Dimension%getRank()
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.inoutT%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.inoutT%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,inoutT%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if	
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call inoutT%empty()
		call inoutT%allocate(dimen,T%getType())
		if(ifDynamic)call inoutT%Dynamic()
		call inoutT%random([-randomNumberScal,randomNumberScal])
		call inoutT%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		inoutT=inoutT.pbi.ithleg
		return
	end subroutine

	subroutine enlargeTensorReal8_NameRoutine(inoutT,Nameleg,newD,randomNumberScal)
		class(Tensor),intent(inout)::inoutT
		character(len=*),intent(in)::Nameleg
		integer,intent(in)::newD
		real*8,intent(in)::randomNumberScal
		integer::D,rank,ithleg
		type(Dimension)::dimen
		class(Tensor),allocatable::T
		logical::ifDynamic
		allocate(T,mold=inoutT)
		ithleg=inoutT%FindOrder(nameleg)
		D=inoutT%dim(ithleg)
		if(D.eq.newD)return
		T=inoutT
		ifDynamic=inoutT%ifDynamic()
		rank=inoutT%Dimension%getRank()
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.inoutT%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.inoutT%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,inoutT%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if		
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call inoutT%empty()
		call inoutT%allocate(dimen,T%getType())
		if(ifDynamic)call inoutT%Dynamic()
		call inoutT%random([-randomNumberScal,randomNumberScal])
		call inoutT%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		inoutT=inoutT.pbi.ithleg
		return
	end subroutine

	subroutine enlargeTensorAllReal8Routine(inoutT,newD,randomNumberScal)
		class(Tensor),intent(inout)::inoutT
		integer,intent(in)::newD(:)
		real*8,intent(in)::randomNumberScal
		integer::rank,i
		rank=inoutT%Dimension%getRank()
		call enlargeTensorReal8Routine(inoutT,1,newD(1),randomNumberScal)
		do i=2,rank
			call enlargeTensorReal8Routine(inoutT,i,newD(i),randomNumberScal)
		end do
		return
	end subroutine
	function enlargeTensorReal4(T,ithleg,newD,randomNumberScal)result(enlargeTensor)
		class(Tensor),allocatable::enlargeTensor
		class(Tensor),intent(in)::T
		integer,intent(in)::newD,ithleg
		real*4,intent(in)::randomNumberScal
		integer::D,rank
		type(Dimension)::dimen
		allocate(enlargeTensor,mold=T)
		rank=T%Dimension%getRank()
		D=T%dim(ithleg)
		if(D.eq.newD)then
			enlargeTensor=T
			return
		end if
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.T%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,T%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if	
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call enlargeTensor%allocate(dimen,T%getType())
		if(T%ifDynamic())call enlargeTensor%Dynamic()
		call enlargeTensor%random([-randomNumberScal,randomNumberScal])
		call enlargeTensor%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		enlargeTensor=enlargeTensor.pbi.ithleg
		return
	end function
	function enlargeTensorReal4_name(T,nameleg,newD,randomNumberScal)result(enlargeTensor)
		class(Tensor),allocatable::enlargeTensor
		class(Tensor),intent(in)::T
		integer,intent(in)::newD
		character(len=*),intent(in)::nameleg
		real*4,intent(in)::randomNumberScal
		integer::D,rank,ithleg
		type(Dimension)::dimen
		allocate(enlargeTensor,mold=T)
		rank=T%Dimension%getRank()
		ithleg=T%FindOrder(nameleg)
		D=T%dim(ithleg)
		if(D.eq.newD)then
			enlargeTensor=T
			return
		end if
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.T%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,T%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if		
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call enlargeTensor%allocate(dimen,T%getType())
		if(T%ifDynamic())call enlargeTensor%Dynamic()
		call enlargeTensor%random([-randomNumberScal,randomNumberScal])
		call enlargeTensor%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		enlargeTensor=enlargeTensor.pbi.ithleg
		return
	end function
	function enlargeTensorAllReal4(T,newD,randomNumberScal)result(enlargeTensor)
		class(Tensor),allocatable::enlargeTensor
		class(Tensor),intent(in)::T
		integer,intent(in)::newD(:)
		real*4,intent(in)::randomNumberScal
		integer::rank,i
		allocate(enlargeTensor,mold=T)
		rank=T%Dimension%getRank()
		enlargeTensor=enlargeTensorReal4(T,1,newD(1),randomNumberScal)
		do i=2,rank
			enlargeTensor=enlargeTensorReal4(enlargeTensor,i,newD(i),randomNumberScal)
		end do
		return
	end function
	subroutine enlargeTensorReal4Routine(inoutT,ithleg,newD,randomNumberScal)
		class(Tensor),intent(inout)::inoutT
		integer,intent(in)::newD,ithleg
		real*4,intent(in)::randomNumberScal
		integer::D,rank
		type(Dimension)::dimen
		class(Tensor),allocatable::T
		logical::ifDynamic
		allocate(T,mold=inoutT)
		D=inoutT%dim(ithleg)
		if(D.eq.newD)return
		T=inoutT
		ifDynamic=inoutT%ifDynamic()
		rank=T%Dimension%getRank()
		
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.inoutT%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.inoutT%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,inoutT%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if	
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call inoutT%empty()
		call inoutT%allocate(dimen,T%getType())
		if(ifDynamic)call inoutT%Dynamic()
		call inoutT%random([-randomNumberScal,randomNumberScal])
		call inoutT%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		inoutT=inoutT.pbi.ithleg
		return
	end subroutine
	subroutine enlargeTensorReal4_nameRoutine(inoutT,nameleg,newD,randomNumberScal)
		class(Tensor),intent(inout)::inoutT
		integer,intent(in)::newD
		character(len=*),intent(in)::nameleg
		real*4,intent(in)::randomNumberScal
		integer::D,rank,ithleg
		type(Dimension)::dimen
		class(Tensor),allocatable::T
		logical::ifDynamic
		allocate(T,mold=inoutT)
		ithleg=inoutT%FindOrder(nameleg)
		D=inoutT%dim(ithleg)
		if(D.eq.newD)return
		T=inoutT
		ifDynamic=inoutT%ifDynamic()
		rank=T%Dimension%getRank()
		
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.inoutT%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.inoutT%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,inoutT%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if		
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call inoutT%empty()
		call inoutT%allocate(dimen,T%getType())
		if(ifDynamic)call inoutT%Dynamic()
		call inoutT%random([-randomNumberScal,randomNumberScal])
		call inoutT%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		inoutT=inoutT.pbi.ithleg
		return
	end subroutine
	subroutine enlargeTensorAllReal4Routine(inoutT,newD,randomNumberScal)
		class(Tensor),intent(inout)::inoutT
		integer,intent(in)::newD(:)
		real*4,intent(in)::randomNumberScal
		integer::rank,i
		rank=inoutT%Dimension%getRank()
		call enlargeTensorReal4Routine(inoutT,1,newD(1),randomNumberScal)
		do i=2,rank
			call enlargeTensorReal4Routine(inoutT,i,newD(i),randomNumberScal)
		end do
		return
	end subroutine
	function enlargeTensorInt(T,ithleg,newD,randomNumberScal)result(enlargeTensor)
		class(Tensor),allocatable::enlargeTensor
		class(Tensor),intent(in)::T
		integer,intent(in)::newD,ithleg
		integer,intent(in)::randomNumberScal
		integer::D,rank
		type(Dimension)::dimen
		rank=T%Dimension%getRank()
		allocate(enlargeTensor,mold=T)
		D=T%dim(ithleg)
		if(D.eq.newD)then
			enlargeTensor=T
			return
		end if
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.T%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,T%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if	
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call enlargeTensor%allocate(dimen,T%getType())
		if(T%ifDynamic())call enlargeTensor%Dynamic()
		call enlargeTensor%random([-randomNumberScal,randomNumberScal])
		call enlargeTensor%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		enlargeTensor=enlargeTensor.pbi.ithleg
		return
	end function
	function enlargeTensorInt_name(T,nameleg,newD,randomNumberScal)result(enlargeTensor)
		class(Tensor),allocatable::enlargeTensor
		class(Tensor),intent(in)::T
		integer,intent(in)::newD
		character(len=*),intent(in)::nameleg
		integer,intent(in)::randomNumberScal
		integer::D,rank,ithleg
		type(Dimension)::dimen
		allocate(enlargeTensor,mold=T)
		rank=T%Dimension%getRank()
		ithleg=T%findOrder(nameleg)
		D=T%dim(ithleg)
		if(D.eq.newD)then
			enlargeTensor=T
			return
		end if
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.T%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,T%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if		
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call enlargeTensor%allocate(dimen,T%getType())
		if(T%ifDynamic())call enlargeTensor%Dynamic()
		call enlargeTensor%random([-randomNumberScal,randomNumberScal])
		call enlargeTensor%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		enlargeTensor=enlargeTensor.pbi.ithleg
		return
	end function
	function enlargeTensorAllInt(T,newD,randomNumberScal)result(enlargeTensor)
		class(Tensor),allocatable::enlargeTensor
		class(Tensor),intent(in)::T
		integer,intent(in)::newD(:)
		integer,intent(in)::randomNumberScal
		integer::rank,i
		allocate(enlargeTensor,mold=T)
		rank=T%Dimension%getRank()
		enlargeTensor=enlargeTensorInt(T,1,newD(1),randomNumberScal)
		do i=2,rank
			enlargeTensor=enlargeTensorInt(enlargeTensor,i,newD(i),randomNumberScal)
		end do
		return
	end function
	subroutine enlargeTensorIntRoutine(inoutT,ithleg,newD,randomNumberScal)
		class(Tensor),intent(inout)::inoutT
		integer,intent(in)::newD,ithleg
		integer,intent(in)::randomNumberScal
		integer::D,rank
		type(Dimension)::dimen
		class(Tensor),allocatable::T
		logical::ifDynamic
		allocate(T,mold=inoutT)
		D=inoutT%dim(ithleg)
		if(D.eq.newD)return
		T=inoutT
		ifDynamic=inoutT%ifDynamic()
		rank=T%Dimension%getRank()
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.inoutT%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.inoutT%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,inoutT%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if	
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call inoutT%empty()
		call inoutT%allocate(dimen,T%getType())
		if(ifDynamic)call inoutT%Dynamic()
		call inoutT%random([-randomNumberScal,randomNumberScal])
		call inoutT%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		inoutT=inoutT.pbi.ithleg
		return
	end subroutine
	subroutine enlargeTensorInt_nameRoutine(inoutT,nameleg,newD,randomNumberScal)
		class(Tensor),intent(inout)::inoutT
		integer,intent(in)::newD
		character(len=*),intent(in)::nameleg
		integer,intent(in)::randomNumberScal
		integer::D,rank,ithleg
		type(Dimension)::dimen
		class(Tensor),allocatable::T
		logical::ifDynamic
		allocate(T,mold=inoutT)
		ithleg=inoutT%findOrder(nameleg)
		D=inoutT%dim(ithleg)
		if(D.eq.newD)return
		T=inoutT
		ifDynamic=inoutT%ifDynamic()
		rank=T%Dimension%getRank()
		if(D.gt.newD)then
			call writemess('The input dimension is smaller the one in the Tensor',-1)
			call writemess('ERROR in enlargeTensor',-1)
			call error_stop
		end if
		if(.not.inoutT%getFlag())then
			call writemess('There is no data in the tensor when enlarge D')
			call error_stop
		end if
		if(.not.inoutT%Dimension%if_simple_dimension())then
			call writemess('enlarge D work only on the simple dimension(you can not fuse the Tensor)')
			call error_stop
		end if
		if(rank.eq.1)then
			dimen= [newD] 
		else
			call dimpermute_backwards(dimen,inoutT%Dimension,ithleg)
			dimen=(dimen.subDim.[1,rank-1]) + [newD] 
		end if		
		if(T%Dimension%getNameFlag())then
			call dimen%setName(rank,T%Dimension%getName(ithleg))
		end if
		call inoutT%empty()
		call inoutT%allocate(dimen,T%getType())
		if(ifDynamic)call inoutT%Dynamic()
		call inoutT%random([-randomNumberScal,randomNumberScal])
		call inoutT%setValue(1,T%TData%getTotalData(),T.pb.ithleg)
		inoutT=inoutT.pbi.ithleg
		return
	end subroutine
	subroutine enlargeTensorAllIntRoutine(inoutT,newD,randomNumberScal)
		class(Tensor),intent(inout)::inoutT
		integer,intent(in)::newD(:)
		integer,intent(in)::randomNumberScal
		integer::rank,i
		rank=inoutT%Dimension%getRank()
		call enlargeTensorIntRoutine(inoutT,1,newD(1),randomNumberScal)
		do i=2,rank
			call enlargeTensorIntRoutine(inoutT,i,newD(i),randomNumberScal)
		end do
		return
	end subroutine


	!A*X=B
		!input A  and B
		!output X
		!if A^{-1} may not exit,input RCOND,
		!perform SVD on A, Only keep  singular values S(i) <= RCOND*S(1) 
		!if RCONDM<0 keep the S(i)>0

	function  linequ(A,B,RCOND)
		class(Tensor),allocatable::linequ
		class(Tensor),intent(in)::A
		class(Tensor),intent(in)::B
		class(*),intent(in),optional::RCOND
		integer::Na,Nb,Na2,classtype
		type(dimension)::Xdim
		allocate(linequ,mold=A)
		Na=A.dim.1
		if(A%Dimension%getRank().ne.2)then
			write(*,*)"error in linequ,A should be a matrix"
			call error_stop()
		end if
		Na2=A.dim.2
		if(B%Dimension%getRank().eq.1)then
			Nb=1
		else
			Nb=B.dim.2
		end if
		if(Na.ne.(B.dim.1)) then
			write(*,*)"error in linequ,dimension of A and B"
			call error_stop()
		end if
		classtype=select_type_in_add_minu(A%TData,B%TData)
		if(classtype.eq.1)classtype=2
		Xdim=(A%Dimension.subdim.1)
		if(B%Dimension%getRank().ne.1)then
			Xdim=Xdim+(B%Dimension.subdim.2)
		end if
		call linequ%allocate(Xdim,classtype)
		if(present(RCOND))then
			call linequ2_routine_TData(linequ%TData,A%TData,B%TData,Na,Na2,Nb,RCOND)
		else
			if(Na.ne.Na2) then
				write(*,*)"error in linequ,dimension of A"
				call error_stop()
			end if
			call linequ_routine_TData(linequ%TData,A%TData,B%TData,Na,Nb)
		end if
		return
	end function		

	!! The inverse of a matrix: the input tensor should be a square matrix 
	!A*X=E ==>X=A^-1

	function inverse(T)
		class(Tensor),allocatable ::inverse
		class(Tensor),intent(in) :: T
		class(Tensor),allocatable:: E
		integer :: M,N
		type(dimension)::Xdim
		allocate(inverse,mold=T)
		allocate(E,mold=T)
		if(T%Dimension%getRank().ne.2) then
			write(*,*)"ERROR in calculating the inverse of a Tensor"
			write(*,*)"input Tensor should be a square matrix"
			call error_stop()
		endif
		M = T.dim.1
		N = T.dim.2
		if(M.ne.N) then
			write(*,*)"ERROR in calculating the inverse of a Tensor"
			write(*,*)"input Tensor should be a square matrix"
			call error_stop()
		endif
		call E%setType(T%getclassType())
		E=eye(M,N)
		inverse=linequ(T,E)
		call inverse%resetdim(T%Dimension)
		return
	end function	
	function inverseTen(T,RCOND)result(inverse)
		class(Tensor) ,allocatable::inverse
		class(Tensor),intent(in) :: T
		class(*),intent(in)::RCOND
		class(Tensor),allocatable:: E
		integer :: M,N
		type(dimension)::Xdim
		allocate(inverse,mold=T)
		allocate(E,mold=T)
		if(T%Dimension%getRank().ne.2) then
			write(*,*)"ERROR in calculating the inverse of a Tensor"
			write(*,*)"input Tensor should be a square matrix"
			call error_stop()
		endif
		M = T.dim.1
		N = T.dim.2
		if(M.ne.N) then
			write(*,*)"ERROR in calculating the inverse of a Tensor"
			write(*,*)"input Tensor should be a square matrix"
			call error_stop()
		endif
		call E%setType(T%getclassType())
		E=eye(M,N)
		inverse=linequ(T,E,RCOND)
		call inverse%resetdim(T%Dimension)
		return
	end function	
	
	function inverseTensor(T,RCOND)result(inverse)
		class(Tensor),allocatable ::inverse
		class(Tensor),intent(in) :: T
		class(*),intent(in),optional::RCOND
		class(Tensor),allocatable:: E
		integer :: M,N
		type(dimension)::Xdim
		allocate(inverse,mold=T)
		allocate(E,mold=T)
		if(T%Dimension%getRank().ne.2) then
			write(*,*)"ERROR in calculating the inverse of a Tensor"
			write(*,*)"input Tensor should be a square matrix"
			call error_stop()
		endif
		M = T.dim.1
		N = T.dim.2
		if(M.ne.N) then
			write(*,*)"ERROR in calculating the inverse of a Tensor"
			write(*,*)"input Tensor should be a square matrix"
			call error_stop()
		endif
		call E%setType(T%getclassType())
		E=eye(M,N)
		inverse=linequ(T,E,RCOND)
		call inverse%resetdim(T%Dimension)
		return
	end function

	!***************************************************************************
	!***************************************************************************
	!
	!                                  permutation
	!
	!***************************************************************************
	!**************************************************************************

	function permute_rank3(T1,index_not_permute)result(Res)
		class(Tensor),allocatable::Res
		class(Tensor),intent(in) :: T1
		integer,intent(in) ::   index_not_permute
		integer::dimen(3),lenD
		allocate(Res,mold=T1)
		if(T1%Dimension%getRank().ne.3) then
			write(*,*)"ERROR in permute_rank3"
			write(*,*)"stop"
			call error_stop()
			return
		end if
		call Res%allocate(T1,T1%getType())
		call permutation_data3(Res%TData,T1%TData,index_not_permute,Res%Dimension)
		select case(index_not_permute)
			case(1)
				Res%Dimension=T1%Dimension%Dimpermute((/1,3,2/))
			case(2)
				Res%Dimension=T1%Dimension%Dimpermute((/3,2,1/))
			case(3)
				Res%Dimension=T1%Dimension%Dimpermute((/2,1,3/))
		end select
		return
	 end function
	function permute_rank2 (T)result(Res)
		class(Tensor),allocatable::Res
		class(Tensor),intent(in) :: T
		integer::dimen(2)
		allocate(Res,mold=T)
		if(T%Dimension%getRank().eq.1) then
			Res=T
			return
		end if
		if(T%Dimension%getRank().gt.2) then
			write(*,*)"ERROR in Dpermute_rank2"
			write(*,*)"stop"
			call error_stop()
			return
		end if
		call Res%allocate(T,T%getType())
		call permutation_data2(Res%TData,T%TData,Res%Dimension)
		Res%Dimension=T%Dimension%Dimpermute((/2,1/))
		return
	end function 
	function permutation(T,newOrder)
		class(Tensor),allocatable::permutation
		class(Tensor),intent(in) :: T
		integer,intent(in)::newOrder(:)
		integer,pointer ::inde(:)
		integer::lenOrder,i,j
		allocate(permutation,mold=T)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutation)",-1)
			call error_stop()
		end if
		call WorkingMemory%check('permutation')
		lenorder=size(newOrder)-1
		!allocate(inde(lenorder))
		call WorkingMemory%get_memory(inde,lenorder,'permutation')
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		permutation=T
		do i=lenorder,1,-1
			call permutefo_data_inout(permutation%TData,inde(i),permutation%Dimension)
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end function
	subroutine permutation_routine(T,newOrder)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::newOrder(:)
		integer,pointer ::inde(:)
		integer::lenOrder,i,j
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutation)",-1)
			call error_stop()
		end if
		call WorkingMemory%check('permutation_routine')
		lenorder=size(newOrder)-1
		!allocate(inde(lenorder))
		call WorkingMemory%get_memory(inde,lenorder,'permutation_routine')
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		do i=lenorder,1,-1
			call permutefo_data_inout(T%TData,inde(i),T%Dimension)
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end subroutine
	function permutation_name(T,newOrderchar)result(permutation)
		class(Tensor),allocatable::permutation
		class(Tensor),intent(in) :: T
		CHARACTER(len=*),intent(in)::newOrderchar(:)
		integer,pointer::newOrder(:)
		integer,pointer ::inde(:)
		integer::lenOrder,i,j
		allocate(permutation,mold=T)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutation)",-1)
			call error_stop()
		end if
		call WorkingMemory%check('permutation_name')
		!allocate(newOrder(size(newOrderchar)))
		lenorder=size(newOrderchar)-1
		call WorkingMemory%allocate(1,lenorder+lenorder+2)
		!allocate(inde(lenorder))
		call WorkingMemory%get_memory(newOrder,lenorder+1,'permutation_name,1')
		call WorkingMemory%get_memory(inde,lenorder,'permutation_name,2')
		newOrder=T%Dimension%FindOrder(newOrderchar)
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		permutation=T
		do i=lenorder,1,-1
			call permutefo_data_inout(permutation%TData,inde(i),permutation%Dimension)
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end function
	
	function permutation_Tensor(T,Order)result(permutation)
		class(Tensor),allocatable::permutation
		class(Tensor),intent(in) :: T
		type(Tensor),intent(in)::Order
		integer,pointer ::inde(:)
		integer::lenOrder,i,j
		allocate(permutation,mold=T)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutation)",-1)
			call error_stop()
		end if
		call WorkingMemory%check('permutation_Tensor')
		lenOrder=Order%TData%getTotalData()-1
		select case(Order%getType())
			case (1)
				allocate(inde(lenorder))
				call WorkingMemory%get_memory(inde,lenorder,'permutation_Tensor,1')
				do i=lenorder,1,-1
						inde(i)=Order%ii(i)
				end do
			case (7)
				allocate(inde(lenorder))
				call WorkingMemory%get_memory(inde,lenorder,'permutation_Tensor,2')
				do i=lenorder,1,-1
						inde(i)=T%FindOrder(Order%ai(i))
				end do
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()
		end select
		permutation=T
		do i=lenorder,1,-1
			call permutefo_data_inout(permutation%TData,inde(i),permutation%Dimension)
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end function
	
	subroutine permutation_name_routine(T,newOrderchar)
		class(Tensor),intent(inout) :: T
		CHARACTER(len=*),intent(in)::newOrderchar(:)
		integer,pointer::newOrder(:)
		integer,pointer ::inde(:)
		integer::lenOrder,i,j
		type(Dimension)::dimen		
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutation)",-1)
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			call writemess("split dimension before calling permuation on name",-1)
			call error_stop()
		end if
		call WorkingMemory%check('permutation_name_routine')
		!allocate(newOrder(size(newOrderchar)))
		lenorder=size(newOrderchar)-1
		call WorkingMemory%allocate(1,lenorder+lenorder+2)
		call WorkingMemory%get_memory(newOrder,lenorder+1,'permutation_name_routine,1')
		!allocate(inde(lenorder))
		call WorkingMemory%get_memory(inde,lenorder,'permutation_name_routine,2')
		newOrder=T%Dimension%FindOrder(newOrderchar)
		do i=lenorder,1,-1
				inde(i)=newOrder(i)
		end do
		do i=lenorder,1,-1
			call permutefo_data_inout(T%TData,inde(i),T%Dimension)
			do j=1,i-1
				if(inde(j).lt.inde(i))then
					inde(j)=inde(j)+1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end subroutine
	subroutine permutation_Tensor_routine(T,Tenorder)
		class(Tensor),intent(inout)::T
		type(Tensor),intent(in)::Tenorder
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutefo_vec_name)",-1)
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			call writemess("split dimension before calling permuation on name",-1)
			call error_stop()
		end if
		select case (Tenorder%getType())
			case (1)
				call permutation_routine(T,Tenorder%ii())
			case (7)
				call permutation_name_routine(T,Tenorder%ai())
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()	
		end select
		return
	end subroutine
	
	function permutefo(T,inde)
		class(Tensor),allocatable::permutefo
		class(Tensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		allocate(permutefo,mold=T)
		rank=T%Dimension%getRank()
		if(inde.gt.rank) then
			call writemess("ERROR in function permutefo",-1)
			call writemess(inde+','+rank,-1)
			call error_stop()
		end if
		if(inde.le.0) then
			call writemess("ERROR in function permutefo",-1)
			call writemess("index="+inde,-1)
			call error_stop()
		end if
		if(inde.eq.1) then
			permutefo=T
			return
		end if
		if(inde.eq.rank) then
			permutefo=T.fuse.(/1,rank-1/)
			permutefo=.p.permutefo
			call Dimpermute_forwards(permutefo%Dimension,T%Dimension,inde)
			return
		end if
		permutefo=T.fuse.(/1,inde-1/)
		call permutefo%fuse(3,rank)
		permutefo=permutefo.p.3
		call Dimpermute_forwards(permutefo%Dimension,T%Dimension,inde)
		return
	end function
	subroutine permutefo_routine(T,inde)
		class(Tensor),intent(inout)::T
		integer,intent(in)::inde
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutefo_vec)",-1)
			call error_stop()
		end if
		call permutefo_data_inout(T%TData,inde,T%Dimension)
		return
	end subroutine
	function permutefo_name(T,indechar)result(permutefo)
		class(Tensor),allocatable::permutefo
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::indechar
		integer::inde
		integer::rank
		allocate(permutefo,mold=T)
		rank=T%Dimension%getRank()
		if(if_long_Name(indechar))then
			inde=T%Dimension%FindOrder(indechar)
			if(inde.gt.rank) then
				call writemess("ERROR in function permutefo",-1)
				call writemess(inde+','+rank,-1)
				call error_stop()
			end if
			if(inde.le.0) then
				call writemess("ERROR in function permutefo_name,Can not find the name",-1)
				call writemess("stop",-1)
				call error_stop()
			end if
			if(.not.T%dimension%if_original_dim())then
				call writemess("split dimension before calling permuation on name",-1)
				call error_stop()
			end if
			if(inde.eq.1) then
				permutefo=T
				return
			end if
			if(inde.eq.rank) then
				permutefo=T.fuse.(/1,rank-1/)
				permutefo=.p.permutefo
				call Dimpermute_forwards(permutefo%Dimension,T%Dimension,inde)
				return
			end if
			permutefo=T.fuse.(/1,inde-1/)
			call permutefo%fuse(3,rank)
			permutefo=permutefo.p.3
			call Dimpermute_forwards(permutefo%Dimension,T%Dimension,inde)
			return
		else
			permutefo=T
			call permutefo_name_routine(permutefo,indechar)
			return
		end if
	end function
	subroutine permutefo_name_routine(T,indechar)
		class(Tensor),intent(inout)::T
		integer::inde,i
		character(len=len_of_Name)::indexname
		character(len=*),intent(in)::indechar
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutefo_vec)",-1)
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			call writemess("split dimension before calling permuation on name",-1)
			call error_stop()
		end if
		if(if_long_Name(indechar))then
			inde=T%Dimension%FindOrder(indechar)
			call permutefo_data_inout(T%TData,inde,T%Dimension)
		else
			i=T%Dimension%getRank()
			do inde=T%Dimension%getRank(),1,-1
				indexname=T%Dimension%getName(i).subl.indexsymbol
				if(indechar.equ.indexname)then
					call permutefo_data_inout(T%TData,i,T%Dimension)
				else
					i=i-1
				end if
			end do
		end if
		return
	end subroutine
	
	function permutefo_vec(T,vec_)
		class(Tensor),allocatable::permutefo_vec
		class(Tensor),intent(in)::T
		integer,intent(in)::vec_(:)
		integer::lenVec,i,j
		integer::vec(size(vec_))
		allocate(permutefo_vec,mold=T)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutefo_vec)",-1)
			call error_stop()
		end if
		lenVec=size(vec_)
		vec=vec_
		permutefo_vec=T
		do i=lenVec,1,-1
			call permutefo_data_inout(permutefo_vec%TData,vec(i),permutefo_vec%Dimension)
			do j=1,i-1
				if(vec(j).lt.vec(i))then
					vec(j)=vec(j)+1
				end if
			end do
		end do
		return
	end function
	subroutine permutefo_vec_routine(T,vec_)
		class(Tensor),intent(inout)::T
		integer,intent(in)::vec_(:)
		integer::lenVec,i,j
		integer,pointer::vec(:)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutefo_vec)",-1)
			call error_stop()
		end if
		call WorkingMemory%check('permutefo_vec_routine')
		lenVec=size(vec_)
		call WorkingMemory%get_memory(vec,lenVec,'permutefo_vec_routine')
		vec=vec_
		do i=lenVec,1,-1
			call permutefo_data_inout(T%TData,vec(i),T%Dimension)
			do j=1,i-1
				if(vec(j).lt.vec(i))then
					vec(j)=vec(j)+1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end subroutine
	function permutefo_vec_name(T,indechar)result(permutefo_vec)
		class(Tensor),allocatable::permutefo_vec
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::indechar(:)
		integer::lenVec,i,j
		integer,pointer::vec(:)
		allocate(permutefo_vec,mold=T)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutefo_vec_name)",-1)
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			call writemess("split dimension before calling permuation on name",-1)
			call error_stop()
		end if
		call WorkingMemory%check('permutefo_vec_name')
		lenVec=size(indechar)
		call WorkingMemory%get_memory(vec,lenVec,'permutefo_vec_name')
		vec=T%Dimension%FindOrder(indechar)
		permutefo_vec=T
		do i=lenVec,1,-1
			call permutefo_data_inout(permutefo_vec%TData,vec(i),permutefo_vec%Dimension)
			do j=1,i-1
				if(vec(j).lt.vec(i))then
					vec(j)=vec(j)+1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end function
	
	function permutefo_Tensor(T,Tenorder)result(permutefo)
		class(Tensor),allocatable::permutefo
		class(Tensor),intent(in)::T
		type(Tensor),intent(in)::Tenorder
		character(len=characterlen),pointer::indechar(:)
		integer,pointer::indeint(:)
		allocate(permutefo,mold=T)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutefo_vec_name)",-1)
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			call writemess("split dimension before calling permuation on name",-1)
			call error_stop()
		end if
		call WorkingMemory%check('permutefo_Tensor')
		select case (Tenorder%getType())
			case (1)
				!allocate(indeint(Tenorder%TData%getTotalData()))
				call WorkingMemory%get_memory(indeint,Tenorder%TData%getTotalData(),'permutefo_Tensor,1')
				indeint=Tenorder
				permutefo=permutefo_vec(T,indeint)
			case (7)
				!allocate(indechar(Tenorder%TData%getTotalData()))
				call WorkingMemory%get_memory(indechar,Tenorder%TData%getTotalData(),'permutefo_Tensor,2')
				indechar=Tenorder
				permutefo=permutefo_vec_name(T,indechar)
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()	
		end select
		call WorkingMemory%free()
		return
	end function
	subroutine permutefo_vec_name_routine(T,indechar)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::indechar(:)
		integer::lenVec,i,j
		integer,pointer::vec(:)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutefo_vec_name)",-1)
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			call writemess("split dimension before calling permuation on name",-1)
			call error_stop()
		end if
		call WorkingMemory%check('permutefo_vec_name_routine')
		lenVec=size(indechar)
		call WorkingMemory%get_memory(vec,lenVec,'permutefo_vec_name_routine')
		vec=T%Dimension%FindOrder(indechar)
		do i=lenVec,1,-1
			call permutefo_data_inout(T%TData,vec(i),T%Dimension)
			do j=1,i-1
				if(vec(j).lt.vec(i))then
					vec(j)=vec(j)+1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end subroutine
	
	subroutine permutefo_Tensor_routine(T,Tenorder)
		class(Tensor),intent(inout)::T
		type(Tensor),intent(in)::Tenorder
		character(len=max_len_of_char_in_TData),pointer::indechar(:)
		integer,pointer::indeint(:)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permutefo_vec_name)",-1)
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			call writemess("split dimension before calling permuation on name",-1)
			call error_stop()
		end if
		call WorkingMemory%check('permutefo_Tensor_routine')
		select case (Tenorder%getType())
			case (1)
				!allocate(indeint(Tenorder%TData%getTotalData()))
				call WorkingMemory%get_memory(indeint,Tenorder%TData%getTotalData(),'permutefo_Tensor_routine')
				indeint=Tenorder
				call permutefo_vec_routine(T,indeint)
			case (7)
				!allocate(indechar(Tenorder%TData%getTotalData()))
				call WorkingMemory%get_memory(indechar,Tenorder%TData%getTotalData(),'permutefo_Tensor_routine')
				indechar=Tenorder
				call permutefo_vec_name_routine(T,indechar)
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()	
		end select
		call WorkingMemory%free()
		return
	end subroutine
	function permuteback(T,inde)
		class(Tensor),allocatable::permuteback
		class(Tensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		allocate(permuteback,mold=T)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permuteback)",-1)
			call error_stop()
		end if
		rank=T%Dimension%getRank()
		if(inde.gt.rank) then
			call writemess("ERROR in function permuteback",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		if(inde.le.0) then
			call writemess("ERROR in function permuteback,Can not find the name",-1)
			call writemess("index"+'='+inde,-1)
			call error_stop()
		end if
		if(inde.eq.rank) then
			permuteback=T
			return
		end if
		if(inde.eq.1) then
			permuteback=T.fuse.(/2,rank/)
			permuteback=.p.permuteback
			call Dimpermute_backwards(permuteback%Dimension,T%Dimension,inde)
			return
		end if
		permuteback=T.fuse.(/1,inde-1/)
		call permuteback%fuse(3,rank)
		permuteback=permuteback.p.1
		call Dimpermute_backwards(permuteback%Dimension,T%Dimension,inde)
		return
	end function
	subroutine permuteback_routine(T,inde)
		class(Tensor),intent(inout)::T
		integer,intent(in)::inde
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permuteback_vec)",-1)
			call error_stop()
		end if
		call permuteback_data_inout(T%TData,inde,T%Dimension)
		return
	end subroutine
	subroutine permuteback_name_routine(T,indechar)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::indechar
		integer::inde,i
		character(len=len_of_Name)::indexname
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permuteback_vec)",-1)
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			call writemess("split dimension before calling permuation on name",-1)
			call error_stop()
		end if
		if(if_long_Name(indechar))then
			inde=T%Dimension%FindOrder(indechar)
			call permuteback_data_inout(T%TData,inde,T%Dimension)
		else
			i=1
			do inde=1,T%Dimension%getRank()
				indexname=T%Dimension%getName(i).subl.indexsymbol
				if(indechar.equ.indexname)then
					call permuteback_data_inout(T%TData,i,T%Dimension)
				else
					i=i+1
				end if
			end do
		end if
		return
	end subroutine
	function permuteback_name(T,indechar)result(permuteback)
		class(Tensor),allocatable::permuteback
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::indechar
		integer::inde
		integer::rank
		allocate(permuteback,mold=T)
		if(.not.T%getflag())then
			call writemess("There is no data in the Tensor,(permuteback_name)",-1)
			call error_stop()
		end if
		if(if_long_Name(indechar))then
			if(.not.T%Dimension%if_original_dim())then
				call writemess("split dimension before calling permuation on name",-1)
				call error_stop()
			end if
			rank=T%Dimension%getRank()
			inde=T%dimension%FindOrder(indechar)
			if(inde.gt.rank) then
				call writemess("ERROR in function permuteback_name",-1)
				call writemess("stop",-1)
				call error_stop()
			end if
			if(inde.le.0) then
				call writemess("ERROR in function permuteback_name,Can not find the name",-1)
				call writemess("stop",-1)
				call error_stop()
			end if
			if(inde.eq.rank) then
				permuteback=T
				return
			end if
			if(inde.eq.1) then
				permuteback=T.fuse.(/2,rank/)
				permuteback=.p.permuteback
				call Dimpermute_backwards(permuteback%dimension,T%dimension,inde)
				return
			end if
			permuteback=T.fuse.(/1,inde-1/)
			call permuteback%fuse(3,rank)
			permuteback=permuteback.p.1
			call Dimpermute_backwards(permuteback%dimension,T%dimension,inde)
			return
		else
			permuteback=T
			call permuteback_name_routine(permuteback,indechar)
			return
		end if
	end function
	function permuteback_vec(T,vec_)
		class(Tensor),allocatable::permuteback_vec
		class(Tensor),intent(in)::T
		integer,intent(in)::vec_(:)
		integer::rank,lenVec,i,j
		integer,pointer::vec(:)
		allocate(permuteback_vec,mold=T)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(permuteback_vec)"
			call error_stop()
		end if
		call WorkingMemory%check('permuteback_vec')
		rank=T%Dimension%getRank()
		lenVec=size(vec_)
		call WorkingMemory%get_memory(vec,lenVec,'permuteback_vec')
		vec=vec_
		permuteback_vec=T
		do i=1,lenVec
			call permuteback_data_inout(permuteback_vec%TData,vec(i),permuteback_vec%Dimension)
			do j=lenVec,i+1,-1
				if(vec(j).gt.vec(i))then
					vec(j)=vec(j)-1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end function
	function permuteback_vec_name(T,indechar)result(permuteback_vec)
		class(Tensor),allocatable::permuteback_vec
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::indechar(:)
		integer::rank,lenVec,i,j
		integer,pointer::vec(:)
		allocate(permuteback_vec,mold=T)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(permuteback_vec_name)"
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			write(*,*)"split dimension before calling permuation on name"
			call error_stop()
		end if
		call WorkingMemory%check('permuteback_vec_name')
		rank=T%Dimension%getRank()
		lenVec=size(indechar)
		call WorkingMemory%get_memory(vec,lenVec,'permuteback_vec_name')
		vec=T%Dimension%FindOrder(indechar)
		permuteback_vec=T
		do i=1,lenVec
			call permuteback_data_inout(permuteback_vec%TData,vec(i),permuteback_vec%Dimension)
			do j=lenVec,i+1,-1
				if(vec(j).gt.vec(i))then
					vec(j)=vec(j)-1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end function
	
	function permuteback_Tensor(T,Tenorder)result(permutefo)
		class(Tensor),allocatable::permutefo
		class(Tensor),intent(in)::T
		type(Tensor),intent(in)::Tenorder
		character(len=characterlen),pointer::indechar(:)
		integer,pointer::indeint(:)
		allocate(permutefo,mold=T)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(permutefo_vec_name)"
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			write(*,*)"split dimension before calling permuation on name"
			call error_stop()
		end if
		call WorkingMemory%check('permuteback_Tensor')
		select case (Tenorder%getType())
			case (1)
				!allocate(indeint(Tenorder%TData%getTotalData()))
				call WorkingMemory%get_memory(indeint,Tenorder%TData%getTotalData(),'permuteback_Tensor')
				indeint=Tenorder
				permutefo=permuteback_vec(T,indeint)
			case (7)
				!allocate(indechar(Tenorder%TData%getTotalData()))
				call WorkingMemory%get_memory(indechar,Tenorder%TData%getTotalData(),'permuteback_Tensor')
				indechar=Tenorder
				permutefo=permuteback_vec_name(T,indechar)
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()	
		end select
		call WorkingMemory%free()
		return
	end function
	
	subroutine permuteback_vec_routine(T,vec_)
		class(Tensor),intent(inout)::T
		integer,intent(in)::vec_(:)
		integer::rank,lenVec,i,j
		integer,pointer::vec(:)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(permuteback_vec)"
			call error_stop()
		end if
		call WorkingMemory%check('permuteback_vec_routine')
		rank=T%Dimension%getRank()
		lenVec=size(vec_)
		call WorkingMemory%get_memory(vec,lenVec,'permuteback_vec_routine')
		vec=vec_
		do i=1,lenVec
			call permuteback_data_inout(T%TData,vec(i),T%Dimension)
			do j=lenVec,i+1,-1
				if(vec(j).gt.vec(i))then
					vec(j)=vec(j)-1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end subroutine
	subroutine permuteback_vec_name_routine(T,indechar)
		class(Tensor),intent(inout)::T
		character(len=*),intent(in)::indechar(:)
		integer::rank,lenVec,i,j
		integer,pointer::vec(:)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(permuteback_vec_name)"
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			write(*,*)"split dimension before calling permuation on name"
			call error_stop()
		end if
		call WorkingMemory%check('permuteback_vec_name_routine')
		rank=T%Dimension%getRank()
		lenVec=size(indechar)
		call WorkingMemory%get_memory(vec,lenVec,'permuteback_vec_name_routine')
		vec=T%Dimension%FindOrder(indechar)
		do i=1,lenVec
			call permuteback_data_inout(T%TData,vec(i),T%Dimension)
			do j=lenVec,i+1,-1
				if(vec(j).gt.vec(i))then
					vec(j)=vec(j)-1
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end subroutine
	
	subroutine permuteback_Tensor_routine(T,Tenorder)
		class(Tensor),intent(inout)::T
		type(Tensor),intent(in)::Tenorder
		character(len=characterlen),pointer::indechar(:)
		integer,pointer::indeint(:)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(permutefo_vec_name)"
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			write(*,*)"split dimension before calling permuation on name"
			call error_stop()
		end if
		call WorkingMemory%check('permuteback_Tensor_routine')
		select case (Tenorder%getType())
			case (1)
				!allocate(indeint(Tenorder%TData%getTotalData()))
				call WorkingMemory%get_memory(indeint,Tenorder%TData%getTotalData(),'permuteback_Tensor_routine')
				indeint=Tenorder
				call permuteback_vec_routine(T,indeint)
			case (7)
				!allocate(indechar(Tenorder%TData%getTotalData()))
				call WorkingMemory%get_memory(indechar,Tenorder%TData%getTotalData(),'permuteback_Tensor_routine')
				indechar=Tenorder
				call permuteback_vec_name_routine(T,indechar)
			case default
				call writemess('error in permutation, the data type of order',-1)
				call error_Stop()	
		end select
		call WorkingMemory%free()
		return
	end subroutine

	!****************** permuteInde**************************
	!		T_{1,2,3,..,i,..,n},permuteInde(T,i)=_{2,3,..,i,1,i+1,..,n}

	function permuteInde(T,inde)
		class(Tensor),allocatable::permuteInde
		class(Tensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		allocate(permuteInde,mold=T)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(permuteInde)"
			call error_stop()
		end if
		rank=T%Dimension%getRank()
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permuteInde",inde,rank
			write(*,*)"stop"
			call error_stop()
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permuteInde"
			write(*,*)"index",inde
			call error_stop()
		end if
		if(inde.eq.1) then
			permuteInde=T
			return
		end if
		if(inde.eq.rank) then
			permuteInde=T.fuse.(/2,rank/)
			permuteInde=.p.permuteInde
			call Dimpermute_forwards_index(permuteInde%Dimension,T%Dimension,inde)
			return
		end if
		permuteInde=T.fuse.(/2,inde/)
		call permuteInde%fuse(3,rank)
		permuteInde=permuteInde.p.3
		call Dimpermute_forwards_index(permuteInde%Dimension,T%Dimension,inde)
		return
	end function
	function permuteInde_name(T,w)result(permuteInde)
		class(Tensor),allocatable::permuteInde
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::w
		integer::inde
		integer::rank
		allocate(permuteInde,mold=T)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(permuteInde)"
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			write(*,*)"split dimension before calling permuation on name"
			call error_stop()
		end if
		rank=T%Dimension%getRank()
		inde=T%Dimension%FindOrder(w)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permuteInde",inde,rank
			write(*,*)"stop"
			call error_stop()
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permuteInde"
			write(*,*)"index",inde
			call error_stop()
		end if
		if(inde.eq.1) then
			permuteInde=T
			return
		end if
		if(inde.eq.rank) then
			permuteInde=T.fuse.(/2,rank/)
			permuteInde=.p.permuteInde
			call Dimpermute_forwards_index(permuteInde%Dimension,T%Dimension,inde)
			return
		end if
		permuteInde=T.fuse.(/2,inde/)
		call permuteInde%fuse(3,rank)
		permuteInde=permuteInde.p.3
		call Dimpermute_forwards_index(permuteInde%Dimension,T%Dimension,inde)
		return
	end function

	!****************** permutebackInde**************************
	!		T_{1,2,3,..,i,..,n},permutebackInde(T,i)=_{1,2,3,..,i-1,n,i,i+1,..,n-1}

	function permutebackInde(T,inde)
		class(Tensor),allocatable::permutebackInde
		class(Tensor),intent(in)::T
		integer,intent(in)::inde
		integer::rank
		allocate(Tensor::permutebackInde)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(permutebackInde)"
			call error_stop()
		end if
		rank=T%Dimension%getRank()
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permutebackInde",inde,rank
			write(*,*)"stop"
			call error_stop()
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permutebackInde"
			write(*,*)"index",inde
			call error_stop()
		end if
		if(inde.eq.rank) then
			permutebackInde=T
			return
		end if
		if(inde.eq.1) then
			permutebackInde=T.fuse.(/1,rank-1/)
			permutebackInde=.p.permutebackInde
			call Dimpermute_backwards_index(permutebackInde%Dimension,T%Dimension,inde)
			return
		end if
		permutebackInde=T.fuse.(/1,inde-1/)
		call permutebackInde%fuse(2,rank-inde+1)
		permutebackInde=permutebackInde.p.1
		call Dimpermute_backwards_index(permutebackInde%Dimension,T%Dimension,inde)
		return
	end function
	function permutebackInde_name(T,w)result(permutebackInde)
		class(Tensor),allocatable::permutebackInde
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::w
		integer::inde
		integer::rank
		allocate(permutebackInde,mold=T)
		if(.not.T%getflag())then
			write(*,*)"There is no data in the Tensor,(permutebackInde)"
			call error_stop()
		end if
		if(.not.T%Dimension%if_original_dim())then
			write(*,*)"split dimension before calling permuation on name"
			call error_stop()
		end if
		rank=T%Dimension%getRank()
		inde=T%Dimension%FindOrder(w)
		if(inde.gt.rank) then
			write(*,*)"ERROR in function permutebackInde",inde,rank
			write(*,*)"stop"
			call error_stop()
		end if
		if(inde.le.0) then
			write(*,*)"ERROR in function permutebackInde"
			write(*,*)"index",inde
			call error_stop()
		end if
		if(inde.eq.rank) then
			permutebackInde=T
			return
		end if
		if(inde.eq.1) then
			permutebackInde=T.fuse.(/1,rank-1/)
			permutebackInde=.p.permutebackInde
			call Dimpermute_backwards_index(permutebackInde%Dimension,T%Dimension,inde)
			return
		end if
		permutebackInde=T.fuse.(/1,inde-1/)
		call permutebackInde%fuse(2,rank-inde+1)
		permutebackInde=permutebackInde.p.1
		call Dimpermute_backwards_index(permutebackInde%Dimension,T%Dimension,inde)
		return
	end function



	!***************************************************************************
	!***************************************************************************
	!
	!                                 sub-Tensor
	!
	!***************************************************************************
	!**************************************************************************
	!	inde=[-1,inde_min,inde_max] output data(inde_min:inde_max,:)
	!	or[-2,inde_min,inde_max],data(:,inde_min:inde_max)
	!	or[-3,inde_min,inde_max],data(inde_min:inde_max)
	!	or [-1,inde_row] [-2,inde_col],output row or col
	!	or [inde1_min,inde1_max,inde2_min,inde2_max] output data(inde1_min:inde1_max,inde2_min:inde2_max)
	! output col or row allow for rank>3 , the Tensor will regard as Matrix

	function subTen(T,inde,keepdim)	
		class(Tensor),allocatable::subTen
		class(Tensor),intent(in) ::T
		integer,intent(in)::inde(:)
		logical,optional,intent(in)::keepdim
		integer::dim1,dim2,i,rank,newDim1,newDim2
		type(dimension)::newdim
		allocate(subTen,mold=T)
		rank=T%Dimension%getRank()
		if(size(inde).eq.4) then
			if(T%Dimension%getRank().gt.2)then
				write(*,*)"error in subTen,only matrix or vector is allowed"
				call error_stop()
			end if
			dim1=inde(2)-inde(1)+1
			dim2=inde(4)-inde(3)+1
			call SubTen%allocate((/dim1,dim2/),T%getType())
			call subTen_TData_routine1(T%TData,subTen%TData,T.dim.1,T.dim.2,dim1,&
					dim2,inde(1),inde(2),inde(3),inde(4))
			return
		end if
		if(size(inde).eq.2) then
			dim1=T.dim.1
			dim2=T.dim.2
			select case (inde(1))
			
				case (-1)!output row
					if(T%Dimension%getRank().gt.2)then
						write(*,*)"error in subTen,only matrix or vector is allowed"
						call error_stop()
					end if
					if(present(keepdim).and.keepdim)then
						call SubTen%allocate([1,dim2],T%getType())
					else
						call SubTen%allocate((/dim2/),T%getType())
					end if
					call subTen_TData_routine2(T%TData,subTen%TData,dim1,dim2,dim2,inde(2),.true.)
				case (-2)!
					dim1=T%dim(1)
					newdim=T.subdim.1
					do i=2,rank-1
						newdim=newdim+(T.subdim.i)
						dim1=dim1*T%dim(i)
					end do
					dim2=T%dim(rank)
					if(present(keepdim).and.keepdim)then
						call SubTen%allocate(newdim+[1],T%getType())
						if(newdim%getNameFlag()) call subTen%setName(rank,T%Dimension%getName(rank))
					else
						call SubTen%allocate(newdim,T%getType())
					end if
					call subTen_TData_routine2(T%TData,subTen%TData,dim1,dim2,dim1,inde(2),.false.)
				case default 
					write(*,*) "no such case in subTen"
					write(*,*)inde
					call error_stop()
			end 	select
			return
		end if
		if(size(inde).ne.3) then
			write(*,*) "no such case in subTen"
			write(*,*) "length of inde is ",size(inde)
			write(*,*)inde
			call error_stop()
		end if
		select case (inde(1))
			case (-1)!output some rows
				newDim1=inde(3)-inde(2)+1
				newdim=(/newDim1/)
				newDim2=1
				do i=2,rank
					newdim=newdim+(T.subdim.i)
					newDim2=newDim2*T%dim(i)
				end do
				dim1=T.dim.1
				dim2=T%TData%getTotalData()/dim1
				call SubTen%allocate(newdim,T%getType())
				call subTen_TData_routine3(T%TData,subTen%TData,dim1,dim2,newDim1,&
					newDim2,inde(2),inde(3),.true.)
			case (-2)!!output some cols
				newDim2=inde(3)-inde(2)+1
				newdim=T.subdim.1
				newDim1=T%dim(1)
				do i=2,rank-1
					newdim=newdim+(T.subdim.i)
					newDim1=newDim1*T%dim(i)
				end do
				newdim=newdim+(/newDim2/)
				dim2=T.dim.rank
				dim1=T%TData%getTotalData()/dim2
				call SubTen%allocate(newdim,T%getType())
				call subTen_TData_routine3(T%TData,subTen%TData,dim1,dim2,newDim1,&
					newDim2,inde(2),inde(3),.false.)
			case (-3)!
				dim1=inde(3)-inde(2)+1
				call SubTen%allocate((/dim1/),T%getType())
				call subTen_TData_routine1_(T%TData,subTen%TData,dim1,inde(2),inde(3))
			case default 
				write(*,*) "no such case in subTen"
				write(*,*)inde
				call error_stop()
		end 	select
		return
	end function
	subroutine subTenRoutine(outT,A,legi,legith,keepdim)
		class(Tensor),intent(inout)::outT
		class(Tensor),intent(in)::A
		integer,intent(in)::legi,legith
		logical,optional,intent(in)::keepdim
		integer::dim1,dim2,i,rank,newDim1,newDim2
		type(Tensor)::T
		type(dimension)::newdim
		if(A%Dimension%getRank().le.1)then
			call writemess('ERROR in subTensor, the rank should be larger than or equal to 2')
			call error_stop
		end if
		T=A.pb.legi
		dim1=T%dim(1)
		rank=T%Dimension%getRank()
		newdim=T.subdim.[1,rank-1]
		do i=2,rank-1
			dim1=dim1*T%dim(i)
		end do
		dim2=T%dim(rank)
		call outT%empty()
		newdim=newdim+[1]
		call outT%allocate(newdim,T%getType())
		if(newdim%getNameFlag()) call outT%setName(rank,T%Dimension%getName(rank))
		call subTen_TData_routine2(T%TData,outT%TData,dim1,dim2,dim1,legith,.false.)
		outT=outT.pbi.legi
		if(present(keepdim).and.keepdim)return
		call outT%killLeg(legi,'kill') 
		return
	end subroutine
	function subTen2(T,legi,legith,keepdim)
		class(Tensor),allocatable::subTen2
		class(Tensor),intent(in)::T
		integer,intent(in)::legi,legith
		logical,optional,intent(in)::keepdim
		allocate(Tensor::subTen2)
		call subTenRoutine(subTen2,T,legi,legith,keepdim)
		return
	end function
	function subTen2_Name(T,legName,legith,keepdim)
		class(Tensor),allocatable::subTen2_Name
		class(Tensor),intent(in)::T
		integer,intent(in)::legith
		character(len=*),intent(in)::legName
		logical,optional,intent(in)::keepdim
		allocate(subTen2_Name,mold=T)
		call subTenRoutine(subTen2_Name,T,T%FindOrder(legName),legith,keepdim)
		return
	end function
	subroutine subTenRoutineLegs(outT,A,legi,legith)
		class(Tensor),intent(inout)::outT
		class(Tensor),intent(in)::A
		integer,intent(in)::legi,legith(2)
		integer::dim1,dim2,i,rank,newDim1,newDim2
		class(Tensor),allocatable::T
		type(dimension)::newdim
		if(A%Dimension%getRank().le.1)then
			call writemess('ERROR in subTensor, the rank should be larger than or equal to 2')
			call error_stop
		end if
		allocate(T,mold=A)
		T=A.pb.legi
		rank=T%Dimension%getRank()
		newDim2=legith(2)-legith(1)+1
		if(newDim2.le.0)then
			call writemess('ERROR in subTensor, legith=('+legith(1)+','+legith(2)+')')
			call error_stop
		end if
		if(newDim2.gt.T%dim(rank))then
			call writemess('ERROR in subTensor, legith=('+legith(1)+','+legith(2)+')')
			call writemess('New dimension is larger than original one, dimen='+T%dim(rank))
			call error_stop
		end if
		newdim=T.subdim.[1,rank-1]
		newDim1=T%dim(1)
		do i=2,rank-1
			newDim1=newDim1*T%dim(i)
		end do
		newdim=newdim+[newDim2]
		dim2=T.dim.rank
		dim1=newDim1
		call outT%allocate(newdim,T%getType())
		if(newdim%getNameFlag()) call outT%setName(rank,T%Dimension%getName(rank))
		call subTen_TData_routine3(T%TData,outT%TData,dim1,dim2,newDim1,&
			newDim2,legith(1),legith(2),.false.)
		outT=outT.pbi.legi
		return
	end subroutine
	function subTen3(T,legi,legith)
		class(Tensor),allocatable::subTen3
		class(Tensor),intent(in)::T
		integer,intent(in)::legi,legith(2)
		allocate(subTen3,mold=T)
		call subTenRoutineLegs(subTen3,T,legi,legith)
		return
	end function
	function subTen3_Name(T,legName,legith)
		class(Tensor),allocatable::subTen3_Name
		class(Tensor),intent(in)::T
		integer,intent(in)::legith(2)
		character(len=*),intent(in)::legName
		allocate(subTen3_Name,mold=T)
		call subTenRoutineLegs(subTen3_Name,T,T%FindOrder(legName),legith)
		return
	end function


	!**************   combinationCol   ********************
		!	T1 :a [...,m,n,l] matrix
		!	T2 :a [...,m,n,l] matrix
		!	combination(T1,T2):a [...,m,n,l,2] matrix
		!or 
		!	T1 :a [...,m,n,l] matrix
		!	T2 :a [...,m,n] matrix
		!	combination(T1,T2):a [...,m,n,l+1] matrix

	function combinationCol(T1,T2)result(combination)
		class(Tensor),allocatable::combination
		class(Tensor),intent(in)::T1
		class(Tensor),intent(in)::T2
		integer,allocatable::dim1(:),dim2(:)
		integer::total,total1,i,dim_n,classtype
		type(Dimension)::dimen1,dimen2
		allocate(combination,mold=T1)
		if(T1%Dimension%getRank().eq.T2%Dimension%getRank()) then
			call copydimension(dim1,T1%Dimension)
			call copydimension(dim2,T2%Dimension)
			if(.not.(dim1.equ.dim2)) then
				write(*,*)"can not conbie two Tensor"
				call T1%diminfo('A')
				call T2%diminfo('B')
				write(*,*)"stop"
				call error_stop()
				return
			end if
			classtype=select_type_in_add_minu(T1%TData,T2%TData)
			call combination%allocate(T1%Dimension+[2],classtype)
			call combinationCol_TData(combination%TData,T1%TData,T2%TData)
			return
		end if
		dimen1=T1%Dimension
		call dimen1%fuse(1,T1%Dimension%getRank()-1)
		dimen2=T2%Dimension
		call dimen2%fuse(1,T2%Dimension%getRank())
		if(.not.((dimen1.subdim.1).equ.dimen2)) then
			write(*,*)"can not conbie two Tensor"
			write(*,*)"stop2"
			call error_stop()
			return
		end if
		dim_n=dimen1.dim.2
		classtype=select_type_in_add_minu(T1%TData,T2%TData)
		call combination%allocate(T2%Dimension+(/dim_n+1/),classtype)
		call combinationCol_TData(combination%TData,T1%TData,T2%TData)
		return
	end function
	!**************   combinationrow   ********************
		!	T1 :a [l,m,n,...] matrix
		!	T2 :a [l,m,n,...] matrix
		!	combinationrow(T1,T2):a [2,l,m,n,...] matrix
		!	or 
		!	T1 :a [l,m,n,...] matrix
		!	T2 :a [m,n,...] matrix
		!	combinationrow(T1,T2):a [l+1,m,n,...] matrix	

	function combinationrow(T1,T2)
		class(Tensor),allocatable::combinationrow
		class(Tensor),intent(in)::T1
		class(Tensor),intent(in)::T2
		integer,allocatable::dim1(:),dim2(:)
		integer::total1,total2,i,dim_n,classtype
		type(Dimension)::dimen1,dimen2
		allocate(combinationrow,mold=T1)
		if(T1%Dimension%getRank().eq.T2%Dimension%getRank()) then
			call copydimension(dim1,T1%Dimension)
			call copydimension(dim2,T2%Dimension)
			if(.not.(dim1.equ.dim2)) then
				write(*,*)"can not combinationrow two Tensor"
				write(*,*)dim1
				write(*,*)dim2
				write(*,*)"stop"
				call error_stop()
				return
			end if
			total1=T1%TData%getTotalData()
			total2=T2%TData%getTotalData()
			classtype=select_type_in_add_minu(T1%TData,T2%TData)
			call combinationrow%allocate([2]+T1%Dimension,classtype)
			call combinationRow_TData(combinationrow%TData,T1%TData,T2%TData,2,total1,1,1)
			return
		end if
		dimen1=T1%Dimension
		call dimen1%fuse(2,T1%Dimension%getRank())
		dimen2=T2%Dimension
		call dimen2%fuse(1,T2%Dimension%getRank())
		if(.not.((dimen1.subdim.2).equ.dimen2)) then
			write(*,*)"can not combinationrow two Tensor"
			write(*,*)"stop2"
			call error_stop()
			return
		end if
		dim_n=dimen1.dim.1
		total1=T1%TData%getTotalData()
		total2=T2%TData%getTotalData()
		classtype=select_type_in_add_minu(T1%TData,T2%TData)
		call combinationrow%allocate([dim_n+1]+T2%Dimension,classtype)
		call combinationRow_TData(combinationrow%TData,T1%TData,T2%TData,dim_n+1,total2,dim_n,1)
		return
	end function		
	
	!pasteTensor:
		!pasteTensor(T1,T2,.true.)
		!			T1 :a [l,m,n,...] matrix
		!			T2 :a [l,m,n,...] matrix
		!			pasteTensor(T1,T2):a [2*l,m,n,...] matrix
		!        [1--->l, m, n,...] is T1
		!			[l+1--->2*l, m, n,...] is T2
		!        /   \
		!        | T1 |
		!        |----|
		!        | T2 |
		!         \    /
		!pasteTensor(T1,T2,.false.)
		!			T1 :a [...,m,n,l] matrix
		!			T2 :a [...,m,n,l] matrix
		!			pasteTensor(T1,T2):a [...,m,n,2*l] matrix
		!        [... , m, n, 1--->l] is T1
		!			[... , m, n, l+1--->2*l] is T2
		!        /         \
		!        | T1 | T2 |
		!        \         /
		!        

	subroutine pasteTensorSubroutine(T1,T2,row)
		class(Tensor),intent(inout)::T1
		class(Tensor),intent(in)::T2
		logical,intent(in)::row
		integer::i,classtype,rank1,rank2,pasteDim1,pasteDim2,collen,newDim_i
		type(Dimension)::newDim
		class(Tensor),allocatable::pasteTensor
		allocate(pasteTensor,mold=T1)
		if(.not.T2%getFlag())return
		rank1=T1%Dimension%getRank()
		rank2=T2%Dimension%getRank()
		if(rank1.ne.rank2) then
			call writemess('can not paste two Tensor,ranks are,'+rank1+','+rank2,-1)
			call error_stop()
		end if
		if(row)then
			if(rank1.eq.1)then
				!call writemess('Do not finsiehd this part,can not paste two Tensor,ranks are,'+rank1+','+rank2)
				!call error_stop()
				pasteDim1=1
				pasteDim2=1
				newDim_i=2
				collen=T1%TData%getTotalData()
				newDim=(/2,collen/)
			else
				pasteDim1=T1%dim(1)
				pasteDim2=T2%dim(1)
				newDim_i=pasteDim1+pasteDim2
				newDim=(/newDim_i/)
				collen=1
				do i=2,rank1
					newDim=newDim+(T1.subDim.i)
					collen=collen*T1%dim(i)
					if(T1%dim(i).ne.T2%dim(i))then
						call writemess('can not paste two Tensor,in the ,'+i+'th dimension',-1)
						call writemess('T1%dim('+i+')='+T1%dim(i),-1)
						call writemess('T2%dim('+i+')='+T2%dim(i),-1)
						call error_stop()
					end if
				end do
			end if
			if(T1%Dimension%getNameFlag()) then
				call newDim%setName(1,T1%Dimension%getName(1))
			end if
			classtype=select_type_in_add_minu(T1%TData,T2%TData)
			call pasteTensor%allocate(newDim,classtype)
			call combinationRow_TData(pasteTensor%TData,T1%TData,T2%TData,newDim_i,collen,pasteDim1,pasteDim2)
		else
			pasteDim1=T1%dim(rank1)
			pasteDim2=T2%dim(rank2)
			newDim_i=pasteDim1+pasteDim2
			collen=1
			do i=1,rank1-1
				newDim=newDim+(T1.subDim.i)
				collen=collen*T1%dim(i)
				if(T1%dim(i).ne.T2%dim(i))then
					call writemess('can not paste two Tensor,in the ,'+i+'th dimension',-1)
					call writemess('T1%dim('+i+')='+T1%dim(i),-1)
					call writemess('T2%dim('+i+')='+T2%dim(i),-1)
					call error_stop()
				end if
			end do
			newDim=newDim+(/newDim_i/)
			if(T1%Dimension%getNameFlag()) then
				call newDim%setName(rank1,T1%Dimension%getName(rank1))
			end if
			classtype=select_type_in_add_minu(T1%TData,T2%TData)
			call pasteTensor%allocate(newDim,classtype)
			call combinationCol_TData(pasteTensor%TData,T1%TData,T2%TData)
		end if
		T1=pasteTensor
		return
	end subroutine
	
	function pasteTensorRow(T1,T2)result(pasteTensor)
		class(Tensor),allocatable::pasteTensor
		class(Tensor),intent(in)::T1
		class(Tensor),intent(in)::T2
		integer::i,classtype,rank1,rank2,pasteDim1,pasteDim2,collen,newDim_i
		type(Dimension)::newDim
		allocate(pasteTensor,mold=T1)
		rank1=T1%Dimension%getRank()
		rank2=T2%Dimension%getRank()
		if(rank1.ne.rank2) then
			call writemess('can not paste two Tensor,ranks are,'+rank1+','+rank2,-1)
			call error_stop()
		end if
		if(rank1.eq.1)then
			call writemess('Do not finsiehd this part,can not paste two Tensor,ranks are,'+rank1+','+rank2,-1)
			call error_stop()
		end if
		pasteDim1=T1%dim(1)
		pasteDim2=T2%dim(1)
		newDim_i=pasteDim1+pasteDim2
		newDim=(/newDim_i/)
		collen=1
		do i=2,rank1
			newDim=newDim+(T1.subDim.i)
			collen=collen*T1%dim(i)
			if(T1%dim(i).ne.T2%dim(i))then
				call writemess('can not paste two Tensor,in the ,'+i+'th dimension',-1)
				call writemess('T1%dim('+i+')='+T1%dim(i))
				call writemess('T2%dim('+i+')='+T2%dim(i))
				call error_stop()
			end if
		end do
		if(T1%Dimension%getNameFlag()) then
			call newDim%setName(1,T1%Dimension%getName(1))
		end if
		classtype=select_type_in_add_minu(T1%TData,T2%TData)
		call pasteTensor%allocate(newDim,classtype)
		call combinationRow_TData(pasteTensor%TData,T1%TData,T2%TData,newDim_i,collen,pasteDim1,pasteDim2)
		return
	end function
	
	function pasteTensorCol(T1,T2)result(pasteTensor)
		class(Tensor),allocatable::pasteTensor
		class(Tensor),intent(in)::T1
		class(Tensor),intent(in)::T2
		integer::i,classtype,rank1,rank2,pasteDim1,pasteDim2,collen,newDim_i
		type(Dimension)::newDim
		allocate(pasteTensor,mold=T1)
		rank1=T1%Dimension%getRank()
		rank2=T2%Dimension%getRank()
		if(rank1.ne.rank2) then
			call writemess('can not paste two Tensor,ranks are,'+rank1+','+rank2,-1)
			call error_stop()
		end if
		pasteDim1=T1%dim(rank1)
		pasteDim2=T2%dim(rank2)
		newDim_i=pasteDim1+pasteDim2
		collen=1
		do i=1,rank1-1
			newDim=newDim+(T1.subDim.i)
			collen=collen*T1%dim(i)
			if(T1%dim(i).ne.T2%dim(i))then
				call writemess('can not paste two Tensor,in the ,'+i+'th dimension',-1)
				call writemess('T1%dim('+i+')='+T1%dim(i),-1)
				call writemess('T2%dim('+i+')='+T2%dim(i),-1)
				call error_stop()
			end if
		end do
		newDim=newDim+(/newDim_i/)
		if(T1%Dimension%getNameFlag()) then
			call newDim%setName(rank1,T1%Dimension%getName(rank1))
		end if
		classtype=select_type_in_add_minu(T1%TData,T2%TData)
		call pasteTensor%allocate(newDim,classtype)
		call combinationCol_TData(pasteTensor%TData,T1%TData,T2%TData)
		return
	end function
	

	!***********************************************************************
	!***********************************************************************
	!
	!                    SVD
	!
	!***********************************************************************
	!***********************************************************************	

	subroutine SVDcutoff(T,U,s,V,Ncut_) 
		class(Tensor),target,intent(in)::T
		class(Tensor),target,intent(inout)::U
		class(Tensor),target,intent(inout)::s
		class(Tensor),target,intent(inout)::V
		integer,optional,intent(in)::Ncut_
		integer::Ncut,info
		integer ::min_MN,m,n,classtype
		type(Dimension) :: T1Dim,T2Dim
		class(Tensor),pointer::Tp
		class(Tensor),pointer::Up,sp,Vp
		if(T%Dimension%getRank().ne.2) then
			if(.not.T%getflag()) then
				write(*,*)"ERROR in svd"
				write(*,*)"There is no data in the Tensor"
				write(*,*)"stop"
				call error_stop()
			end if
			write(*,*) "Input Tensor should be 2 dimension in svd"
			write(*,*)"stop"
			call error_stop()
		end if
		Tp=>T
		Up=>U
		sp=>s
		Vp=>V
		if(associated(Tp,Up).or.associated(Tp,sp).or.associated(Tp,Vp).or.associated(Up,sp).or.&
				associated(Up,Vp).or.associated(sp,Vp))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%SVDroutine(U,S,V,Ncut)',-1)
			call writemess('T, U, s and V can not be a same variable',-1)
			call error_stop
		end if
		Tp=>null()
		Up=>null()
		sp=>null()
		Vp=>null()
		
		m=T.dim.1
		n=T.dim.2
		min_MN=min(M,N)
		classtype=T%getType()
		call U%empty()
		call s%empty()
		call V%empty()
		if(present(Ncut_))then
			Ncut=Ncut_
			if(Ncut.gt.min_MN) Ncut=min_MN
			T1Dim=(/Ncut/)
			T1Dim=(T%dimension.subdim.1)+T1Dim
			T2Dim=(/Ncut/)
			T2Dim=T2Dim+(T%dimension.subdim.2)
			if(classtype.eq.1)then
				call U%allocate(T1Dim,2)
				call V%allocate(T2Dim,2)
			else
				call U%allocate(T1Dim,classtype)
				call V%allocate(T2Dim,classtype)
			end if
			select case(classtype)
				case(5)
					call s%allocate((/Ncut/),3)
				case(3)
					call s%allocate((/Ncut/),3)
				case default
					call s%allocate((/Ncut/),2)
			end select
		else
			T1Dim=(/min_MN/)
			T1Dim=(T%dimension.subdim.1)+T1Dim
			T2Dim=(/min_MN/)
			T2Dim=T2Dim+(T%dimension.subdim.2)
			if(classtype.eq.1)then
				call U%allocate(T1Dim,2)
				call V%allocate(T2Dim,2)
			else
				call U%allocate(T1Dim,classtype)
				call V%allocate(T2Dim,classtype)
			end if
			select case(classtype)
				case(5)
					call s%allocate((/min_MN/),3)
				case(3)
					call s%allocate((/min_MN/),3)
				case default
					call s%allocate((/min_MN/),2)
			end select
		end if
		call SVD_TData_routine(T%TData,U%TData,S%TData,V%TData,m,n,min_MN,Ncut_,info)
		if(info.ne.0) then
			call writemess('Error in svd ,info='+info,-1)
			call writemess('output The data in ./_SVD_ERROR_LOG.err',-1)
			open(unit=9991,file='./_SVD_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
			call T%writeinfo('The Matrix in SVD',9991)
			call U%writeinfo('The Matrix U in SVD',9991)
			call S%writeinfo('The Matrix S in SVD',9991)
			call V%writeinfo('The Matrix V in SVD',9991)
			close(9991)
			call error_stop()
		end if
		if(S%isnan())then!The number in S is less than U,V and T
			call writemess('Error in svd ,NAN ERROR',-1)
			call writemess('output The data in ./_SVD_ERROR_LOG.err',-1)
			open(unit=9991,file='./_SVD_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
			call T%writeinfo('The Matrix in SVD',9991)
			call U%writeinfo('The Matrix U in SVD',9991)
			call S%writeinfo('The Matrix S in SVD',9991)
			call V%writeinfo('The Matrix V in SVD',9991)
			close(9991)
			call error_stop()
		end if
		return
	end subroutine
	subroutine SVDcutoff_kill_inData(T,U,s,V,Ncut_) 
		class(Tensor),target,intent(inout)::T
		class(Tensor),target,intent(inout)::U
		class(Tensor),target,intent(inout)::s
		class(Tensor),target,intent(inout)::V
		integer,optional,intent(in)::Ncut_
		integer::Ncut,info
		integer ::min_MN,m,n,classtype
		type(Dimension) :: T1Dim,T2Dim
		class(Tensor),pointer::Tp
		type(Tensor),pointer::Up,sp,Vp
		if(deallocate_memory_flag)then
			call writemess(' The subroutine of SVDcutoff_kill_inData Can not use when deallocate_memory_flag=.true.',-1)
			call writemess(' One should use "call unset_deallocate_memory_flag()" and than use this subroutine',-1)
			call error_stop
		end if
		if(T%Dimension%getRank().ne.2) then
			if(.not.T%getflag()) then
				write(*,*)"ERROR in svd"
				write(*,*)"There is no data in the Tensor"
				write(*,*)"stop"
				call error_stop()
			end if
			write(*,*) "Input Tensor should be 2 dimension in svd"
			write(*,*)"stop"
			call error_stop()
		end if
		Tp=>T
		Up=>U
		sp=>s
		Vp=>V
		if(associated(Tp,Up).or.associated(Tp,sp).or.associated(Tp,Vp).or.associated(Up,sp).or.&
				associated(Up,Vp).or.associated(sp,Vp))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%SVDroutine(U,S,V,Ncut)',-1)
			call writemess('T, U, s and V can not be a same variable',-1)
			call error_stop
		end if
		Tp=>null()
		Up=>null()
		sp=>null()
		Vp=>null()
		
		call U%empty()
		call s%empty()
		call V%empty()
		m=T.dim.1
		n=T.dim.2
		min_MN=min(M,N)

		classtype=T%getType()
		T1Dim=(/min_MN/)
		T1Dim=(T%Dimension.subdim.1)+T1Dim
		T2Dim=(/min_MN/)
		T2Dim=T2Dim+(T%Dimension.subdim.2)
		if(classtype.eq.1)then
			call U%allocate(T1Dim,2)
			call V%allocate(T2Dim,2)
		else
			call U%allocate(T1Dim,classtype)
			call V%allocate(T2Dim,classtype)
		end if
		select case(classtype)
			case(5)
				call s%allocate((/min_MN/),3)
			case(3)
				call s%allocate((/min_MN/),3)
			case default
				call s%allocate((/min_MN/),2)
		end select

		!if deallocate_memory_flag=.false. call U%empty() will not deallocate the momery of U
		! when trucate the data, just reset the info in dimension
		if(present(Ncut_))then
			Ncut=Ncut_
			if(Ncut.gt.min_MN) Ncut=min_MN
			T1Dim=(/Ncut/)
			T1Dim=(T%Dimension.subdim.1)+T1Dim
			T2Dim=(/Ncut/)
			T2Dim=T2Dim+(T%Dimension.subdim.2)
			call U%reset_dim_no_check(T1Dim)
			call V%reset_dim_no_check(T2Dim)
			call S%reset_dim_no_check([Ncut])
			!call U%empty()
			!call S%empty()
			!call V%empty()
			!if(classtype.eq.1)then
			!	call U%allocate(T1Dim,2)
			!	call V%allocate(T2Dim,2)
			!else
			!	call U%allocate(T1Dim,classtype)
			!	call V%allocate(T2Dim,classtype)
			!end if
			!select case(classtype)
			!	case(5)
			!		call S%allocate([Ncut],3)
			!	case(3)
			!		call S%allocate([Ncut],3)
			!	case default
			!		call S%allocate([Ncut],2)
			!end select
		end if
		call SVD_TData_routine_Kill_inData(T%TData,U%TData,S%TData,V%TData,m,n,min_MN,Ncut_,info)
		if(info.ne.0) then
			call writemess('Error in svd ,info='+info,-1)
			call writemess('output The data in ./_SVD_ERROR_LOG.err',-1)
			open(unit=9991,file='./_SVD_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
			call T%writeinfo('The Matrix in SVD',9991)
			call U%writeinfo('The Matrix U in SVD',9991)
			call S%writeinfo('The Matrix S in SVD',9991)
			call V%writeinfo('The Matrix V in SVD',9991)
			close(9991)
			call error_stop()
		end if
		if(S%isnan())then!The number in S is less than U,V and T
			call writemess('Error in svd ,NAN ERROR',-1)
			call writemess('output The data in ./_SVD_ERROR_LOG.err',-1)
			open(unit=9991,file='./_SVD_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
			call T%writeinfo('The Matrix in SVD',9991)
			call U%writeinfo('The Matrix U in SVD',9991)
			call S%writeinfo('The Matrix S in SVD',9991)
			call V%writeinfo('The Matrix V in SVD',9991)
			close(9991)
			call error_stop()
		end if
		return
	end subroutine
	subroutine SVDcutoff_name(inputT,U,s,V,nameU,nameV,Ncut) 
		class(Tensor),target,intent(in)::inputT
		class(Tensor),target,intent(inout)::U
		class(Tensor),target,intent(inout)::s
		class(Tensor),target,intent(inout)::V
		class(Tensor),allocatable::T
		character(len=*),intent(in)::nameU,nameV
		integer,optional,intent(in)::Ncut
		integer::rankU,rankV,rank,i
		if(.not.inputT%getflag()) then
			write(*,*)"ERROR in svd"
			write(*,*)"There is no data in the Tensor"
			write(*,*)"stop"
			call error_stop()
		end if
		allocate(T,mold=inputT)
		T=inputT
		rank=T%Dimension%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if((T%Dimension%getName(i).subl.indexsymbol).equ.nameU) rankU=rankU+1
			if((T%Dimension%getName(i).subl.indexsymbol).equ.nameV) rankV=rankV+1
		end do
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in SVDcutoff_name",-1)
			call writemess(nameU+','+nameV,-1)
			call inputT%diminfo()
			call writemess(rankU+','+rankV+','+rank,-1)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameU,-1)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameV,-1)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			call T%forward(nameU)
		else
			call T%backward(nameV)
		end if
		call T%fuse(1,rankU)
		call T%fuse(2,rankV+1)
		call SVDcutoff(T,U,s,V,Ncut)
		call U%split()
		call V%split()
		call U%setName(U%Dimension%getRank(),'SVD.U')
		call V%setName(1,'SVD.V')
		return
	end subroutine
	subroutine SVDcutoff_name_kill_inData(T,U,s,V,nameU,nameV,Ncut) 
		class(Tensor),target,intent(inout)::T
		class(Tensor),target,intent(inout)::U
		class(Tensor),target,intent(inout)::s
		class(Tensor),target,intent(inout)::V
		character(len=*),intent(in)::nameU,nameV
		integer,optional,intent(in)::Ncut
		integer::rankU,rankV,rank,i
		if(.not.T%getflag()) then
			write(*,*)"ERROR in svd"
			write(*,*)"There is no data in the Tensor"
			write(*,*)"stop"
			call error_stop()
		end if
		
		rank=T%Dimension%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if((T%Dimension%getName(i).subl.indexsymbol).equ.nameU) rankU=rankU+1
			if((T%Dimension%getName(i).subl.indexsymbol).equ.nameV) rankV=rankV+1
		end do
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in SVDcutoff_name",-1)
			call writemess(nameU+','+nameV,-1)
			call T%diminfo()
			call writemess(rankU+','+rankV+','+rank,-1)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameU,-1)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in SVDcutoff_name,no such name",-1)
			call writemess(nameV,-1)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			call T%forward(nameU)
		else
			call T%backward(nameV)
		end if
		call T%fuse(1,rankU)
		call T%fuse(2,rankV+1)
		call SVDcutoff_kill_inData(T,U,s,V,Ncut)
		call U%split()
		call V%split()
		call U%setName(U%Dimension%getRank(),'SVD.U')
		call V%setName(1,'SVD.V')
		return
	end subroutine
	
	
	function SVDTensor_noname(T,Ncut)	result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),intent(in)::T
		integer,optional,intent(in)::Ncut
		allocate(res(3),mold=T)
		call SVDcutoff(T,res(1),res(2),res(3),Ncut) 
		return
	end function
	
	function SVDTensor_name(T,nameU,nameV,Ncut)	result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),intent(in)::T
		integer,optional,intent(in)::Ncut
		character(len=*),intent(in)::nameU,nameV
		allocate(res(3),mold=T)
		call SVDcutoff_name(T,res(1),res(2),res(3),nameU,nameV,Ncut) 
		return
	end function
	
	function SVDNameLeft(T,LegName,Ncut,Left)result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),intent(in)::T
		integer,intent(in)::Ncut
		logical,intent(in)::Left
		character(len=*),intent(in)::LegName(:)
		class(Tensor),allocatable::Temp
		integer::lenName,rank
		allocate(res(3),mold=T)
		allocate(Temp,mold=T)
		lenName=size(LegName)
		rank=T%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			Temp=T.pf.LegName
			call Temp%fuse(1,lenName)
			call Temp%fuse(2,rank)
			call SVDcutoff(Temp,res(1),res(2),res(3),Ncut)
		else 
			Temp=T.pb.LegName
			call Temp%fuse(1,rank-lenName)
			call Temp%fuse(2,Temp%Dimension%getRank())
			call SVDcutoff(Temp,res(1),res(2),res(3),Ncut)
		end if
		call res(1)%split()
		call res(3)%split()

		call res(1)%setName(res(1)%Dimension%getRank(),'SVD.U')
		call res(3)%setName(1,'SVD.V')
		return
	end function
	subroutine SVDNameLeftRoutine(T,U,s,V,LegName,Ncut,Left)
		class(Tensor),intent(in)::T
		class(Tensor),intent(inout)::U
		class(Tensor),intent(inout)::s
		class(Tensor),intent(inout)::V
		logical,intent(in)::Left
		integer,intent(in)::Ncut
		character(len=*),intent(in)::LegName(:)
		class(Tensor),allocatable::Temp
		integer::lenName,rank
		allocate(Temp,mold=T)
		lenName=size(LegName)
		rank=T%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			Temp=T.pf.LegName
			call Temp%fuse(1,lenName)
			call Temp%fuse(2,rank)
			call SVDcutoff(Temp,U,s,V,Ncut)
		else 
			Temp=T.pb.LegName
			call Temp%fuse(1,rank-lenName)
			call Temp%fuse(2,Temp%Dimension%getRank())
			call SVDcutoff(Temp,U,s,V,Ncut)
		end if
		call U%split()
		call V%split()
		call U%setName(U%Dimension%getRank(),'SVD.U')
		call V%setName(1,'SVD.V')
		return
	end subroutine
	subroutine SVDNameLeftRoutine_kill_indata(T,U,s,V,LegName,Ncut,Left)
		class(Tensor),intent(inout)::T
		type(Tensor),intent(inout)::U,s,V
		integer,intent(in)::Ncut
		logical,intent(in)::Left
		character(len=*),intent(in)::LegName(:)
		integer::lenName,rank
		lenName=size(LegName)
		rank=T%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			call T%forward(LegName)
			call T%fuse(1,lenName)
			call T%fuse(2,rank)
			call SVDcutoff_kill_inData(T,U,s,V,Ncut)
		else
			call T%Backward(LegName)
			call T%fuse(1,rank-lenName)
			call T%fuse(2,T%Dimension%getRank())
			call SVDcutoff_kill_inData(T,U,s,V,Ncut)
		end if
		call U%split()
		call V%split()
		call U%setName(U%Dimension%getRank(),'SVD.U')
		call V%setName(1,'SVD.V')
		return
	end subroutine

	function SVDNameLeft_no_cut(T,LegName,Left)result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),intent(in)::T
		logical,intent(in)::Left
		character(len=*),intent(in)::LegName(:)
		class(Tensor),allocatable::Temp
		integer::lenName,rank
		allocate(res(3),mold=T)
		allocate(Temp,mold=T)
		lenName=size(LegName)
		rank=T%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			Temp=T.pf.LegName
			call Temp%fuse(1,lenName)
			call Temp%fuse(2,rank)
			call SVDcutoff(Temp,res(1),res(2),res(3))
		else 
			Temp=T.pb.LegName
			call Temp%fuse(1,rank-lenName)
			call Temp%fuse(2,Temp%Dimension%getRank())
			call SVDcutoff(Temp,res(1),res(2),res(3))
		end if
		call res(1)%split()
		call res(3)%split()

		call res(1)%setName(res(1)%Dimension%getRank(),'SVD.U')
		call res(3)%setName(1,'SVD.V')
		return
	end function
	subroutine SVDNameLeftRoutine_no_cut(T,U,s,V,LegName,Left)
		class(Tensor),intent(in)::T
		class(Tensor),intent(inout)::U
		class(Tensor),intent(inout)::s
		class(Tensor),intent(inout)::V
		logical,intent(in)::Left
		character(len=*),intent(in)::LegName(:)
		class(Tensor),allocatable::Temp
		integer::lenName,rank
		allocate(Temp,mold=T)
		lenName=size(LegName)
		rank=T%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			Temp=T.pf.LegName
			call Temp%fuse(1,lenName)
			call Temp%fuse(2,rank)
			call SVDcutoff(Temp,U,s,V)
		else 
			Temp=T.pb.LegName
			call Temp%fuse(1,rank-lenName)
			call Temp%fuse(2,Temp%Dimension%getRank())
			call SVDcutoff(Temp,U,s,V)
		end if
		call U%split()
		call V%split()
		call U%setName(U%Dimension%getRank(),'SVD.U')
		call V%setName(1,'SVD.V')
		return
	end subroutine
	subroutine SVDNameLeftRoutine_kill_indata_no_cut(T,U,s,V,LegName,Left)
		class(Tensor),intent(inout)::T
		type(Tensor),intent(inout)::U,s,V
		logical,intent(in)::Left
		character(len=*),intent(in)::LegName(:)
		integer::lenName,rank
		lenName=size(LegName)
		rank=T%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			call T%forward(LegName)
			call T%fuse(1,lenName)
			call T%fuse(2,rank)
			call SVDcutoff_kill_inData(T,U,s,V)
		else
			call T%Backward(LegName)
			call T%fuse(1,rank-lenName)
			call T%fuse(2,T%Dimension%getRank())
			call SVDcutoff_kill_inData(T,U,s,V)
		end if
		call U%split()
		call V%split()
		call U%setName(U%Dimension%getRank(),'SVD.U')
		call V%setName(1,'SVD.V')
		return
	end subroutine

	function SVDTensor_Leg2(T,LegNameRow,LegNameCol,Ncut)result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),intent(in)::T
		integer,optional,intent(in)::Ncut
		character(len=*),intent(in)::LegNameRow(:),LegNameCol(:)
		class(Tensor),allocatable::Temp
		integer::lenName1,lenName2,rank
		allocate(res(3),mold=T)
		allocate(Temp,mold=T)
		lenName1=size(LegNameRow)
		lenName2=size(LegNameCol)
		rank=T%Dimension%getRank()
		if((lenName1+lenName2).ne.rank)then
			call writemess('size(LegNameRow)='+lenName1+'size(LegNameCol)='+lenName2,-1)
			call writemess('rank='+rank,-1)
			call writemess('rank should be equal to size(LegNameRow)+size(LegNameCol)')
			call error_stop
		end if
		if(lenName1.le.lenName2)then
			Temp=T.pf.LegNameRow
		else
			Temp=T.pb.LegNameCol
		end if
		call Temp%fuse(1,lenName1)
		call Temp%fuse(2,lenName2+1)
		call SVDcutoff(Temp,res(1),res(2),res(3),Ncut)
		call res(1)%split()
		call res(3)%split()

		call res(1)%setName(res(1)%Dimension%getRank(),'SVD.U')
		call res(3)%setName(1,'SVD.V')
		return
	end function
	subroutine SVDTensor_Leg2_Routine(T,U,s,V,LegNameRow,LegNameCol,Ncut)
		class(Tensor),intent(in)::T
		class(Tensor),intent(inout)::U
		class(Tensor),intent(inout)::s
		class(Tensor),intent(inout)::V
		integer,optional,intent(in)::Ncut
		character(len=*),intent(in)::LegNameRow(:),LegNameCol(:)
		class(Tensor),allocatable::Temp
		integer::lenName1,lenName2,rank
		allocate(Temp,mold=T)
		lenName1=size(LegNameRow)
		lenName2=size(LegNameCol)
		rank=T%Dimension%getRank()
		if((lenName1+lenName2).ne.rank)then
			call writemess('size(LegNameRow)='+lenName1+'size(LegNameCol)='+lenName2,-1)
			call writemess('rank='+rank,-1)
			call writemess('rank should be equal to size(LegNameRow)+size(LegNameCol)')
			call error_stop
		end if
		if(lenName1.le.lenName2)then
			Temp=T.pf.LegNameRow
		else
			Temp=T.pb.LegNameCol
		end if
		call Temp%fuse(1,lenName1)
		call Temp%fuse(2,lenName2+1)
		call SVDcutoff(Temp,U,s,V,Ncut)
		call U%split()
		call V%split()
		call U%setName(U%Dimension%getRank(),'SVD.U')
		call V%setName(1,'SVD.V')
		return
	end subroutine
	subroutine SVDTensor_Leg2_Routine_kill_inData(T,U,s,V,LegNameRow,LegNameCol,Ncut)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(inout)::U
		class(Tensor),intent(inout)::s
		class(Tensor),intent(inout)::V
		integer,optional,intent(in)::Ncut
		character(len=*),intent(in)::LegNameRow(:),LegNameCol(:)
		integer::lenName1,lenName2,rank
		lenName1=size(LegNameRow)
		lenName2=size(LegNameCol)
		rank=T%Dimension%getRank()
		if((lenName1+lenName2).ne.rank)then
			call writemess('size(LegNameRow)='+lenName1+'size(LegNameCol)='+lenName2,-1)
			call writemess('rank='+rank,-1)
			call writemess('rank should be equal to size(LegNameRow)+size(LegNameCol)')
			call error_stop
		end if
		if(lenName1.le.lenName2)then
			call T%forward(LegNameRow)
		else
			call T%backWard(LegNameCol)
		end if
		call T%fuse(1,lenName1)
		call T%fuse(2,lenName2+1)
		call SVDcutoff_kill_inData(T,U,s,V,Ncut)
		call U%split()
		call V%split()
		call U%setName(U%Dimension%getRank(),'SVD.U')
		call V%setName(1,'SVD.V')
		return
	end subroutine

	!***********************************************************************
	!***********************************************************************
	!
	!                    LQ
	!
	!***********************************************************************
	!***********************************************************************	


	subroutine LQdecomposition1(T,L,Q)
		class(Tensor),target,intent(in)::T
		class(Tensor),target,intent(inout)::L
		class(Tensor),target,intent(inout)::Q
		class(Tensor),allocatable::v,vv,identity,Tau
		type(Dimension)::dimen
		integer :: i,j,M,N,min_MN,classtype,INFO
		class(Tensor),pointer::Tp
		class(Tensor),pointer::Lp,Qp
		if(T%Dimension%getRank().ne.2) then
			write(*,*)"ERROR in LQ decomposition"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		endif
		Tp=>T
		Lp=>L
		Qp=>Q
		if(associated(Tp,Lp).or.associated(Tp,Qp).or.associated(Lp,Qp))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%LQTensor(L,Q)')
			call writemess('T, L and V can not be a same variable')
			call error_stop
		end if
		Tp=>null()
		Lp=>null()
		Qp=>null()

		allocate(v,mold=T)
		allocate(vv,mold=T)
		allocate(identity,mold=T)
		allocate(Tau,mold=T)
		M = T.dim.1
		N = T.dim.2
		min_MN=min(M,N)
		
		call L%empty()
		call Q%empty()
		classtype=max(2,T%getType())
		call L%setType(classtype)
		L=T
		call Tau%allocate((/min_MN/),classtype)
		INFO=999
		call TData_LQ(Tau%TData,L%TData,M,N,INFO)
		
		if(info.ne.0) then
			call writemess('Error in LQ decomposition ,info='+info,-1)
			call writemess('output The data in ./_LQ_ERROR_LOG.err',-1)
			open(unit=9991,file='./_LQ_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
			call T%writeinfo('The Matrix in LQ',9991)
			close(9991)
			call error_stop()
		end if
		call v%setType(classtype)
		call vv%setType(classtype)
		call identity%setType(classtype)
		if(M.gt.N)then
			!compute Q
			call v%empty()
			call v%allocate([N],classtype)
			call v%zero()
			identity=eye(N,N)
			do i=1,min_MN
				if(i.ne.1)then
					call v%setValue(i-1,0)
				end if
				call v%setValue(i,dcmplx(1d0,0d0))
				do j=i+1,N
					call v%setValue((/j/),L.i.(/i,j/))
				end do
				vv=tau%i(i)*((.h.v).xx.v)
				if(i.eq.1)then
					Q=identity-vv
				else
					Q=Q*(identity-vv)
				end if
			end do
			Q=.h.Q
			Q=Q%subTensor((/-1,1,min_MN/))
			dimen=(/min_MN/)+(T%dimension.subdim.2)
			Q%Dimension=dimen
			!compute L
			do i=1,M
				do j=i+1,min_MN
					call L%setvalue((/i,j/),0)
				end do
			end do
			L=L%subTensor((/-2,1,min_MN/))
			dimen=(T%dimension.subdim.1)+(/min_MN/)
			L%dimension=dimen
		else
			Q=L
			!compute L
			do i=1,M
				do j=i+1,min_MN
					call L%setvalue((/i,j/),0)
				end do
			end do
			L=L%subTensor((/-2,1,min_MN/))
			dimen=(T%Dimension.subdim.1)+(/min_MN/)
			L%Dimension=dimen
			!compute Q
			INFO=999
			call TData_ORGLQ(Tau%TData,Q%TData,M,N,min_MN,INFO)
			dimen=(/min_MN/)+(T.subdim.2)
			Q%Dimension=dimen
			if(info.ne.0) then
				call writemess('Error in LQ decomposition ,info='+info,-1)
				call writemess('output The data in ./_LQ_ERROR_LOG.err',-1)
				open(unit=9991,file='./_LQ_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
				call T%writeinfo('The Matrix in LQ',9991)
				call Q%writeinfo('The Q Matrix in LQ',9991)
				call L%writeinfo('The L Matrix in LQ',9991)
				close(9991)
				call error_stop()
			end if
		end if
		return
	end subroutine

	subroutine LQdecomposition_kill_inData(Q,L)
		class(Tensor),target,intent(inout)::Q
		class(Tensor),target,intent(inout)::L
		class(Tensor),allocatable::v,vv,identity,Tau
		type(Dimension)::dimen
		integer :: i,j,M,N,min_MN,classtype,INFO
		class(Tensor),pointer::Qp
		class(Tensor),pointer::Lp
		real*8,pointer::dp(:,:),Qdp(:,:)
		real*4,pointer::sp(:,:),Qsp(:,:)
		complex(kind=4),pointer::cp(:,:),Qcp(:,:)
		complex(kind=8),pointer::zp(:,:),Qzp(:,:)
		type(Dimension)::TenDim
		if(Q%Dimension%getRank().ne.2) then
			write(*,*)"ERROR in LQ decomposition"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		endif
		Lp=>L
		Qp=>Q
		if(associated(Lp,Qp))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%LQTensor(L,Q)')
			call writemess('T, L and V can not be a same variable')
			call error_stop
		end if
		Lp=>null()
		Qp=>null()
		allocate(v,mold=Q)
		allocate(vv,mold=Q)
		allocate(identity,mold=Q)
		allocate(Tau,mold=Q)

		TenDim=Q%Dimension
		M = Q.dim.1
		N = Q.dim.2
		min_MN=min(M,N)
		
		classtype=max(2,Q%getType())
		call Tau%allocate((/min_MN/),classtype)
		INFO=999
		call TData_LQ(Tau%TData,Q%TData,M,N,INFO)
		
		if(info.ne.0) then
			call writemess('Error in LQ decomposition ,info='+info,-1)
			call error_stop()
		end if
		
		call L%empty()
		call L%allocate((TenDim.subdim.1)+(/min_MN/),classtype)
		


		!compute L
		select case(classtype)
			case(2)
				call L%pointer(sp)
				call Q%pointer(Qsp,[1,M],[1,min_MN])
				if(M.le.N)then
					do i=1,min_MN
						sp(i:M,i)=Qsp(i:M,i)
						if(i.gt.1)sp(1:i-1,i)=0
					end do
					INFO=999
					call TData_ORGLQ(Tau%TData,Q%TData,M,N,min_MN,INFO)
					call Q%reset_dim_no_check((/min_MN/)+(Tendim.subdim.2))
				else
					call scopy(Q%TData%getTotalData(),Qsp,1,sp,1)
					do i=2,min_MN
						sp(1:i-1,i)=0
					end do
					INFO=999
					call TData_ORGLQ(Tau%TData,Q%TData,M,N,min_MN,INFO)
					call WorkingMemory%check()
					call WorkingMemory%get_memory(sp,min_MN,N)
					sp=Qsp(1:min_MN,1:N)
					call Q%reset_dim_no_check((/min_MN/)+(Tendim.subdim.2))
					call Q%pointer(Qsp)
					Qsp=sp
					call WorkingMemory%free()
				end if
			case(3)
				call L%pointer(dp)
				call Q%pointer(Qdp,[1,M],[1,min_MN])
				if(M.le.N)then
					do i=1,min_MN
						dp(i:M,i)=Qdp(i:M,i)
						if(i.gt.1)dp(1:i-1,i)=0
					end do
					INFO=999
					call TData_ORGLQ(Tau%TData,Q%TData,M,N,min_MN,INFO)
					call Q%reset_dim_no_check((/min_MN/)+(Tendim.subdim.2))
				else
					call dcopy(Q%TData%getTotalData(),Qdp,1,dp,1)
					do i=2,min_MN
						dp(1:i-1,i)=0
					end do
					INFO=999
					call TData_ORGLQ(Tau%TData,Q%TData,M,N,min_MN,INFO)
					call WorkingMemory%check()
					call WorkingMemory%get_memory(dp,min_MN,N)
					dp=Qdp(1:min_MN,1:N)
					call Q%reset_dim_no_check((/min_MN/)+(Tendim.subdim.2))
					call Q%pointer(Qdp)
					Qdp=dp
					call WorkingMemory%free()
				end if
			case(4)
				call L%pointer(cp)
				call Q%pointer(Qcp,[1,M],[1,min_MN])
				if(M.le.N)then
					do i=1,min_MN
						cp(i:M,i)=Qcp(i:M,i)
						if(i.gt.1)cp(1:i-1,i)=0
					end do
					INFO=999
					call TData_ORGLQ(Tau%TData,Q%TData,M,N,min_MN,INFO)
					call Q%reset_dim_no_check((/min_MN/)+(Tendim.subdim.2))
				else
					call ccopy(Q%TData%getTotalData(),Qcp,1,cp,1)
					do i=2,min_MN
						cp(1:i-1,i)=0
					end do
					INFO=999
					call TData_ORGLQ(Tau%TData,Q%TData,M,N,min_MN,INFO)
					call WorkingMemory%check()
					call WorkingMemory%get_memory(cp,min_MN,N)
					cp=Qcp(1:min_MN,1:N)
					call Q%reset_dim_no_check((/min_MN/)+(Tendim.subdim.2))
					call Q%pointer(Qcp)
					Qcp=cp
					call WorkingMemory%free()
				end if
			case(5)
				call L%pointer(zp)
				call Q%pointer(Qzp,[1,M],[1,min_MN])
				if(M.le.N)then
					do i=1,min_MN
						zp(i:M,i)=Qzp(i:M,i)
						if(i.gt.1)zp(1:i-1,i)=0
					end do
					INFO=999
					call TData_ORGLQ(Tau%TData,Q%TData,M,N,min_MN,INFO)
					call Q%reset_dim_no_check((/min_MN/)+(Tendim.subdim.2))
				else
					call zcopy(Q%TData%getTotalData(),Qzp,1,zp,1)
					do i=2,min_MN
						zp(1:i-1,i)=0
					end do
					INFO=999
					call TData_ORGLQ(Tau%TData,Q%TData,M,N,min_MN,INFO)
					call WorkingMemory%check()
					call WorkingMemory%get_memory(zp,min_MN,N)
					zp=Qzp(1:min_MN,1:N)
					call Q%reset_dim_no_check((/min_MN/)+(Tendim.subdim.2))
					call Q%pointer(Qzp)
					Qzp=zp
					call WorkingMemory%free()
				end if
		end select
		if(info.ne.0) then
			call writemess('Error in LQ decomposition ,info='+info,-1)
			call error_stop()
		end if
		return
	end subroutine
	
	function LQTensor_noName(T)	result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),target,intent(in)::T
		allocate(res(2),mold=T)
		call LQdecomposition1(T,Res(1),Res(2))
		return
	end function
	
	function LQTensor_name(inputT,nameU,nameV)	result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),target,intent(in)::inputT
		character(len=*),intent(in)::nameU,nameV
		type(Tensor)::T
		integer::rank,rankU,rankV,i
		T=inputT
		rank=T%Dimension%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if((T%Dimension%getName(i).subl.indexsymbol).equ.nameU) rankU=rankU+1
			if((T%Dimension%getName(i).subl.indexsymbol).equ.nameV) rankV=rankV+1
		end do
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in LQTensor_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in LQTensor_name,no such name",-1)
			call writemess(nameU,-1)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in LQTensor_name,no such name",-1)
			call writemess(nameV,-1)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			call T%forward(nameU)
		else
			call T%backward(nameV)
		end if
		call T%fuse(1,rankU)
		call T%fuse(2,rankV+1)
		allocate(res(2),mold=T)
		call LQdecomposition1(T,Res(1),Res(2))
		call Res(1)%split()
		call Res(2)%split()

		call Res(1)%setName(Res(1)%Dimension%getRank(),'LQ.L')
		call Res(2)%setName(1,'LQ.Q')
		return
	end function
	
	function LQTensorNameLeft(T,LegName,Left_)result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),intent(in)::T
		logical,optional,intent(in)::Left_
		character(len=*),intent(in)::LegName(:)
		class(Tensor),allocatable::Temp
		logical::Left
		integer::lenName,rank
		if(present(Left_))then
			Left=Left_
		else
			Left=.true.
		end if
		allocate(res(2),mold=T)
		allocate(Temp,mold=T)
		lenName=size(LegName)
		rank=T%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			Temp=T.pf.LegName
			call Temp%fuse(1,lenName)
			call Temp%fuse(2,rank)
			call LQdecomposition1(Temp,res(1),res(2))
		else
			Temp=T.pb.LegName
			call Temp%fuse(1,rank-lenName)
			call Temp%fuse(2,Temp%Dimension%getRank())
			call LQdecomposition1(Temp,res(1),res(2))
		end if
		call res(1)%split()
		call res(2)%split()
		call Res(1)%setName(Res(1)%Dimension%getRank(),'LQ.L')
		call Res(2)%setName(1,'LQ.Q')
		return
	end function
	
	subroutine LQTensorNameLeftRoutine(T,L,Q,LegName,Left_)
		class(Tensor),intent(in)::T
		class(Tensor),intent(inout)::L
		class(Tensor),intent(inout)::Q
		logical,optional,intent(in)::Left_
		character(len=*),intent(in)::LegName(:)
		logical::Left
		class(Tensor),allocatable::Temp
		integer::lenName,rank
		if(present(Left_))then
			Left=Left_
		else
			Left=.true.
		end if
		allocate(temp,mold=T)
		lenName=size(LegName)
		rank=T%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			Temp=T.pf.LegName
			call Temp%fuse(1,lenName)
			call Temp%fuse(2,rank)
			call LQdecomposition1(Temp,L,Q)
		else
			Temp=T.pb.LegName
			call Temp%fuse(1,rank-lenName)
			call Temp%fuse(2,Temp%Dimension%getRank())
			call LQdecomposition1(Temp,L,Q)
		end if
		call L%split()
		call Q%split()
		call L%setName(L%Dimension%getRank(),'LQ.L')
		call Q%setName(1,'LQ.Q')
		return
	end subroutine

	subroutine LQTensorNameLeftRoutine_kill_inData(Q,L,LegName,Left_)
		class(Tensor),intent(inout)::Q
		class(Tensor),intent(inout)::L
		logical,optional,intent(in)::Left_
		character(len=*),intent(in)::LegName(:)
		logical::Left
		integer::lenName,rank
		if(present(Left_))then
			Left=Left_
		else
			Left=.true.
		end if
		lenName=size(LegName)
		rank=Q%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			call Q%forward(LegName)
			call Q%fuse(1,lenName)
			call Q%fuse(2,rank)
			call LQdecomposition_kill_inData(Q,L)
		else
			call Q%backward(LegName)
			call Q%fuse(1,rank-lenName)
			call Q%fuse(2,Q%Dimension%getRank())
			call LQdecomposition_kill_inData(Q,L)
		end if
		call L%split()
		call Q%split()

		call L%setName(L%Dimension%getRank(),'LQ.L')
		call Q%setName(1,'LQ.Q')
		return
	end subroutine


	function LQTensorLegName(T,LegNameRow,LegNameCol)result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::LegNameRow(:),LegNameCol(:)
		class(Tensor),allocatable::Temp
		integer::lenName1,lenName2,rank
		allocate(res(2),mold=T)
		allocate(Temp,mold=T)
		lenName1=size(LegNameRow)
		lenName2=size(LegNameCol)
		rank=T%Dimension%getRank()
		if((lenName1+lenName2).ne.rank)then
			call writemess('size(LegNameRow)='+lenName1+'size(LegNameCol)='+lenName2,-1)
			call writemess('rank='+rank,-1)
			call writemess('rank should be equal to size(LegNameRow)+size(LegNameCol)')
			call error_stop
		end if
		if(lenName1.le.lenName2)then
			Temp=T.pf.LegNameRow
		else
			Temp=T.pb.LegNameCol
		end if
		call Temp%fuse(1,lenName1)
		call Temp%fuse(2,lenName2+1)
		call LQdecomposition1(Temp,res(1),res(2))
		call res(1)%split()
		call res(2)%split()

		call Res(1)%setName(Res(1)%Dimension%getRank(),'LQ.L')
		call Res(2)%setName(1,'LQ.Q')
		return
	end function
	
	subroutine LQTensorLegNameRoutine(T,L,Q,LegNameRow,LegNameCol)
		class(Tensor),intent(inout)::L,Q
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::LegNameRow(:),LegNameCol(:)
		class(Tensor),allocatable::Temp
		integer::lenName1,lenName2,rank
		allocate(Temp,mold=T)
		lenName1=size(LegNameRow)
		lenName2=size(LegNameCol)
		rank=T%Dimension%getRank()
		if((lenName1+lenName2).ne.rank)then
			call writemess('size(LegNameRow)='+lenName1+'size(LegNameCol)='+lenName2,-1)
			call writemess('rank='+rank,-1)
			call writemess('rank should be equal to size(LegNameRow)+size(LegNameCol)')
			call error_stop
		end if
		if(lenName1.le.lenName2)then
			Temp=T.pf.LegNameRow
		else
			Temp=T.pb.LegNameCol
		end if
		call Temp%fuse(1,lenName1)
		call Temp%fuse(2,lenName2+1)
		call LQdecomposition1(Temp,L,Q)
		call L%split()
		call Q%split()

		call L%setName(L%Dimension%getRank(),'LQ.L')
		call Q%setName(1,'LQ.Q')
		return
	end subroutine
	subroutine LQTensorLegNameRoutine_kill_inData(Q,L,LegNameRow,LegNameCol)
		class(Tensor),intent(inout)::Q
		class(Tensor),intent(inout)::L
		character(len=*),intent(in)::LegNameRow(:),LegNameCol(:)
		integer::lenName1,lenName2,rank
		lenName1=size(LegNameRow)
		lenName2=size(LegNameCol)
		rank=Q%Dimension%getRank()
		if((lenName1+lenName2).ne.rank)then
			call writemess('size(LegNameRow)='+lenName1+'size(LegNameCol)='+lenName2,-1)
			call writemess('rank='+rank,-1)
			call writemess('rank should be equal to size(LegNameRow)+size(LegNameCol)')
			call error_stop
		end if
		if(lenName1.lt.lenName2)then
			call Q%forward(LegNameRow)
		else
			call Q%backward(LegNameCol)
		end if
		call Q%fuse(1,lenName1)
		call Q%fuse(2,lenName2+1)
		call LQdecomposition_kill_inData(Q,L)
		call L%split()
		call Q%split()

		call L%setName(L%Dimension%getRank(),'LQ.L')
		call Q%setName(1,'LQ.Q')
		return
	end subroutine

	!***********************************************************************
	!***********************************************************************
	!
	!                    LQ
	!
	!***********************************************************************
	!***********************************************************************	

	subroutine QRdecomposition1(T,Q,R) 
		class(Tensor),target,intent(in)::T
		class(Tensor),target,intent(inout)::Q,R
		class(Tensor),allocatable::v,vv,identity,Tau
		type(Dimension)::dimen
		integer :: i,j,M,N,min_MN,classtype,INFO
		class(Tensor),pointer::Tp
		class(Tensor),pointer::Qp,Lp
		if(T%Dimension%getRank().ne.2) then
			write(*,*)"ERROR in QR decomposition"
			write(*,*)"input Tensor should be a matrix"
			call error_stop()
		endif
		Tp=>T
		Lp=>Q
		Qp=>R
		if(associated(Tp,Lp).or.associated(Tp,Qp).or.associated(Lp,Qp))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%QRTensor(Q,R)')
			call writemess('T, Q and R can not be a same variable')
			call error_stop
		end if
		Tp=>null()
		Lp=>null()
		Qp=>null()

		allocate(v,mold=T)
		allocate(vv,mold=T)
		allocate(identity,mold=T)
		allocate(Tau,mold=T)
		M = T.dim.1
		N = T.dim.2
		min_MN=min(M,N)
		INFO=999
		classtype=max(2,T%getType())
		call Q%empty()
		call R%empty()
		call R%setType(classtype)
		R=T
		call Tau%allocate((/min_MN/),classtype)
		call TData_QR(Tau%TData,R%TData,M,N,INFO)
		if(info.ne.0) then
			call writemess('Error in QR decomposition ,info='+info,-1)
			call writemess('output The data in ./_QR_ERROR_LOG.err',-1)
			open(unit=9991,file='./_QR_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
			call T%writeinfo('The Matrix in QR',9991)
			close(9991)
			call error_stop()
		end if
		
		call v%setType(classtype)
		call vv%setType(classtype)
		call identity%setType(classtype)
		
		if((M.lt.N)) then
			!compute Q
			call v%empty()
			call v%allocate([M],classtype)
			call v%zero()
			identity=eye(M,M)
			do i=1,min_MN
				if(i.ne.1)then
					call v%setValue(i-1,0)
				end if
				call v%setValue(i,1)
				do j=i+1,M
					call v%setValue((/j/),R.i.(/j,i/))
				end do
				vv=tau%i(i)*(v.xx.(.h.v))
				if(i.eq.1)then
					Q=identity-vv
				else
					Q=Q*(identity-vv)
				end if
			end do
			Q=Q%subTensor((/-2,1,min_MN/))
			dimen=(T%Dimension.subdim.1)+(/min_MN/)
			Q%Dimension=dimen
			!compute R
			do i=1,min_MN
				do j=1,i-1
					call R%setvalue((/i,j/),0)
				end do
			end do
			R=R%subTensor((/-1,1,min_MN/))
			dimen=(/min_MN/)+(T%Dimension.subdim.2)
			R%Dimension=dimen
		else
			!compute R
			Q=R
			do i=1,min_MN
				do j=1,i-1
					call R%setvalue((/i,j/),0)
				end do
			end do
			R=R%subTensor((/-1,1,min_MN/))
			dimen=(/min_MN/)+(T%Dimension.subdim.2)
			R%Dimension=dimen
			!compute Q
			INFO=999
			call TData_ORGQR(Tau%TData,Q%TData,M,N,min_MN,INFO)
			dimen=(T%Dimension.subdim.1)+(/min_MN/)
			Q%Dimension=dimen
			if(info.ne.0) then
				call writemess('Error in QR decomposition ,info='+info,-1)
				call writemess('output The data in ./_QR_ERROR_LOG.err',-1)
				open(unit=9991,file='./_QR_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
				call T%writeinfo('The Matrix in QR',9991)
				call Q%writeinfo('The Q Matrix in QR',9991)
				call R%writeinfo('The R Matrix in QR',9991)
				close(9991)
				call error_stop()
			end if
		end if
		RETURN
	end subroutine
	subroutine QRdecomposition_kill_inData(Q,R) 
		class(Tensor),target,intent(inout)::Q
		class(Tensor),target,intent(inout)::R
		class(Tensor),allocatable::v,vv,identity,Tau
		type(Dimension)::dimen
		integer :: i,j,M,N,min_MN,MM,classtype,INFO
		class(Tensor),pointer::Rp
		class(Tensor),pointer::Qp
		real*8,pointer::dp(:,:),Qdp(:,:),dvp(:),dworkingvp(:),dtaup(:)
		real*4,pointer::sp(:,:),Qsp(:,:),svp(:),sworkingvp(:),staup(:)
		complex(kind=4),pointer::cp(:,:),Qcp(:,:),cvp(:),cworkingvp(:),ctaup(:)
		complex(kind=8),pointer::zp(:,:),Qzp(:,:),zvp(:),zworkingvp(:),ztaup(:)
		type(Dimension)::TenDim
		if(Q%Dimension%getRank().ne.2) then
			call writemess("ERROR in QR decomposition",-1)
			call writemess("input Tensor should be a matrix",-1)
			call error_stop()
		endif
		TenDim=Q%Dimension
		Qp=>Q
		Rp=>R
		if(associated(Qp,Rp))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%QRTensor(Q,R)',-1)
			call writemess('T, Q and R can not be a same variable',-1)
			call error_stop
		end if
		Rp=>null()
		Qp=>null()

		allocate(v,mold=Q)
		allocate(vv,mold=Q)
		allocate(identity,mold=Q)
		allocate(Tau,mold=Q)
		M = Q.dim.1
		N = Q.dim.2
		min_MN=min(M,N)
		INFO=999
		classtype=max(2,Q%getType())
		call Tau%allocate((/min_MN/),classtype)
		call TData_QR(Tau%TData,Q%TData,M,N,INFO)
		if(info.ne.0) then
			call writemess('Error in QR decomposition ,info='+info,-1)
			call error_stop()
		end if
		
		call R%empty()
		call R%allocate((/min_MN/)+(TenDim.subdim.2),classtype)
		
		select case(classtype)
			case(2)
				call R%pointer(sp)
				call Q%pointer(Qsp,[1,min_MN],[1,N])
				if(M.lt.N)then
					call scopy(Q%TData%getTotalData(),Qsp,1,sp,1)
					do i=1,min_MN-1
						sp(i+1:min_MN,i)=0
					end do
				else
					do i=1,N
						sp(:i,i)=Qsp(:i,i)
						if(i.lt.N)sp(i+1:,i)=0
					end do
				end if
				
			case(3)
				call R%pointer(dp)
				call Q%pointer(Qdp,[1,min_MN],[1,N])
				if(M.lt.N)then
					call dcopy(Q%TData%getTotalData(),Qdp,1,dp,1)
					do i=1,min_MN-1
						dp(i+1:min_MN,i)=0
					end do
				else
					do i=1,min_MN
						dp(:i,i)=Qdp(:i,i)
						if(i.lt.min_MN)dp(i+1:min_MN,i)=0
					end do
				end if
			case(4)
				call R%pointer(cp)
				call Q%pointer(Qcp,[1,min_MN],[1,N])
				if(M.lt.N)then
					call ccopy(Q%TData%getTotalData(),Qcp,1,cp,1)
					do i=1,min_MN-1
						cp(i+1:min_MN,i)=0
					end do
				else
					do i=1,min_MN
						cp(:i,i)=Qcp(:i,i)
						if(i.lt.min_MN)cp(i+1:min_MN,i)=0
					end do
				end if
			case(5)
				call R%pointer(zp)
				call Q%pointer(Qzp,[1,min_MN],[1,N])
				if(M.lt.N)then
					call zcopy(Q%TData%getTotalData(),Qzp,1,zp,1)
					do i=1,min_MN-1
						zp(i+1:min_MN,i)=0
					end do
				else
					do i=1,N
						zp(:i,i)=Qzp(:i,i)
						if(i.lt.N)zp(i+1:N,i)=0
					end do
				end if
				
		end select
		INFO=999
		call TData_ORGQR(Tau%TData,Q%TData,M,N,min_MN,INFO)
		call Q%reset_dim_no_check((TenDim.subdim.1)+[min_MN])
		if(info.ne.0) then
			call writemess('Error in QR decomposition ,info='+info,-1)
			call error_stop()
		end if
		RETURN
	end subroutine


	function QRTensor_noName(T)	result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),target,intent(in)::T
		allocate(res(2),mold=T)
		call QRdecomposition1(T,Res(1),Res(2))
		return
	end function
	function QRTensor_name(inputT,nameU,nameV)	result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),target,intent(in)::inputT
		character(len=*),intent(in)::nameU,nameV
		class(Tensor),allocatable::T
		integer::rank,rankU,rankV,i
		allocate(T,mold=inputT)
		T=inputT
		rank=T%Dimension%getRank()
		rankU=0
		rankV=0
		do i=1,rank
			if((T%Dimension%getName(i).subl.indexsymbol).equ.nameU) rankU=rankU+1
			if((T%Dimension%getName(i).subl.indexsymbol).equ.nameV) rankV=rankV+1
		end do
		if(rankU+rankV.ne.rank) then
			call writemess("ERROR in QRTensor_name",-1)
			call writemess(rankU+','+rankV+','+rank,-1)
			call error_stop()
		end if
		if(rankU.eq.0) then
			call writemess("ERROR in QRTensor_name,no such name",-1)
			call writemess(nameU,-1)
			call error_stop()
		end if
		if(rankV.eq.0) then
			call writemess("ERROR in QRTensor_name,no such name",-1)
			call writemess(nameV,-1)
			call error_stop()
		end if
		if(rankU.le.rankV)then
			call T%forward(nameU)
		else
			call T%backward(nameV)
		end if
		call T%fuse(1,rankU)
		call T%fuse(2,rankV+1)
		allocate(res(2),mold=inputT)
		call QRdecomposition1(T,Res(1),Res(2))
		call Res(1)%split()
		call Res(2)%split()

		call Res(1)%setName(Res(1)%Dimension%getRank(),'QR.Q')
		call Res(2)%setName(1,'QR.R')
		return
	end function
	
	
	function QRTensorNameLeft(T,LegName,Left_)result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),intent(in)::T
		logical,optional,intent(in)::Left_
		character(len=*),intent(in)::LegName(:)
		class(Tensor),allocatable::Temp
		logical::Left
		integer::lenName,rank
		if(present(Left_))then
			Left=Left_
		else
			Left=.true.
		end if
		allocate(res(2),mold=T)
		allocate(Temp,mold=T)
		lenName=size(LegName)
		rank=T%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			Temp=T.pf.LegName
			call Temp%fuse(1,lenName)
			call Temp%fuse(2,rank)
			call QRdecomposition1(Temp,res(1),res(2))
		else 
			Temp=T.pb.LegName
			call Temp%fuse(1,rank-lenName)
			call Temp%fuse(2,Temp%Dimension%getRank())
			call QRdecomposition1(Temp,res(1),res(2))
		end if
		call res(1)%split()
		call res(2)%split()

		call Res(1)%setName(Res(1)%Dimension%getRank(),'QR.Q')
		call Res(2)%setName(1,'QR.R')
		return
	end function

	
	subroutine QRTensorNameLeftRoutine(T,Q,R,LegName,Left_)
		class(Tensor),intent(in)::T
		class(Tensor),intent(inout)::Q,R
		logical,optional,intent(in)::Left_
		character(len=*),intent(in)::LegName(:)
		logical::Left
		class(Tensor),allocatable::Temp
		integer::lenName,rank
		if(present(Left_))then
			Left=Left_
		else
			Left=.true.
		end if
		allocate(temp,mold=T)
		lenName=size(LegName)
		rank=T%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			Temp=T.pf.LegName
			call Temp%fuse(1,lenName)
			call Temp%fuse(2,rank)
			call QRdecomposition1(Temp,Q,R)
		else 
			Temp=T.pb.LegName
			call Temp%fuse(1,rank-lenName)
			call Temp%fuse(2,Temp%Dimension%getRank())
			call QRdecomposition1(Temp,Q,R)
		end if
		call Q%split()
		call R%split()

		call Q%setName(Q%Dimension%getRank(),'QR.Q')
		call R%setName(1,'QR.R')
		return
	end subroutine

	subroutine QRTensorNameLeftRoutine_kill_indaTa(Q,R,LegName,Left_)
		class(Tensor),intent(inout)::Q
		class(Tensor),intent(inout)::R
		logical,optional,intent(in)::Left_
		character(len=*),intent(in)::LegName(:)
		logical::Left
		integer::lenName,rank
		if(present(Left_))then
			Left=Left_
		else
			Left=.true.
		end if
		lenName=size(LegName)
		rank=Q%Dimension%getRank()
		if(lenName.ge.rank)then
			call writemess('ERROR in the number of legs')
			call writemess('number of leg can not larger than of equal to the rank of the tensor')
			call error_stop
		end if
		if(Left)then
			call Q%forward(LegName)
			call Q%fuse(1,lenName)
			call Q%fuse(2,rank)
			call QRdecomposition_kill_inData(Q,R)
		else
			call Q%backWard(LegName)
			call Q%fuse(1,rank-lenName)
			call Q%fuse(2,Q%Dimension%getRank())
			call QRdecomposition_kill_inData(Q,R)
		end if
		call Q%split()
		call R%split()

		call Q%setName(Q%Dimension%getRank(),'QR.Q')
		call R%setName(1,'QR.R')
		return
	end subroutine


	function QRTensorLegName(T,LegNameRow,LegNameCol)result(res)
		class(Tensor),allocatable::res(:)
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::LegNameRow(:),LegNameCol(:)
		class(Tensor),allocatable::Temp
		integer::lenName1,lenName2,rank
		allocate(res(2),mold=T)
		allocate(Temp,mold=T)
		lenName1=size(LegNameRow)
		lenName2=size(LegNameCol)
		rank=T%Dimension%getRank()
		if((lenName1+lenName2).ne.rank)then
			call writemess('size(LegNameRow)='+lenName1+'size(LegNameCol)='+lenName2,-1)
			call writemess('rank='+rank,-1)
			call writemess('rank should be equal to size(LegNameRow)+size(LegNameCol)')
			call error_stop
		end if
		if(lenName1.le.lenName2)then
			Temp=T.pf.LegNameRow
		else
			Temp=T.pb.LegNameCol
		end if
		call Temp%fuse(1,lenName1)
		call Temp%fuse(2,lenName2+1)
		call QRdecomposition1(Temp,res(1),res(2))
		call res(1)%split()
		call res(2)%split()

		call Res(1)%setName(Res(1)%Dimension%getRank(),'QR.Q')
		call Res(2)%setName(1,'QR.R')
		return
	end function
	
	subroutine QRTensorLegNameRoutine(T,Q,R,LegNameRow,LegNameCol)
		class(Tensor),intent(inout)::Q,R
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::LegNameRow(:),LegNameCol(:)
		class(Tensor),allocatable::Temp
		integer::lenName1,lenName2,rank
		allocate(temp,mold=T)
		lenName1=size(LegNameRow)
		lenName2=size(LegNameCol)
		rank=T%Dimension%getRank()
		if((lenName1+lenName2).ne.rank)then
			call writemess('size(LegNameRow)='+lenName1+'size(LegNameCol)='+lenName2,-1)
			call writemess('rank='+rank,-1)
			call writemess('rank should be equal to size(LegNameRow)+size(LegNameCol)')
			call error_stop
		end if
		if(lenName1.le.lenName2)then
			Temp=T.pf.LegNameRow
		else
			Temp=T.pb.LegNameCol
		end if
		call Temp%fuse(1,lenName1)
		call Temp%fuse(2,lenName2+1)
		call QRdecomposition1(Temp,Q,R)
		call Q%split()
		call R%split()

		call Q%setName(Q%Dimension%getRank(),'QR.Q')
		call R%setName(1,'QR.R')

		return
	end subroutine

	subroutine QRTensorLegNameRoutine_kill_inData(Q,R,LegNameRow,LegNameCol)
		class(Tensor),intent(inout)::Q
		class(Tensor),intent(inout)::R
		character(len=*),intent(in)::LegNameRow(:),LegNameCol(:)
		integer::lenName1,lenName2,rank
		lenName1=size(LegNameRow)
		lenName2=size(LegNameCol)
		rank=Q%Dimension%getRank()
		if((lenName1+lenName2).ne.rank)then
			call writemess('size(LegNameRow)='+lenName1+'size(LegNameCol)='+lenName2,-1)
			call writemess('rank='+rank,-1)
			call writemess('rank should be equal to size(LegNameRow)+size(LegNameCol)')
			call error_stop
		end if
		if(lenName1.lt.lenName2)then
			call Q%forward(LegNameRow)
		else
			call Q%backward(LegNameCol)
		end if
		call Q%fuse(1,lenName1)
		call Q%fuse(2,lenName2+1)
		call QRdecomposition_kill_inData(Q,R)
		call Q%split()
		call R%split()
		call Q%setName(Q%Dimension%getRank(),'QR.Q')
		call R%setName(1,'QR.R')
		return
	end subroutine




	!**************************************************************************************************************
	!**************************************************************************************************************
	!
	!                                  contract
	!
	!**************************************************************************************************************
	!**************************************************************************************************************	

	!******************  contract  *********************
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
		!	if present len_of_contract, len_of_contract specify the length of  i1, and i2

	function contract_noName(T1_,i1,T2_,i2,len_of_contract) result(T)
		class(Tensor),allocatable::T
		class(Tensor),intent(in) :: T1_
		class(Tensor),intent(in) :: T2_
		integer,intent(in) :: i1(:),i2(:)
		integer,optional,intent(in)::len_of_contract(2)
		type(Tensor),pointer :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		allocate(T,mold=T1_)
		T1=>WorkingTensor1
		T2=>WorkingTensor2
		!if(.not.(if_original_dim(T1_%Dimension).and.if_original_dim(T2_%Dimension))) then
		!	write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
		!	write(*,*)"stop"
		!	call error_stop()
		!end if
		rank1=T1_%Dimension%getRank()
		rank2=T2_%Dimension%getRank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(i1))
			leni2=min(len_of_contract(2),size(i2))
		else
			leni1=size(i1)
			leni2=size(i2)
		end if
		T1=T1_.pb.i1(1:leni1)
		call T1%fuse(rank1-leni1+1,rank1)
		T2=T2_.pf.i2(1:leni2)
		call T2%fuse(1,leni2)
		T=T1 * T2
		return
	end function
	function contract_noName2(T1_,i1,T2_,i2) result(T)
		class(Tensor),allocatable::T
		class(Tensor),intent(in) :: T1_
		class(Tensor),intent(in) :: T2_
		integer,intent(in) :: i1,i2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		allocate(T,mold=T1_)
		T = (T1_.pb.i1) * (T2_.pf.i2)
		return
	end function
	function contract_name(T1_,name1,T2_,name2,len_of_contract) result(T)
		class(Tensor),allocatable::T
		class(Tensor),intent(in) :: T1_
		class(Tensor),intent(in) :: T2_
		character(len=*),intent(in)::name1(:),name2(:)
		integer :: i1(size(name1)),i2(size(name2))
		integer,optional,intent(in)::len_of_contract(2)
		type(Tensor),pointer :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(T1_%Dimension%if_original_dim().and.T2_%Dimension%if_original_dim())) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		allocate(T,mold=T1_)
		T1=>WorkingTensor1
		T2=>WorkingTensor2
		i1=T1_%Dimension%FindOrder(name1)
		i2=T2_%Dimension%FindOrder(name2)
		rank1=T1_%Dimension%getRank()
		rank2=T2_%Dimension%getRank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(i1))
			leni2=min(len_of_contract(2),size(i2))
		else
			leni1=size(i1)
			leni2=size(i2)
		end if
		T1=T1_.pb.i1(1:leni1)
		call T1%fuse(rank1-leni1+1,rank1)
		T2=T2_.pf.i2(1:leni2)
		call T1%fuse(1,leni2)
		T=T1 * T2
		return
	end function
	function contract_name2(T1_,name1,T2_,name2) result(T)
		class(Tensor),allocatable::T
		class(Tensor),intent(in) :: T1_
		class(Tensor),intent(in) :: T2_
		character(len=*),intent(in)::name1,name2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(T1_%Dimension%if_original_dim().and.T2_%Dimension%if_original_dim())) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		allocate(T,mold=T1_)
		T= (T1_.pb.name1) * (T2_.pf.name2)
		return
	end function
	
	
!******************  contract the same names  *********************	

	subroutine find_same_name(inname1,inname2,SameName,lenofname)
		character(len=*),intent(in)::inname1(:),inname2(:)
		character(len=*),intent(inout)::SameName(:)
		integer,intent(inout)::lenofname
		integer::i,j,k,sizeSameName
		k=0
		lenofname=0
		sizeSameName=size(SameName)
		do i=1,size(inname1)
			do j=1,size(inname2)
				if(inname1(i).equ.inname2(j))then
					k=k+1
					if(k.gt.sizeSameName)then
						call writemess('ERROR in finding the same name, maybe something wrong in the tensorname',-1)
						call writemess('You can diminfo to see if the tensor names are right')
						call writemess('It is not allow to have two or more same name in one tensor')
						call error_stop()
					end if
					SameName(k)=inname1(i)
					lenofname=lenofname+1
				end if
			end do
		end do
		return
	end subroutine
	
	function contract_Same_name(T1_,T2_) result(T)
		class(Tensor),allocatable::T
		class(Tensor),intent(in) :: T1_
		class(Tensor),intent(in) :: T2_
		type(Tensor),pointer::T1,T2
		character(len=characterlen),allocatable::Samename(:),name1(:),name2(:)
		integer::lenofname,rank1,rank2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(T1_%Dimension%if_original_dim().and.T2_%Dimension%if_original_dim())) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		allocate(T,mold=T1_)
		T1=>WorkingTensor1
		T2=>WorkingTensor2
		rank1=T1_%Dimension%getRank()
		rank2=T2_%Dimension%getRank()
		allocate(name1(rank1))
		allocate(name2(rank2))
		allocate(Samename(max(rank1,rank2) ))
		name1=T1_%Dimension%getName()
		name2=T2_%Dimension%getName()
		call find_same_name(name1,name2,SameName,lenofname) 
		T1=T1_.pb.SameName(1:lenofname)
		call T1%fuse(rank1-lenofname+1,rank1)
		T2=T2_.pf.SameName(1:lenofname)
		call T2%fuse(1,lenofname)
		T=T1 * T2
		return
	end function
	
	
	
	function contract_int_name(T1_,i1,T2_,name2,len_of_contract) result(T)
		class(Tensor),allocatable::T
		class(Tensor),intent(in) :: T1_
		class(Tensor),intent(in) :: T2_
		character(len=*),intent(in)::name2(:)
		integer,intent(in)::i1(:)
		integer,optional,intent(in)::len_of_contract(2)
		integer :: i2(size(name2))
		type(Tensor),pointer :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%Dimension%if_original_dim()) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		allocate(T,mold=T1_)
		T1=>WorkingTensor1
		T2=>WorkingTensor2
		i2=T2_%Dimension%FindOrder(name2)
		rank1=T1_%Dimension%getRank()
		rank2=T2_%Dimension%getRank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(i1))
			leni2=min(len_of_contract(2),size(i2))
		else
			leni1=size(i1)
			leni2=size(i2)
		end if
		T1=T1_.pb.i1(1:leni1)
		call T1%fuse(rank1-leni1+1,rank1)
		T2=T2_.pf.i2(1:leni2)
		call T2%fuse(1,leni2)
		T=T1 * T2
		return
	end function
	function contract_int_name2(T1_,i1,T2_,name2) result(T)
		class(Tensor),allocatable::T
		class(Tensor),intent(in) :: T1_
		class(Tensor),intent(in) :: T2_
		character(len=*),intent(in)::name2
		integer,intent(in)::i1
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%Dimension%if_original_dim()) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		allocate(T,mold=T1_)
		T = (T1_.pb.i1) * (T2_.pf.name2)
		return
	end function
	function contract_name_int(T1_,name1,T2_,i2,len_of_contract) result(T)
		class(Tensor),allocatable::T
		class(Tensor),intent(in) :: T1_
		class(Tensor),intent(in) :: T2_
		character(len=*),intent(in)::name1(:)
		integer,intent(in)::i2(:)
		integer :: i1(size(name1))
		integer,optional,intent(in)::len_of_contract(2)
		type(Tensor),pointer :: T1,T2
		integer::leni1,leni2,rank1,rank2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T1_%Dimension%if_original_dim()) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		allocate(T,mold=T1_)
		T1=>WorkingTensor1
		T2=>WorkingTensor2
		i1=T1_%Dimension%FindOrder(name1)
		rank1=T1_%Dimension%getRank()
		rank2=T2_%Dimension%getRank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(i1))
			leni2=min(len_of_contract(2),size(i2))
		else
			leni1=size(i1)
			leni2=size(i2)
		end if
		T1=T1_.pb.i1(1:leni1)
		call T1%fuse(rank1-leni1+1,rank1)
		T2=T2_.pf.i2(1:leni2)
		call T2%fuse(1,leni2)
		T=T1 * T2
		return
	end function
	function contract_name_int2(T1_,name1,T2_,i2) result(T)
		class(Tensor),allocatable::T
		class(Tensor),intent(in) :: T1_
		class(Tensor),intent(in) :: T2_
		character(len=*),intent(in)::name1
		integer,intent(in)::i2
		if(.not.T1_%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2_%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T1_%Dimension%if_original_dim()) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		allocate(T,mold=T1_)
		T = (T1_.pb.name1) * (T2_.pf.i2)
		return
	end function

	subroutine contract_name_routine(T,T1,name1,T2,name2,len_of_contract) 
		class(Tensor),target :: T
		class(Tensor),target :: T1
		class(Tensor),target :: T2
		character(len=*),intent(in)::name1(:),name2(:)
		integer :: i1(size(name1)),i2(size(name2))
		integer,optional,intent(in)::len_of_contract(2)
		integer::leni1,leni2,rank1,rank2
		class(Tensor),pointer::pT,pT1,pT2
		if(.not.T1%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(T1%Dimension%if_original_dim().and.T2%Dimension%if_original_dim())) then
			call writemess("ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		pT=>T
		pT1=>T1
		pT2=>T2
		if(associated(pT,pT1).or.associated(pT,pT2).or.associated(pT1,pT2))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%contract(A,[names],B,[names])')
			call writemess('T, A and B can not be a same variable')
			call error_stop
		end if
		pT=>null()
		pT1=>null()
		pT2=>null()
		i1=T1%Dimension%FindOrder(name1)
		i2=T2%Dimension%FindOrder(name2)
		rank1=T1%Dimension%getRank()
		rank2=T2%Dimension%getRank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(i1))
			leni2=min(len_of_contract(2),size(i2))
		else
			leni1=size(i1)
			leni2=size(i2)
		end if
		call T1%backward(i1(1:leni1))
		call T1%fuse(rank1-leni1+1,rank1)
		call T2%forward(i2(1:leni2))
		call T2%fuse(1,leni2)
		call T%ProductTensorRoutine(T1 , T2 , 1 )
		call T1%split( )
		call T2%split( )
		return
	end subroutine
	subroutine contract_name_routine1(T,name1,T2,name2,len_of_contract) 
		class(Tensor),target::T
		class(Tensor),target :: T2
		character(len=*),intent(in)::name1(:),name2(:)
		integer :: i1(size(name1)),i2(size(name2))
		integer,optional,intent(in)::len_of_contract(2)
		integer::leni1,leni2,rank1,rank2
		class(Tensor),pointer::pT,pT2
		if(.not.T%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(T%Dimension%if_original_dim().and.T2%Dimension%if_original_dim())) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		pT=>T
		pT2=>T2
		if(associated(pT,pT2))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%contract([names],B,[names])')
			call writemess('T and B can not be a same variable')
			call error_stop
		end if
		pT=>null()
		pT2=>null()
		i1=T%Dimension%FindOrder(name1)
		i2=T2%Dimension%FindOrder(name2)
		rank1=T%Dimension%getRank()
		rank2=T2%Dimension%getRank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(i1))
			leni2=min(len_of_contract(2),size(i2))
		else
			leni1=size(i1)
			leni2=size(i2)
		end if
		call T%backward(i1(1:leni1))
		call T%fuse(rank1-leni1+1,rank1)
		call T2%forward(i2(1:leni2))
		call T2%fuse(1,leni2)
		T=T*T2
		call T2%split()
		return
	end subroutine
	subroutine contract_name_routine2(T,T1,name1,name2,len_of_contract) 
		class(Tensor),target::T
		class(Tensor),target:: T1
		character(len=*),intent(in)::name1(:),name2(:)
		integer :: i1(size(name1)),i2(size(name2))
		integer,optional,intent(in)::len_of_contract(2)
		integer::leni1,leni2,rank1,rank2
		class(Tensor),pointer::pT,pT1
		if(.not.T1%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(T1%Dimension%if_original_dim().and.T%Dimension%if_original_dim())) then
			write(*,*)"ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function"
			write(*,*)"stop"
			call error_stop()
		end if
		pT=>T
		pT1=>T1
		if(associated(pT,pT1))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%contract(A,[names],[names])')
			call writemess('T and A can not be a same variable')
			call error_stop
		end if
		pT=>null()
		pT1=>null()
		i1=T1%Dimension%FindOrder(name1)
		i2=T%Dimension%FindOrder(name2)
		rank1=T1%Dimension%getRank()
		rank2=T%Dimension%getRank()
		if(present(len_of_contract))then
			leni1=min(len_of_contract(1),size(i1))
			leni2=min(len_of_contract(2),size(i2))
		else
			leni1=size(i1)
			leni2=size(i2)
		end if
		call T1%backward(i1(1:leni1))
		call T1%fuse(rank1-leni1+1,rank1)
		call T%forward(i2(1:leni2))
		call T%fuse(1,leni2)
		T=T1*T
		call T1%split( )
		return
	end subroutine
	subroutine contract_name_routine4(T,T1,name1,T2,name2) 
		class(Tensor),target::T
		class(Tensor),target :: T1
		class(Tensor),target :: T2
		character(len=*),intent(in)::name1,name2
		class(Tensor),pointer::pT,pT1,pT2
		if(.not.T1%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(T1%Dimension%if_original_dim().and.T2%dimension%if_original_dim())) then
			call writemess("ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		pT=>T
		pT1=>T1
		pT2=>T2
		if(associated(pT,pT1).or.associated(pT,pT2).or.associated(pT1,pT2))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%contract(A,name1,B,name2)')
			call writemess('T, A and B can not be a same variable')
			call error_stop
		end if
		pT=>null()
		pT1=>null()
		pT2=>null()
		call T1%backward(name1)
		call T2%forward(name2)
		call T%ProductTensorRoutine(T1 , T2 , 1 )
		return
	end subroutine
	subroutine contract_name_routine5(T,name1,T2,name2) 
		class(Tensor),target::T
		class(Tensor),target :: T2
		character(len=*),intent(in)::name1,name2
		class(Tensor),pointer::pT,pT2
		if(.not.T%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T2%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(T2%dimension%if_original_dim())) then
			call writemess("ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		pT=>T
		pT2=>T2
		if(associated(pT,pT2))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%contract(name1,B,name2)')
			call writemess('T, A and B can not be a same variable')
			call error_stop
		end if
		pT=>null()
		pT2=>null()
		call T%backward(name1)
		call T2%forward(name2)
		T=T*T2
		return
	end subroutine
	subroutine contract_name_routine6(T,T1,name1,name2) 
		class(Tensor),target::T
		class(Tensor),target:: T1
		character(len=*),intent(in)::name1,name2
		class(Tensor),pointer::pT,pT1
		if(.not.T1%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(T1%Dimension%if_original_dim())) then
			call writemess("ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		pT=>T
		pT1=>T1
		if(associated(pT,pT1))then
			call writemess('input Tensors can not be the same variable',-1)
			call writemess('error in call T%contract(A,name1,name2)')
			call writemess('T, A and B can not be a same variable')
			call error_stop
		end if
		pT=>null()
		pT1=>null()
		call T1%backward(name1)
		call T%forward(name2)
		T=T1*T
		return
	end subroutine

	subroutine contract_name_ownlegs_routine(T,name1,name2) 
		class(Tensor)::T
		type(Tensor),pointer::pT
		character(len=*),intent(in)::name1,name2
		type(Dimension)::NewDimen
		integer::rank,classtype,dim1,i,k
		integer,pointer::idata(:,:,:),newidata(:)
		real*4,pointer::sdata(:,:,:),newsdata(:)
		real*8,pointer::ddata(:,:,:),newddata(:)
		complex*8,pointer::cdata(:,:,:),newcdata(:)
		complex*16,pointer::zdata(:,:,:),newzdata(:)
		if(.not.T%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(T%Dimension%if_original_dim())) then
			call writemess("ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		rank=T%Dimension%getRank()
		if(rank.eq.2)then
			T=T%trace()
			return
		end if
		pT=>WorkingTensor1
		pT=T.pf.name1
		call pT%forward(name2)
		dim1=pT%dim(1)
		if(dim1.ne.pT%dim(2))then
			call writemess(' ERROR in contract(name1,name2), dimension')
			call error_stop
		end if
		NewDimen=pT%Dimension.subdim.[3,rank]
		classtype=T%getType()
		call T%empty()
		call T%allocate(NewDimen,classtype)
		call pT%fuse(3,rank)

		select case (classtype)
			case (1)
				call pT%pointer(idata)
				call T%pointer(newidata)
				do k=1,T%TData%getTotalData()
					newidata(k)=0
					do i=1,dim1
						newidata(k)=newidata(k)+idata(i,i,k)
					end do
				end do
			case (2)
				call pT%pointer(sdata)
				call T%pointer(newsdata)
				do k=1,T%TData%getTotalData()
					newsdata(k)=0
					do i=1,dim1
						newsdata(k)=newsdata(k)+sdata(i,i,k)
					end do
				end do
			case(3)
				call pT%pointer(ddata)
				call T%pointer(newddata)
				do k=1,T%TData%getTotalData()
					newddata(k)=0
					do i=1,dim1
						newddata(k)=newddata(k)+ddata(i,i,k)
					end do
				end do
			case(4)
				call pT%pointer(cdata)
				call T%pointer(newcdata)
				do k=1,T%TData%getTotalData()
					newcdata(k)=0
					do i=1,dim1
						newcdata(k)=newcdata(k)+cdata(i,i,k)
					end do
				end do
			case(5)
				call pT%pointer(zdata)
				call T%pointer(newzdata)
				do k=1,T%TData%getTotalData()
					newzdata(k)=0
					do i=1,dim1
						newzdata(k)=newzdata(k)+zdata(i,i,k)
					end do
				end do
			case default
				call writemess(' ERROR in contract(name1,name2), clasee type')
				call error_stop
		end select
		return
	end subroutine
	subroutine contract_ownlegs_routine(T,ith1,ith2) 
		class(Tensor)::T
		type(Tensor),pointer::pT
		integer,intent(in)::ith1,ith2
		type(Dimension)::NewDimen
		integer::rank,classtype,dim1,i,k
		integer,pointer::idata(:,:,:),newidata(:)
		real*4,pointer::sdata(:,:,:),newsdata(:)
		real*8,pointer::ddata(:,:,:),newddata(:)
		complex*8,pointer::cdata(:,:,:),newcdata(:)
		complex*16,pointer::zdata(:,:,:),newzdata(:)
		if(.not.T%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		rank=T%Dimension%getRank()
		if(rank.eq.2)then
			T=T%trace()
			return
		end if
		pT=>WorkingTensor1
		pT=T.pf.ith1
		call pT%forward(ith2)
		dim1=pT%dim(1)
		if(dim1.ne.pT%dim(2))then
			call writemess(' ERROR in contract(ith1,ith2), dimension')
			call error_stop
		end if
		NewDimen=pT.subdim.[3,rank]
		classtype=T%getType()
		call T%empty()
		call T%allocate(NewDimen,classtype)
		call pT%fuse(3,rank)

		select case (classtype)
			case (1)
				call pT%pointer(idata)
				call T%pointer(newidata)
				do k=1,T%TData%getTotalData()
					newidata(k)=0
					do i=1,dim1
						newidata(k)=newidata(k)+idata(i,i,k)
					end do
				end do
			case (2)
				call pT%pointer(sdata)
				call T%pointer(newsdata)
				do k=1,T%TData%getTotalData()
					newsdata(k)=0
					do i=1,dim1
						newsdata(k)=newsdata(k)+sdata(i,i,k)
					end do
				end do
			case(3)
				call pT%pointer(ddata)
				call T%pointer(newddata)
				do k=1,T%TData%getTotalData()
					newddata(k)=0
					do i=1,dim1
						newddata(k)=newddata(k)+ddata(i,i,k)
					end do
				end do
			case(4)
				call pT%pointer(cdata)
				call T%pointer(newcdata)
				do k=1,T%TData%getTotalData()
					newcdata(k)=0
					do i=1,dim1
						newcdata(k)=newcdata(k)+cdata(i,i,k)
					end do
				end do
			case(5)
				call pT%pointer(zdata)
				call T%pointer(newzdata)
				do k=1,T%TData%getTotalData()
					newzdata(k)=0
					do i=1,dim1
						newzdata(k)=newzdata(k)+zdata(i,i,k)
					end do
				end do
			case default
				call writemess(' ERROR in contract(name1,name2), clasee type')
				call error_stop
		end select
		return
	end subroutine
	function contract_ownlegs(Tin,ith1,ith2) Result(Res)
		class(Tensor),allocatable::Res
		class(Tensor),intent(in)::Tin
		type(Tensor),pointer::pT
		integer,intent(in)::ith1,ith2
		type(Dimension)::NewDimen
		integer::rank,classtype,dim1,i,k
		integer,pointer::idata(:,:,:),newidata(:)
		real*4,pointer::sdata(:,:,:),newsdata(:)
		real*8,pointer::ddata(:,:,:),newddata(:)
		complex*8,pointer::cdata(:,:,:),newcdata(:)
		complex*16,pointer::zdata(:,:,:),newzdata(:)
		if(.not.Tin%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		allocate(Res,mold=Tin)
		rank=Tin%Dimension%getRank()
		if(rank.eq.2)then
			Res=Tin%trace()
			return
		end if
		pT=>WorkingTensor1
		pT=Tin.pf.ith1
		call pT%forward(ith2)
		dim1=pT%dim(1)
		if(dim1.ne.pT%dim(2))then
			call writemess(' ERROR in contract(ith1,ith2), dimension')
			call error_stop
		end if
		NewDimen=pT.subdim.[3,rank]
		classtype=Tin%getType()
		call Res%empty()
		call Res%allocate(NewDimen,classtype)
		call pT%fuse(3,rank)

		select case (classtype)
			case (1)
				call pT%pointer(idata)
				call Res%pointer(newidata)
				do k=1,Res%TData%getTotalData()
					newidata(k)=0
					do i=1,dim1
						newidata(k)=newidata(k)+idata(i,i,k)
					end do
				end do
			case (2)
				call pT%pointer(sdata)
				call Res%pointer(newsdata)
				do k=1,Res%TData%getTotalData()
					newsdata(k)=0
					do i=1,dim1
						newsdata(k)=newsdata(k)+sdata(i,i,k)
					end do
				end do
			case(3)
				call pT%pointer(ddata)
				call Res%pointer(newddata)
				do k=1,Res%TData%getTotalData()
					newddata(k)=0
					do i=1,dim1
						newddata(k)=newddata(k)+ddata(i,i,k)
					end do
				end do
			case(4)
				call pT%pointer(cdata)
				call Res%pointer(newcdata)
				do k=1,Res%TData%getTotalData()
					newcdata(k)=0
					do i=1,dim1
						newcdata(k)=newcdata(k)+cdata(i,i,k)
					end do
				end do
			case(5)
				call pT%pointer(zdata)
				call Res%pointer(newzdata)
				do k=1,Res%TData%getTotalData()
					newzdata(k)=0
					do i=1,dim1
						newzdata(k)=newzdata(k)+zdata(i,i,k)
					end do
				end do
			case default
				call writemess(' ERROR in contract(name1,name2), clasee type')
				call error_stop
		end select
		return
	end function

	function contract_name_ownlegs(Tin,name1,name2) Result(Res)
		class(Tensor),allocatable::Res
		class(Tensor),intent(in)::Tin
		type(Tensor),pointer::pT
		character(len=*),intent(in)::name1,name2
		type(Dimension)::NewDimen
		integer::rank,classtype,dim1,i,k
		integer,pointer::idata(:,:,:),newidata(:)
		real*4,pointer::sdata(:,:,:),newsdata(:)
		real*8,pointer::ddata(:,:,:),newddata(:)
		complex*8,pointer::cdata(:,:,:),newcdata(:)
		complex*16,pointer::zdata(:,:,:),newzdata(:)
		if(.not.Tin%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.(Tin%Dimension%if_original_dim())) then
			call writemess("ERROR in contract with TensorName, one can not fuse any legs of the Tensor to use this function",-1)
			call writemess("stop",-1)
			call error_stop()
		end if
		allocate(Res,mold=Tin)
		rank=Tin%Dimension%getRank()
		if(rank.eq.2)then
			Res=Tin%trace()
			return
		end if
		pT=>WorkingTensor1
		pT=Tin.pf.name1
		call pT%forward(name2)
		dim1=pT%dim(1)
		if(dim1.ne.pT%dim(2))then
			call writemess(' ERROR in contract(ith1,ith2), dimension')
			call error_stop
		end if
		NewDimen=pT.subdim.[3,rank]
		classtype=Tin%getType()
		call Res%empty()
		call Res%allocate(NewDimen,classtype)
		call pT%fuse(3,rank)

		select case (classtype)
			case (1)
				call pT%pointer(idata)
				call Res%pointer(newidata)
				do k=1,Res%TData%getTotalData()
					newidata(k)=0
					do i=1,dim1
						newidata(k)=newidata(k)+idata(i,i,k)
					end do
				end do
			case (2)
				call pT%pointer(sdata)
				call Res%pointer(newsdata)
				do k=1,Res%TData%getTotalData()
					newsdata(k)=0
					do i=1,dim1
						newsdata(k)=newsdata(k)+sdata(i,i,k)
					end do
				end do
			case(3)
				call pT%pointer(ddata)
				call Res%pointer(newddata)
				do k=1,Res%TData%getTotalData()
					newddata(k)=0
					do i=1,dim1
						newddata(k)=newddata(k)+ddata(i,i,k)
					end do
				end do
			case(4)
				call pT%pointer(cdata)
				call Res%pointer(newcdata)
				do k=1,Res%TData%getTotalData()
					newcdata(k)=0
					do i=1,dim1
						newcdata(k)=newcdata(k)+cdata(i,i,k)
					end do
				end do
			case(5)
				call pT%pointer(zdata)
				call Res%pointer(newzdata)
				do k=1,Res%TData%getTotalData()
					newzdata(k)=0
					do i=1,dim1
						newzdata(k)=newzdata(k)+zdata(i,i,k)
					end do
				end do
			case default
				call writemess(' ERROR in contract(name1,name2), clasee type')
				call error_stop
		end select
		return
	end function



end module