include '/home/sjdong/Tensor/MPI/Tensor.f90' 
program aaa
	use Tensor_complex
	implicit none
	real*8::time1,time2,num1,num2
	real*8,external::omp_get_wtime
	type(Dimension)::dimen,dimen2
	CHARACTER*10::w
	type(Tensor)::T,T2,T3
	integer::ind(2,2),inde,i,j
	complex*16,allocatable::TData(:,:),T2Data(:,:)
	integer::name11(2,1),name12(2,2),name21(2,1),name22(2,2)
	call set_max_len_of_cha(10)
	T2=generate((/2,2,10,10/))
	T3=generate((/2,2,10,10/))
	call setTensorName(T2,'A')
	call setTensorName(T3,'B')
	write(*,*)"Test TensorName for product input array index"
	time1=omp_get_wtime()
	do i=1,10000
		T=Tenproduct(T2,(/'A.'+3,'A.'+4/),T3,(/'B.'+3,'B.'+4/))
	end do
	time2=omp_get_wtime()
	write(*,*)"Char Name",time2-time1
	call cleanTensorName(T2)
	call cleanTensorName(T3)
	call setTensorName(T2,3,(/1/),(/1,1/))
	call setTensorName(T3,3,(/2/),(/1,1/))
	call setTensorName(T2,4,(/1/),(/1,2/))
	call setTensorName(T3,4,(/2/),(/1,2/))
	time1=omp_get_wtime()
	do i=1,10000
		name11=1
		name12(1,:)=(/1,1/)
		name12(2,:)=(/1,2/)
		name21=2
		name22(1,:)=(/1,1/)
		name22(2,:)=(/1,2/)
		!T=Tenproduct(T2,(/1/),(/1,1/),T3,(/2/),(/1,1/))
		T=Tenproduct(T2,name11,name12,T3,name21,name22)
	end do
	time2=omp_get_wtime()
	write(*,*)"Int name",time2-time1
	time1=omp_get_wtime()
	do i=1,10000
		T=Tenproduct(T2,(/3,4/),T3,(/3,4/))
	end do
	time2=omp_get_wtime()
	write(*,*)"No name",time2-time1
	
	write(*,*)"Test TensorName for product input single index"
	call cleanTensorName(T2)
	call cleanTensorName(T3)
		call setTensorName(T2,'A')
	call setTensorName(T3,'B')
	time1=omp_get_wtime()
	do i=1,10000
		T=Tenproduct(T2,'A.'+3,T3,'B.'+3)
	end do
	time2=omp_get_wtime()
	write(*,*)"Char Name",time2-time1
	call cleanTensorName(T2)
	call cleanTensorName(T3)
	call setTensorName(T2,3,(/1/),(/1,1/))
	call setTensorName(T3,3,(/2/),(/1,1/))
	time1=omp_get_wtime()
	do i=1,10000
		T=Tenproduct(T2,(/1/),(/1,1/),T3,(/2/),(/1,1/))
	end do
	time2=omp_get_wtime()
	write(*,*)"Int name",time2-time1
	time1=omp_get_wtime()
	do i=1,10000
		T=Tenproduct(T2,3,T3,3)
	end do
	time2=omp_get_wtime()
	write(*,*)"No name",time2-time1
end
	
	
