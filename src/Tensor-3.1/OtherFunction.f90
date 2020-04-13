module OtherFunction
	use Tensor_type
	use usefull_function
	implicit none



contains

	subroutine ALLREDUCE_Tensor(inTensor,outTensor,OP,ierr,MPIcommon)
		type(Tensor),target,intent(in)::inTensor
		type(Tensor),target,intent(inout)::outTensor
		integer,intent(in)::ierr
		integer,intent(in)::op
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),mpi_comm
		type(Tensor),pointer::p1,p2
		integer::length,classtype
		integer,pointer::idata(:),idata2(:)
		real*4,pointer::sdata(:),sdata2(:)
		real*8,pointer::ddata(:),ddata2(:)
		complex*8,pointer::cdata(:),cdata2(:)
		complex*16,pointer::zdata(:),zdata2(:)
		logical,pointer::ldata(:),ldata2(:)
		character(len=characterlen),pointer::adata(:),adata2(:)
		logical::ALLgoonFlag,goonFlag

		tag=1
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		p1=>inTensor
		p2=>outTensor
		if(associated(p1,p2))then
			call writemess('ERROR in ALLREDUCE_Tensor,input Tensor and output Tensor can not be the same one',-1)
			call error_stop
		end if
		p1=>null()
		p2=>null()
			! if the Tensor is empty
		goonFlag=inTensor%getFlag()
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in ALLREDUCE_Tensor,the is no date in one or some Tensors',-1)
			call error_stop
		end if
			! if the Tensor is the same data type
		classtype=inTensor%getType()
		call MPI_BCAST(classtype,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getType().ne.classtype)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in ALLREDUCE_Tensor,the Data type in the Tensors are not the sames',-1)
			call error_stop
		end if
			! if the length of the Tensor is the same
		length=inTensor%getTotalData()
		call MPI_BCAST(length,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getTotalData().ne.length)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in ALLREDUCE_Tensor,the length od the Tensor is not the same',-1)
			call error_stop
		end if
		call outTensor%empty()
		call outTensor%allocate(inTensor)

		select case(inTensor%getType())
			case(1)
				call inTensor%pointer(idata)
				call outTensor%pointer(idata2)
				call MPI_ALLREDUCE(idata,idata2,length,MPI_INTEGER,OP,mpi_comm,ierr)
			case(2)
				call inTensor%pointer(sdata)
				call outTensor%pointer(sdata2)
				call MPI_ALLREDUCE(sdata,sdata2,length,MPI_real,OP,mpi_comm,ierr)
			case(3)
				call inTensor%pointer(ddata)
				call outTensor%pointer(ddata2)
				call MPI_ALLREDUCE(ddata,ddata2,length,MPI_double_precision,OP,mpi_comm,ierr)
			case(4)
				call inTensor%pointer(cdata)
				call outTensor%pointer(cdata2)
				call MPI_ALLREDUCE(cdata,cdata2,length,MPI_complex,OP,mpi_comm,ierr)
			case(5)
				call inTensor%pointer(zdata)
				call outTensor%pointer(zdata2)
				call MPI_ALLREDUCE(zdata,zdata2,length,MPI_double_complex,OP,mpi_comm,ierr)
			case(6)
				call inTensor%pointer(ldata)
				call outTensor%pointer(ldata2)
				call MPI_ALLREDUCE(ldata,ldata2,length,MPI_logical,OP,mpi_comm,ierr)
			case(7)
				call inTensor%pointer(adata)
				call outTensor%pointer(adata2)
				call MPI_ALLREDUCE(adata,adata2,length,MPI_character,OP,mpi_comm,ierr)
		end  select
		return
	end subroutine


	subroutine REDUCE_Tensor(inTensor,outTensor,OP,root,ierr,MPIcommon)
		type(Tensor),target,intent(in)::inTensor
		type(Tensor),target,intent(inout)::outTensor
		integer,intent(in)::ierr,root
		integer,intent(in)::op
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),mpi_comm
		type(Tensor),pointer::p1,p2
		integer::length,classtype
		integer,pointer::idata(:),idata2(:)
		real*4,pointer::sdata(:),sdata2(:)
		real*8,pointer::ddata(:),ddata2(:)
		complex*8,pointer::cdata(:),cdata2(:)
		complex*16,pointer::zdata(:),zdata2(:)
		logical,pointer::ldata(:),ldata2(:)
		character(len=characterlen),pointer::adata(:),adata2(:)
		logical::ALLgoonFlag,goonFlag

		tag=1
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		p1=>inTensor
		p2=>outTensor
		if(associated(p1,p2))then
			call writemess('ERROR in REDUCE_Tensor,input Tensor and output Tensor can not be the same one',-1)
			call error_stop
		end if
		p1=>null()
		p2=>null()
			! if the Tensor is empty
		goonFlag=inTensor%getFlag()
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in REDUCE_Tensor,the is no date in one or some Tensors',-1)
			call error_stop
		end if
			! if the Tensor is the same data type
		classtype=inTensor%getType()
		call MPI_BCAST(classtype,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getType().ne.classtype)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in REDUCE_Tensor,the Data type in the Tensors are not the sames',-1)
			call error_stop
		end if
			! if the length of the Tensor is the same
		length=inTensor%getTotalData()
		call MPI_BCAST(length,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getTotalData().ne.length)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in REDUCE_Tensor,the length od the Tensor is not the same',-1)
			call error_stop
		end if
		call outTensor%empty()
		call outTensor%allocate(inTensor)

		select case(inTensor%getType())
			case(1)
				call inTensor%pointer(idata)
				call outTensor%pointer(idata2)
				call MPI_REDUCE(idata,idata2,length,MPI_INTEGER,OP,root,mpi_comm,ierr)
			case(2)
				call inTensor%pointer(sdata)
				call outTensor%pointer(sdata2)
				call MPI_REDUCE(sdata,sdata2,length,MPI_real,OP,root,mpi_comm,ierr)
			case(3)
				call inTensor%pointer(ddata)
				call outTensor%pointer(ddata2)
				call MPI_REDUCE(ddata,ddata2,length,MPI_double_precision,OP,root,mpi_comm,ierr)
			case(4)
				call inTensor%pointer(cdata)
				call outTensor%pointer(cdata2)
				call MPI_REDUCE(cdata,cdata2,length,MPI_complex,OP,root,mpi_comm,ierr)
			case(5)
				call inTensor%pointer(zdata)
				call outTensor%pointer(zdata2)
				call MPI_REDUCE(zdata,zdata2,length,MPI_double_complex,OP,root,mpi_comm,ierr)
			case(6)
				call inTensor%pointer(ldata)
				call outTensor%pointer(ldata2)
				call MPI_REDUCE(ldata,ldata2,length,MPI_logical,OP,root,mpi_comm,ierr)
			case(7)
				call inTensor%pointer(adata)
				call outTensor%pointer(adata2)
				call MPI_REDUCE(adata,adata2,length,MPI_character,OP,root,mpi_comm,ierr)
		end  select
		return
	end subroutine

end module