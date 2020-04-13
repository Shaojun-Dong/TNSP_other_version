	subroutine ExternalcheckSymmetryRule(Rule1,Rule2,lenname1,lenname2)
		use Tools
		integer,intent(in)::Rule1,Rule2
		character(len=*),intent(in)::lenname1,lenname2
		call writemess('You have not written the external fucntion for symmetry tensors yet')
		call error_stop
	end subroutine


	subroutine ExternaldefaultreverseSymmetryRule(T,dimen,LD1,LD2)
		use Tools
		use Tensor_type
		use SymDimension_typede
		integer,intent(in)::LD1,LD2
		Type(Tensor),intent(inout)::T(LD1,LD2)
		Type(Symdimension),intent(inout)::dimen
		call writemess('You have not written the external fucntion for symmetry tensors yet')
		call error_stop
	end subroutine

	type(Tensor) function ExternalNonZeroElement(dimen)
		use Tools
		use SymDimension_typede
		use Tensor_type
		type(SymDimension),intent(in)::dimen
		call writemess('You have not written the external fucntion for symmetry tensors yet')
		call error_stop	
	end function

	logical function ExternalRule(dimen,indices,length)
		use Tools
		use SymDimension_typede
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::length
		integer,intent(in)::indices(length)
		call writemess('You have not written the external fucntion for symmetry tensors yet')
		call error_stop	
		ExternalRule=.false.
	end function

	type(QuanNum) function ExternalfuseOrder(order,Q1,Q2,newRule_)
		use Tools
		use QuantumNumber_Type
		use Tensor_type
		type(QuanNum),intent(in)::Q1,Q2
		type(Tensor),intent(inout)::order
		integer,intent(in)::newRule_
		call writemess('You have not written the external fucntion for symmetry tensors yet')
		call error_stop
	end function

	logical function ExternalifParity(dimen,vec,ith,jth,rank)
		use Tools
		use SymDimension_typede
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::rank
		integer,intent(in)::vec(rank)
		integer,intent(in)::ith,jth
		call writemess('You have not written the external fucntion for symmetry tensors yet')
		call error_stop
		ExternalifParity=.false.
	end function

	integer function ExternalQaunNumParity(dimen,ith,jth)
		use Tools
		use SymDimension_typede
		type(SymDimension),intent(in)::dimen
		integer,intent(in)::ith,jth
		call writemess('You have not written the external fucntion for symmetry tensors yet')
		call error_stop
		ExternalQaunNumParity=0
	end function
