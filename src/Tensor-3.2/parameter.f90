module parameter_type! Do not use pointer
	use Tensor_type
	use usefull_function
	implicit none
	private
	public::List
	type List
		type(Tensor)::parameter
		type(Tensor)::name
		integer::length=0
	contains
			procedure,public::deallocate=>clean 
			procedure,public::getLength
			procedure,public::getClassType
			procedure,public::allocate=>allocateParameter
			generic,public::print=>printParameter1,printParameter2,printParameter3
			procedure::printParameter1,printParameter2,printParameter3
			procedure::isetValue,ssetValue,dsetValue,csetValue,zsetValue,lsetValue,asetValue
			procedure::isetValue2,ssetValue2,dsetValue2,csetValue2,zsetValue2,lsetValue2,asetValue2
			generic,public::setValue=>isetValue,ssetValue,dsetValue,csetValue,zsetValue,lsetValue,asetValue,&
				isetValue2,ssetValue2,dsetValue2,csetValue2,zsetValue2,lsetValue2,asetValue2
			procedure::initialNameAll,initialNamei
			procedure,public::index=>elementindex
			generic,public::setName=>initialNameAll,initialNamei
			generic,public::write=>writeout1,writeout2
			procedure::writeout1,writeout2
			procedure,public::read=>readout
			procedure::ielement1,selement1,delement1,celement1,zelement1,lelement1,aelement1,Telement1
			procedure::ielement2,selement2,delement2,celement2,zelement2,lelement2,aelement2,Telement2
			generic,public::ii=>ielement1,ielement2
			generic,public::si=>selement1,selement2 
			generic,public::di=>delement1,delement2
			generic,public::ci=>celement1,celement2
			generic,public::zi=>zelement1,zelement2
			generic,public::li=>lelement1,lelement2
			generic,public::ai=>aelement1,aelement2
			generic,public::i=>Telement1,Telement2
			procedure,public::check=>check_paramters
			generic,public::subList=>subList1,subList2
			procedure::subList1,subList2
	end type
	
	public::writemess
	interface writemess
		module procedure writeoparametemess
	end interface
	
	
	public::operator(+)
	interface operator(+)
		module procedure AddParameter
	end interface
	
	
	public::assignment(=)
	interface assignment(=)
		module procedure copyParapemter
		module procedure initial
	end interface
	
	public::MPI_BCAST_parameter,MPI_Send_parameter
contains

	subroutine clean(p)
		class(List),intent(inout)::p
		call p%parameter%deallocate()
		call p%name%deallocate()
		p%length=0
		return
	end subroutine
	
	character(len=characterlen) function getClassType(L)
		class(List),intent(in)::L
		getClassType=L%Parameter%getClassType()
		return
	end function
	
	integer function getLength(L)
		class(List),intent(in)::L
		getLength=L%length
		return
	end function
	
	subroutine copyParapemter(pinout,inp)
		type(List),intent(inout)::pinout
		type(List),intent(in)::inp
		if(inp%length.eq.0)then
			call pinout%deallocate()
			return
		end if
		pinout%parameter=inp%parameter
		pinout%name=inp%name
		pinout%length=inp%length
		return
	end subroutine
	subroutine initial(pinout,arrayTensor)
		type(List),intent(inout)::pinout
		type(Tensor),intent(in)::arrayTensor(:)
		if(size(arrayTensor).ne.2)then
			call writemess('ERROR in (=) for type(List)',-1)
			call error_stop
		end if
		if(arrayTensor(1)%getType().ne.7)then
			call writemess('ERROR in (=) for type(List), the name should be character',-1)
			call error_stop
		end if
		pinout%length=arrayTensor(1)%getTotalData()
		if(pinout%length.ne.arrayTensor(2)%getTotalData())then
			call writemess('ERROR in (=) for type(List), ERROR LENGTH',-1)
			call error_stop
		end if
		pinout%parameter=arrayTensor(2)
		pinout%name=arrayTensor(1)
		return
	end subroutine
	
	
	subroutine allocateParameter(p,DataTpye,length)
		class(List),intent(inout)::p
		character(len=*),intent(in)::DataTpye
		integer,intent(in)::length
		call p%parameter%allocate([length],DataTpye)
		call p%name%allocate([length],'character')
		p%length=length
		return
	end subroutine
	
	subroutine isetValue(p,ith,val)
		class(List),intent(inout)::p
		integer,intent(in)::val
		integer,intent(in)::ith
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	subroutine ssetValue(p,ith,val)
		class(List),intent(inout)::p
		real*4,intent(in)::val
		integer,intent(in)::ith
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	subroutine dsetValue(p,ith,val)
		class(List),intent(inout)::p
		real*8,intent(in)::val
		integer,intent(in)::ith
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	subroutine csetValue(p,ith,val)
		class(List),intent(inout)::p
		complex*8,intent(in)::val
		integer,intent(in)::ith
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	subroutine zsetValue(p,ith,val)
		class(List),intent(inout)::p
		complex*16,intent(in)::val
		integer,intent(in)::ith
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	subroutine lsetValue(p,ith,val)
		class(List),intent(inout)::p
		logical,intent(in)::val
		integer,intent(in)::ith
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	subroutine asetValue(p,ith,val)
		class(List),intent(inout)::p
		character(len=*),intent(in)::val
		integer,intent(in)::ith
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	
	
	subroutine isetValue2(p,namei,val)
		class(List),intent(inout)::p
		integer,intent(in)::val
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	
	subroutine ssetValue2(p,namei,val)
		class(List),intent(inout)::p
		real*4,intent(in)::val
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	
	subroutine dsetValue2(p,namei,val)
		class(List),intent(inout)::p
		real*8,intent(in)::val
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	
	subroutine csetValue2(p,namei,val)
		class(List),intent(inout)::p
		complex*8,intent(in)::val
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	
	subroutine zsetValue2(p,namei,val)
		class(List),intent(inout)::p
		complex*16,intent(in)::val
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	
	subroutine lsetValue2(p,namei,val)
		class(List),intent(inout)::p
		logical,intent(in)::val
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	
	subroutine asetValue2(p,namei,val)
		class(List),intent(inout)::p
		character(len=*),intent(in)::val
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		call p%parameter%setValue(ith,val)
		return
	end subroutine
	
	subroutine initialNameAll(p,Allnamei)
		class(List),intent(inout)::p
		character(len=*),intent(in)::Allnamei(:)
		if(p%length.eq.0)then
			call writemess('ERROR Name length in initial parameter,Do not allocate parameter yet',-1)
			call error_stop
		end if
		if(size(Allnamei).ne.p%length)then
			call writemess('ERROR Name length in initial parameter',-1)
			call error_stop
		end if
		p%name=Allnamei
		return
	end subroutine
	subroutine initialNamei(p,ith,Allnamei)
		class(List),intent(inout)::p
		integer,intent(in)::ith
		character(len=*),intent(in)::Allnamei
		if(ith.gt.p%length)then
			call writemess('ERROR Name length in initial parameter',-1)
			call error_stop
		end if
		call p%name%setValue(ith,Allnamei)
		return
	end subroutine
	
	
	integer function ielement1(p,ith)result(Res)
		class(List),intent(in)::p
		integer,intent(in)::ith
		Res=p%parameter%ii(ith)
		return
	end function 
	integer function ielement2(p,namei)result(Res)
		class(List),intent(in)::p
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		Res=p%parameter%ii(ith)
		return
	end function 
	real*4 function selement1(p,ith)result(Res)
		class(List),intent(in)::p
		integer,intent(in)::ith
		Res=p%parameter%si(ith)
		return
	end function 
	real*4 function selement2(p,namei)result(Res)
		class(List),intent(in)::p
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		Res=p%parameter%si(ith)
		return
	end function 
	real*8 function delement1(p,ith)result(Res)
		class(List),intent(in)::p
		integer,intent(in)::ith
		Res=p%parameter%di(ith)
		return
	end function 
	real*8 function delement2(p,namei)result(Res)
		class(List),intent(in)::p
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		Res=p%parameter%di(ith)
		return
	end function 
	complex*8 function celement1(p,ith)result(Res)
		class(List),intent(in)::p
		integer,intent(in)::ith
		Res=p%parameter%ci(ith)
		return
	end function 
	complex*8 function celement2(p,namei)result(Res)
		class(List),intent(in)::p
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		Res=p%parameter%ci(ith)
		return
	end function 
	complex*16 function zelement1(p,ith)result(Res)
		class(List),intent(in)::p
		integer,intent(in)::ith
		Res=p%parameter%zi(ith)
		return
	end function 
	complex*16 function zelement2(p,namei)result(Res)
		class(List),intent(in)::p
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		Res=p%parameter%zi(ith)
		return
	end function 
	logical function lelement1(p,ith)result(Res)
		class(List),intent(in)::p
		integer,intent(in)::ith
		Res=p%parameter%li(ith)
		return
	end function 
	logical function lelement2(p,namei)result(Res)
		class(List),intent(in)::p
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		Res=p%parameter%li(ith)
		return
	end function 
	character(len=max_len_of_char_in_TData) function aelement1(p,ith)result(Res)
		class(List),intent(in)::p
		integer,intent(in)::ith
		Res=p%parameter%ai(ith)
		return
	end function 
	character(len=max_len_of_char_in_TData) function aelement2(p,namei)result(Res)
		class(List),intent(in)::p
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		Res=p%parameter%ai(ith)
		return
	end function 
	
	type(Tensor) function Telement1(p,ith)result(Res)
		class(List),intent(in)::p
		integer,intent(in)::ith
		Res=p%parameter%i(ith)
		return
	end function 
	type(Tensor) function Telement2(p,namei)result(Res)
		class(List),intent(in)::p
		character(len=*),intent(in)::namei
		integer::ith
		ith=p%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		Res=p%parameter%i(ith)
		return
	end function 
	
	
	
	subroutine writeout1(p,w,uni)
		class(List),intent(in)::p
		character(len=*),intent(in)::w
		integer,intent(in)::uni
		integer::i,ptype
		ptype=p%parameter%getType()
		write(uni,*)'**********************************************************'
		write(uni,*)'***length_and_type',p%length,p%parameter%getClassType()
		write(uni,*)trim('***'+w)
		write(uni,*)'**********************************************************'
		select case(ptype)
			case(1)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),p%parameter%ii(i)
				end do
			case(2)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),p%parameter%si(i)
				end do
			case(3)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),p%parameter%di(i)
				end do
			case(4)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),real(p%parameter%ci(i),kind=4),aimag(p%parameter%ci(i))
				end do
			case(5)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),real(p%parameter%zi(i),kind=8),dimag(p%parameter%zi(i))
				end do
			case(6)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),p%parameter%li(i)
				end do
			case(7)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),"    ",p%parameter%ai(i)
				end do
		end select
	end subroutine
	
	subroutine writeout2(p,uni)
		class(List),intent(in)::p
		integer,intent(in)::uni
		integer::i,ptype
		ptype=p%parameter%getType()
		write(uni,*)'**********************************************************'
		write(uni,*)'***length_and_type',p%length,p%parameter%getClassType()
		write(uni,*)'***     '
		write(uni,*)'**********************************************************'
		select case(ptype)
			case(1)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),p%parameter%ii(i)
				end do
			case(2)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),p%parameter%si(i)
				end do
			case(3)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),p%parameter%di(i)
				end do
			case(4)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),real(p%parameter%ci(i),kind=4),aimag(p%parameter%ci(i))
				end do
			case(5)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),real(p%parameter%zi(i),kind=8),dimag(p%parameter%zi(i))
				end do
			case(6)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),p%parameter%li(i)
				end do
			case(7)
				do i=1,p%length
					write(uni,*)"    ",trim(p%name%ai(i)),"    ",p%parameter%ai(i)
				end do
		end select
	end subroutine
	
	subroutine readout(p,uni)
		class(List),intent(inout)::p
		integer,intent(in)::uni
		integer::i,ptype,length
		CHARACTER(len=50)::notused,classtype
		integer,allocatable::idata(:)
		real*4,allocatable::sdata(:),sdata2(:)
		real*8,allocatable::ddata(:),ddata2(:)
		character(len=max_len_of_char_in_TData),allocatable::adata(:)
		logical,allocatable::ldata(:)
		character(len=max_len_of_char_in_TData),allocatable::nameadata(:)
		type(Tensor)::temp
		read(uni,*)notused
		read(uni,*)notused,length,classtype
		if(length.eq.0)then
			read(uni,*)notused
			read(uni,*)notused
			call p%deallocate()
			return
		end if
		read(uni,*)notused
		read(uni,*)notused
		call temp%setType(classtype)
		ptype=temp%getType()
		allocate(nameadata(length))
		select case(ptype)
			case(1)
				allocate(idata(length))
				do i=1,length
					read(uni,*)nameadata(i),idata(i)
				end do
				p%parameter=idata
				p%name=nameadata
				p%length=length
			case(2)
				allocate(sdata(length))
				do i=1,length
					read(uni,*)nameadata(i),sdata(i)
				end do
				p%parameter=sdata
				p%name=nameadata
				p%length=length
			case(3)
				allocate(ddata(length))
				do i=1,length
					read(uni,*)nameadata(i),ddata(i)
				end do
				p%parameter=ddata
				p%name=nameadata
				p%length=length
			case(4)
				allocate(sdata(length))
				allocate(sdata2(length))
				do i=1,length
					read(uni,*)nameadata(i),sdata(i),sdata2
				end do
				p%parameter=cmplx(sdata,sdata2,kind=4)
				p%name=nameadata
				p%length=length
			case(5)
				allocate(ddata(length))
				allocate(ddata2(length))
				do i=1,length
					read(uni,*)nameadata(i),ddata(i),ddata2
				end do
				p%parameter=dcmplx(ddata,ddata2)
				p%name=nameadata
				p%length=length
			case(6)
				allocate(ldata(length))
				do i=1,length
					read(uni,*)nameadata(i),ldata(i)
				end do
				p%parameter=ldata
				p%name=nameadata
				p%length=length
			case(7)
				allocate(adata(length))
				do i=1,length
					read(uni,*)nameadata(i),adata(i)
				end do
				p%parameter=adata
				p%name=nameadata
				p%length=length
		end select
	end subroutine
	
	subroutine writeoparametemess(p,cpuInfo)
		class(List),intent(in)::p
		integer,intent(in),optional::cpuInfo
		integer::i,ptype
		ptype=p%parameter%getType()
		select case(ptype)
			case(1)
				do i=1,p%length
					call writemess(p%name%ai(i)+'='+p%parameter%ii(i),cpuInfo )
				end do
			case(2)
				do i=1,p%length
					call writemess(p%name%ai(i)+'='+p%parameter%si(i),cpuInfo )
				end do
			case(3)
				do i=1,p%length
					call writemess(p%name%ai(i)+'='+p%parameter%di(i),cpuInfo )
				end do
			case(4)
				do i=1,p%length
					call writemess(p%name%ai(i)+'='+p%parameter%ci(i),cpuInfo )
				end do
			case(5)
				do i=1,p%length
					call writemess(p%name%ai(i)+'='+p%parameter%zi(i),cpuInfo )
				end do
			case(6)
				do i=1,p%length
					call writemess(p%name%ai(i)+'='+p%parameter%li(i),cpuInfo )
				end do
			case(7)
				do i=1,p%length
					call writemess(p%name%ai(i)+'='+p%parameter%ai(i),cpuInfo )
				end do
		end select
	end subroutine
	
	
		
	subroutine printParameter1(p,lenA,uni)
		class(List),intent(in)::p
		integer,intent(in)::lenA
		integer,intent(in),optional::uni
		character(len=20)::F
		integer::i,ptype
		ptype=p%parameter%getType()
		F='A'+lenA+''
		F='('+F+',A2,'+F+')'
		if(present(uni))then
			select case(ptype)
				case(1)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%ii(i)
					end do
				case(2)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%si(i)
					end do
				case(3)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%di(i)
					end do
				case(4)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%ci(i)
					end do
				case(5)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%zi(i)
					end do
				case(6)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%li(i)
					end do
				case(7)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%ai(i)
					end do
			end select
			return
		end if
		select case(ptype)
			case(1)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%ii(i)
				end do
			case(2)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%si(i)
				end do
			case(3)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%di(i)
				end do
			case(4)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%ci(i)
				end do
			case(5)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%zi(i)
				end do
			case(6)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%li(i)
				end do
			case(7)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%ai(i)
				end do
		end select
		return
	end subroutine
		
	subroutine printParameter2(p,lenA,uni)
		class(List),intent(in)::p
		integer,intent(in)::lenA(2)
		integer,intent(in),optional::uni
		character(len=20)::F
		integer::i,ptype
		ptype=p%parameter%getType()
		F='('+'A'+lenA(1)+',A2,'+'A'+lenA(2)+')'
		if(present(uni))then
			select case(ptype)
				case(1)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%ii(i)
					end do
				case(2)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%si(i)
					end do
				case(3)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%di(i)
					end do
				case(4)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%ci(i)
					end do
				case(5)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%zi(i)
					end do
				case(6)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%li(i)
					end do
				case(7)
					do i=1,p%length
						write(uni,F)' '+p%name%ai(i),'=',' '+p%parameter%ai(i)
					end do
			end select
			return
		end if
		select case(ptype)
			case(1)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%ii(i)
				end do
			case(2)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%si(i)
				end do
			case(3)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%di(i)
				end do
			case(4)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%ci(i)
				end do
			case(5)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%zi(i)
				end do
			case(6)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%li(i)
				end do
			case(7)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%ai(i)
				end do
		end select
		return
	end subroutine
	subroutine printParameter3(p)
		class(List),intent(in)::p
		integer::lenA(2)
		character(len=20)::F
		integer::i,ptype
		lenA=[15,20]
		ptype=p%parameter%getType()
		F='('+'A'+lenA(1)+',A2,'+'A'+lenA(2)+')'
		select case(ptype)
			case(1)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%ii(i)
				end do
			case(2)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%si(i)
				end do
			case(3)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%di(i)
				end do
			case(4)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%ci(i)
				end do
			case(5)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%zi(i)
				end do
			case(6)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%li(i)
				end do
			case(7)
				do i=1,p%length
					write(*,F)' '+p%name%ai(i),'=',' '+p%parameter%ai(i)
				end do
		end select
		return
	end subroutine
	
	integer function elementindex(p,ith)
		class(List),intent(in)::p
		character(len=*),intent(in)::ith
		elementindex=p%name%which(ith)
		return
	end function
		
		
	
	subroutine check_paramters(p2,p1)
		class(List),intent(in)::P2
		type(List),intent(in)::p1
		integer::i,ith,ptype
		character(len=max_len_of_char_in_TData)::namei
		if(p1%length.ne.p2%length)then
			call writemess('The paramters have diferent length')
			return
		end if
		ptype=p1%parameter%getType()
		if(ptype.ne.p2%parameter%getType())then
			call writemess('The paramters have diferent data type')
			return
		end if
		do i=1,p1%length
			namei=p1%name%ai(i)
			ith=p2%index(namei)
			if(ith.eq.0)then
				call writemess('Can not find'+(' '+namei)+'in the secend parameter')
			else	
				select case(ptype)
					case(1)
						if(p1%ii(i).ne.p2%ii(ith))then
							call writemess('The parameter change:'+(' '+namei)+' ='+p1%ii(i)+'-->'+p2%ii(ith))
						end if
					case(2)
						if(p1%si(i).ne.p2%si(ith))then
							call writemess('The parameter change:'+(' '+namei)+' ='+p1%si(i)+'-->'+p2%si(ith))
						end if
					case(3)
						if(p1%di(i).ne.p2%di(ith))then
							call writemess('The parameter change:'+(' '+namei)+' ='+p1%di(i)+'-->'+p2%di(ith))
						end if
					case(4)
						if(p1%ci(i).ne.p2%ci(ith))then
							call writemess('The parameter change:'+(' '+namei)+' ='+p1%ci(i)+'-->'+p2%ci(ith))
						end if
					case(5)
						if(p1%zi(i).ne.p2%zi(ith))then
							call writemess('The parameter change:'+(' '+namei)+' ='+p1%zi(i)+'-->'+p2%zi(ith))
						end if
					case(6)
						if(p1%li(i).neqv.p2%li(ith))then
							call writemess('The parameter change:'+(' '+namei)+' ='+p1%li(i)+'-->'+p2%li(ith))
						end if
					case(7)
						if(p1%ai(i).nequ.p2%ai(ith))then
							call writemess('The parameter change:'+(' '+namei)+' ='+p1%ai(i)+'-->'+p2%ai(ith))
						end if
				end select
			end if
		end do
	end subroutine
	
	type(List) function AddParameter(List1,List2)
		type(List),intent(in)::List1,List2
		integer::length,Classtype
		integer::total1,total2
		length=List1%length+List2%length
		Classtype=List1%parameter%getType()
		if(Classtype.ne.List2%parameter%getType())then
			call writemess('Can not connect the List')
			call error_stop
		end if
		call AddParameter%allocate(List2%parameter%getClassType(),length )
		total1=List1%parameter%getTotalData()
		total2=List2%parameter%getTotalData()
		call AddParameter%Parameter%setValue([1,total1],List1%parameter,[1,total1])
		call AddParameter%name%setValue([1,total1],List1%name,[1,total1])
		call AddParameter%Parameter%setValue([total1+1,length],List2%parameter,[1,total2])
		call AddParameter%name%setValue([total1+1,length],List2%name,[1,total2])
		return
	end function
	
	type(List) function subList1(L,ith,jth)result(Res)
		class(List),intent(in)::L
		integer,intent(in)::ith,jth
		integer::subLength
		subLength=jth-ith+1
		if(subLength.le.0)then
			call writemess('ERROR in subList, length<0,ith='+ith+',jth='+jth)
			call error_stop
		end if
		call Res%allocate(L%getClassType(),subLength)
		call Res%parameter%setValue([1,subLength],L%parameter,[ith,jth])
		call Res%name%setValue([1,subLength],L%name,[ith,jth])
		return
	end function
	
	type(List) function subList2(L,namei,namej)result(Res)
		class(List),intent(inout)::L
		character(len=*),intent(in)::namei,namej
		integer::ith,jth
		integer::subLength
		ith=L%name%which(namei)
		if(ith.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namei))
			call error_stop
		end if
		jth=L%name%which(namej)
		if(jth.eq.0)then
			call writemess('Do not Find the name in parameter',-1)
			call writemess('The Name is'+(' '+namej))
			call error_stop
		end if
		if(ith.gt.ith)then
			call writemess('ERROR in subList, length<0,ith='+ith+',jth='+jth)
			call writemess('namei='+namei+'namej='+namej)
			call writemess('The data in list are')
			call writemess(L)
			call error_stop
		end if
		subLength=jth-ith+1
		call Res%allocate(L%getClassType(),subLength)
		call Res%parameter%setValue([1,subLength],L%parameter,[ith,jth])
		call Res%name%setValue([1,subLength],L%name,[ith,jth])
		return
	end function 
		
		
	subroutine MPI_BCAST_parameter(P,ID,ierr,MPIcommon)
		type(List),intent(inout)::p
		integer,intent(in)::ID,ierr
		integer,optional,intent(in)::MPIcommon
		call MPI_BCAST_Tensor(p%parameter,ID,ierr,MPIcommon)
		call MPI_BCAST_Tensor(p%name,ID,ierr,MPIcommon)
		p%length=p%parameter%getTotalData()
		return
	end subroutine
	
	subroutine MPI_send_parameter(p1,p2,ID1,ID2,ierr,MPIcommon)
		type(List),intent(in)::p1
		type(List),intent(inout)::p2
		integer,intent(in)::ID1,ID2,ierr
		integer,optional,intent(in)::MPIcommon
		call MPI_send_Tensor(p1%parameter,p2%parameter,ID1,ID2,ierr,MPIcommon)
		call MPI_send_Tensor(p1%name,p2%name,ID1,ID2,ierr,MPIcommon)
		p2%length=p2%parameter%getTotalData()
		return
	end subroutine


end module
