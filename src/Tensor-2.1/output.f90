
module output
	implicit none
	integer,save,private::output_cpu_number=0
	logical,save,private::log_flag=.false.!if false,there is no log,create a new file
	CHARACTER*100,private::log_address="outputdata/log"
	integer,private::		log_address_unit=999
!****************************************************************************
	!	MPI parameter
	integer,save,private::TimeEvProID=0,TimeEvProNum=1,TimeEvIerr=1
	integer,private,parameter::IDmin=0
	interface writemess
		module procedure writemess_char
		module procedure writemess_real
	end interface	
contains
	subroutine set_cpu_info(TimeEvProID_,TimeEvProNum_,TimeEvIerr_)
		integer,intent(in)::TimeEvProID_,TimeEvProNum_,TimeEvIerr_
		TimeEvProID=TimeEvProID_
		TimeEvProNum=TimeEvProNum_
		TimeEvIerr=TimeEvIerr_
		log_flag=.false.
		return
	end subroutine
	subroutine set_log_address(address)
		CHARACTER(len=*),intent(in)::address
		log_address=address
		return
	end subroutine
	subroutine set_log_unit(logunit)
		integer,intent(in)::logunit
		log_address_unit=logunit
		return
	end subroutine
	subroutine set_output_cpu(cpu)
		integer,intent(in)::cpu
		output_cpu_number=cpu
		return
	end subroutine
	
	
	
	subroutine openlog()
		if(log_flag)then
			open(unit=log_address_unit,file=log_address,STATUS='old',POSITION='APPEND')
		else
			open(unit=log_address_unit,file=log_address,STATUS='REPLACE',POSITION='APPEND')
			log_flag=.true.
		end if
		return
	end subroutine
	subroutine closelog()
		close(unit=log_address_unit)
		return
	end subroutine

	subroutine printInfo(time1,i,step_i,Nstep)
		real*8,intent(in)::time1(0:1)
		integer,intent(in)::i,step_i,Nstep
		real*8::runningtime
		CHARACTER*100::w1,w2
		if(TimeEvProID.eq.output_cpu_number) then
			runningtime=time1(1)-time1(0)
			write(w1,*)i-step_i
			write(w2,*)runningtime
			w1="i,running time (s)"//trim(adjustl(w1))//","//trim(adjustl(w2))
			call writemess(trim(adjustl(w1)))
			runningtime=runningtime/real((i-step_i))*(Nstep-i)
			call outputmess(runningtime)
		end if
	end subroutine
	
	subroutine printtime(Ttype)
		integer,intent(in)::Ttype
		if(TimeEvProID.eq.output_cpu_number) then
			if(Ttype.eq.1) then
				call openlog()
				call writemess("Start time is")
				call system("date '+%D%n%c'")
				call system("date '+%D%n%c' >> outputdata/log")
				call closelog()
			else 
			  call system("date '+%D%n%c'")
			  call system("date '+%D%n%c' >> outputdata/log")
		 endif
		end if
	end subroutine
	
!ccccccccccccccc   outputmess ccccccccccccccc				
	subroutine outputmess(cputime)
		real*8,intent(in)::cputime
		integer::times,timem,timeh,timed,temp
		CHARACTER(10) :: cput
		Character(8) :: cpud 
		CHARACTER(5) :: cpuz  
		CHARACTER*50::w1,w2,w3
		if(TimeEvProID.eq.output_cpu_number) then
			if(cputime.lt.60) then
				times=cputime
				timem=0
				timeh=0
				timed=0
			else if( (cputime.ge.60).and.(cputime.lt.3600)) then
				timem=cputime/60
				times=cputime-timem*60
				timeh=0
				timed=0
			else if((cputime.ge.3600).and.(cputime.lt.86400)) then
				timeh=cputime/3600
				temp=cputime-timeh*3600
				timem=temp/60
				times=temp-timem*60
				timed=0
			else
				timed=cputime/86400
				temp=cputime-timed*86400
				timeh=temp/3600
				temp=temp-timeh*3600
				timem=temp/60
				times=temp-timem*60
			end if
			CALL DATE_AND_TIME(DATE=cpud,TIME=cput,ZONE=cpuz) 
			call writemess("now the time is :")
			write(w1,*)cpud
			write(w2,*)cput
			w3=trim(adjustl(w1))//" "//trim(adjustl(w2))
			call writemess(trim(adjustl(w3)))
			call system("date '+%D%n%c' ")
			call writemess("The time to finish the job will be")
			
			write(w1,*)timed
			w3=" "//trim(adjustl(w1))//"day,"
			write(w1,*)timeh
			w3=trim(adjustl(w3))//trim(adjustl(w1))//"hour,"
			write(w1,*)timem
			w3=trim(adjustl(w3))//trim(adjustl(w1))//"minute,"
			write(w1,*)times
			w3=trim(adjustl(w3))//trim(adjustl(w1))//"second."
			
			call writemess(trim(adjustl(w3)))
			call writemess("--------------------------------")
		end if
		return
	end subroutine	
	subroutine outputmess_totaltime(cputime)
		real*8,intent(in)::cputime
		integer::times,timem,timeh,timed,temp
		CHARACTER(10) :: cput
		Character(8) :: cpud 
		CHARACTER(5) :: cpuz  
		CHARACTER*50::w1,w2,w3
		if(TimeEvProID.eq.output_cpu_number) then
			if(cputime.lt.60) then
				times=cputime
				timem=0
				timeh=0
				timed=0
			else if( (cputime.ge.60).and.(cputime.lt.3600)) then
				timem=cputime/60
				times=cputime-timem*60
				timeh=0
				timed=0
			else if((cputime.ge.3600).and.(cputime.lt.86400)) then
				timeh=cputime/3600
				temp=cputime-timeh*3600
				timem=temp/60
				times=temp-timem*60
				timed=0
			else
				timed=cputime/86400
				temp=cputime-timed*86400
				timeh=temp/3600
				temp=temp-timeh*3600
				timem=temp/60
				times=temp-timem*60
			end if
			CALL DATE_AND_TIME(DATE=cpud,TIME=cput,ZONE=cpuz) 
			call writemess("*****************************************")
			call writemess("************* Finish ! ******************")
			call writemess("*****************************************")
			call writemess("now the time is :")
			write(w1,*)cpud
			write(w2,*)cput
			w3=trim(adjustl(w1))//" "//trim(adjustl(w2))
			call writemess(trim(adjustl(w3)))
			call system("date '+%D%n%c' ")
			
			call writemess("The time to run the job is")
			write(w1,*)timed
			w3=" "//trim(adjustl(w1))//"day,"
			write(w1,*)timeh
			w3=trim(adjustl(w3))//trim(adjustl(w1))//"hour,"
			write(w1,*)timem
			w3=trim(adjustl(w3))//trim(adjustl(w1))//"minute,"
			write(w1,*)times
			w3=trim(adjustl(w3))//trim(adjustl(w1))//"second."
			call writemess(trim(adjustl(w3)))
		end if
		return
	end subroutine	
	
	subroutine writemess_char(mess,cpu_number)
		CHARACTER(len=*),intent(in)::mess
		integer,optional,intent(in)::cpu_number
		if(present(cpu_number))then
			if(TimeEvproID.eq.cpu_number)then
				call openlog()
				write(log_address_unit,*) trim(adjustl(mess))
				call closelog()
				write(*,*)trim(adjustl( mess))
			end if
		else
			if(TimeEvproID.eq.output_cpu_number)then
				call openlog()
				write(log_address_unit,*) trim(adjustl(mess))
				call closelog()
				write(*,*)trim(adjustl( mess))
			end if
		end if
		return
	end subroutine
	subroutine writemess_real(mess,cpu_number)
		real*8,intent(in)::mess
		integer,optional,intent(in)::cpu_number
		if(present(cpu_number))then
			if(TimeEvproID.eq.cpu_number)then
				call openlog()
				write(log_address_unit,*) mess
				call closelog()
				write(*,*)mess
			end if
		else
			if(TimeEvproID.eq.output_cpu_number)then
				call openlog()
				write(log_address_unit,*) mess
				call closelog()
				write(*,*) mess
			end if
		end if
		return
	end subroutine
end module
