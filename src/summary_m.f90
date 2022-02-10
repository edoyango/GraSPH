module summary_m

	use globvar, only: ntotal,nvirt,niac,pairs,parts,maxtimestep,print_step,save_step,cputime,output_time,itimestep,time
	
	public:: time_print,preamble,print_summary,print_update
	
contains
	
	!==============================================================================================================================
	subroutine time_print
	!***************************************************************************************************
	!   TIME_PRINT:      Print out the current date and time.
	!   The standard Fortran 90 routine DATE_AND_TIME is used to get the current
	!   date and time strings.
	!***************************************************************************************************
		implicit none
		integer, parameter :: output = 6
		
		!Local scalars.
		character ( len =  8 ) :: datstr
		character ( len = 10 ) :: timstr
		
		!Get the current date and time.
		call date_and_time ( datstr, timstr )
		
		!Write out the date and time.
		write ( output, "(/A)"  ) "                  Date = " // datstr(7:8) // "/" // &
											datstr(5:6) // "/" // &
											datstr(1:4)
		write ( output, "(A)"   ) "                  Time = " // timstr(1:2) // ":" // &
											timstr(3:4) // ":" // &
											timstr(5:10)
		write ( output, *)
		
	end subroutine time_print

	!==================================================================================================================================
	subroutine preamble
		
		implicit none
		character(len=100):: args(3)
		integer:: ierr
		
		if (command_argument_count().lt.3) then
			write(*,'(A92)') '!!!ERROR!!! --- 3 commandline arguments must be provided: maxtimestep, print_step, save_step'
			write(*,'(A33)') '                Program ending...'
			stop
		else
			call get_command_argument(1,args(1))
			call get_command_argument(2,args(2))
			call get_command_argument(3,args(3))
			read(args,*) maxtimestep,print_step,save_step
		end if
		
		call time_print
		
		write(*,'(A36)') 'Serial code - running on one core...'
		
		write(*,'(A8,I7,A9)') 'Running ',maxtimestep,' step(s).' 
		write(*,'(A33,I7,A9)') 'Printing summary to screen every ',print_step,' step(s).'
		write(*,'(A29,I7,A9)') 'Writing output to disc every ',save_step,' step(s).'
		
	end subroutine preamble
	
	!==================================================================================================================================
	subroutine print_summary
	
		implicit none
		integer:: i,j,k,d
		
		write (*,*)'================================= TIME SUMMARY ================================'
		write (*,*)'Average CPU compute time = ', cputime
		write (*,*)'Average Output time =      ', output_time
	
	end subroutine print_summary
	
	!==================================================================================================================================
	subroutine print_update
	
		implicit none
		integer:: ns(ntotal+nvirt)
		integer:: maxp,minp,sumiac,maxiac,miniac,noiac,i,j,k,d
		
		write(*,*)'_______________________________________________________________________________'
		write(*,*)'  current number of time step =', itimestep,'     current time=', real(time)
		write(*,*)'                                                 Walltime    =', real(cputime)
		write(*,*)'_______________________________________________________________________________'
		
		!Statistics for the interaction, Print information to screen
		ns(:) = 0
		do k = 1,niac
			ns(pairs(k)%i%ind) = ns(pairs(k)%i%ind) + 1
			ns(pairs(k)%j%ind) = ns(pairs(k)%j%ind) + 1
		end do
		
		sumiac = 0
		maxiac = -huge(1)
		noiac  = 0
		miniac = huge(1)
		maxp = 1
		do i=1,ntotal+nvirt
			sumiac = sumiac + ns(parts(i)%ind)
			if (ns(parts(i)%ind) > maxiac) then
				maxiac = ns(parts(i)%ind)
				maxp = i
			endif
			if (ns(parts(i)%ind) < miniac) then 
				miniac = ns(parts(i)%ind)
				minp = i
			endif
			if (ns(parts(i)%ind).eq.0) noiac = noiac + 1
		enddo
		
		print *,' >> Statistics: interactions per particle:'
		print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
		print *,'**** Particle:',minp, ' minimal interactions:',miniac
		print *,'**** Average :',real(sumiac)/real(ntotal)
		print *,'**** Total pairs : ',niac
		print *,'**** Particles with no interactions:',noiac
		
	end subroutine print_update
		
end module summary_m