!==================================================================================================================================
subroutine preamble
	
	use globvar
	use mpi
	
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

	use globvar
	use param
	
	implicit none
	integer:: i,j,k,d
	
	write (*,*)'================================= TIME SUMMARY ================================'
	write (*,*)'Average Total CPU time = ', cputime
	write (*,*)'Average Output time =    ', output_time

end subroutine print_summary

!==================================================================================================================================
subroutine print_update

	use globvar
	use param
	
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
	