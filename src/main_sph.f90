program SPH

	use globvar
	use param
	
	use summary_m

	implicit none
	
	!Printing preamble to screen
	call preamble
	
	! setting k parameter for kernel radius (r = kh)
	select case (skf)
		case (1)
			scale_k = 2d0
		case (2)
			scale_k = 3d0
		case (3)
			scale_k = 2d0
		case (4)
			scale_k = 2d0
	end select
	
	call input(.false.)
	call virt_part(.false.)
	
	write(*,'(A24,1x,I9,1x,A19)') 'Total simulation size of',ntotal,'physical particles.'
	write(*,'(A24,1x,I9,1x,A19)') '                        ',nvirt,'virtual particles.'
	
	call allocateGlobalArrays( )
	
	!Creat physical and virtual boundary particles
	call input(.true.)
	call virt_part(.true.)
	
	!Entering discretized time-integration loop
	call time_integration
	
	!Printing post-amble to terminal
	call time_print
	call print_summary
	
	!Cleaning up
	call deallocateGlobalArrays

end