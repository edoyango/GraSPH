program SPH

	use globvar
	use param
	
	use input_m
	use kernel_m
	use summary_m
	use time_integration_m

	implicit none
	
	!Printing preamble to screen
	call preamble
	
	! setting k parameter for kernel radius (r = kh)
	scale_k = kernel_k(skf)
	
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