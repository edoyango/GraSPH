module globvar
	
	use derived_types
	implicit none
	
	! particle array
	type(particles),allocatable,target:: parts(:)
	
	! interaction array
	type(interactions),allocatable:: pairs(:)
	
	!global 1D variables
	integer:: ntotal,nvirt
	integer:: maxn,maxinter
	integer:: niac
	integer:: itimestep,maxtimestep,save_step,print_step
	real(8):: time
	
	real(8):: scale_k
	
	!timing
	real(8):: cputime,output_time,test_time
	
	contains
	
		subroutine allocateGlobalArrays
		
			implicit none
			
			maxn = ntotal+nvirt
			maxinter = 12*maxn
			
			allocate( parts(maxn) )
			allocate( pairs(maxinter) )
			
		end subroutine allocateGlobalArrays
		
		subroutine deallocateGlobalArrays
		
			implicit none
			
			deallocate( parts )
			deallocate( pairs )
			
		end subroutine deallocateGlobalArrays
	
end module globvar