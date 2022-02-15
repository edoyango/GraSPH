module globvar
	
	use datatypes
	
	implicit none
	
	! particle array
	type(particles),allocatable,target,public:: parts(:)
	
	! interaction array
	type(interactions),allocatable,public:: pairs(:)
	
	!global 1D variables
	integer,public:: ntotal,nvirt
	integer,public:: maxn,maxinter
	integer,public:: niac
	integer,public:: itimestep,maxtimestep,save_step,print_step
	real(f),public:: time
	
	real(f),public:: scale_k
	
	!timing
	real(f),public:: cputime,output_time,test_time
	
	public:: allocateGlobalArrays,deallocateGlobalArrays

! subroutines to allocate and deallocate global arrays
contains
	
	!==============================================================================================================================
	subroutine allocateGlobalArrays
	
		implicit none
		
		maxn = ntotal+nvirt
		maxinter = 50*maxn
		
		allocate( parts(maxn) )
		allocate( pairs(maxinter) )
		
	end subroutine allocateGlobalArrays
	
	!==============================================================================================================================
	subroutine deallocateGlobalArrays
	
		implicit none
		
		deallocate( parts )
		deallocate( pairs )
		
	end subroutine deallocateGlobalArrays
	
end module globvar