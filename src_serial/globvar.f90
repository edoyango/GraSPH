module globvar
    
    use datatypes
    
    implicit none
    
    ! particle array
    type(particles),allocatable,target,public:: parts(:)
    
    ! interaction array
    type(interactions),allocatable:: pairs(:)
    
    !global 1D variables
    integer,public:: ntotal,nvirt,nghos
    integer,public:: maxn,maxinter
    integer,public:: niac
    integer,public:: itimestep,maxtimestep,save_step,print_step
    real(f),public:: time=0_f
    
    real(f),public:: scale_k
    
    !timing
    real(f),public:: cputime=0_f,output_time=0_f,test_time=0_f
    
    public:: allocateGlobalArrays,deallocateGlobalArrays

! subroutines to allocate and deallocate global arrays
contains
    
    !==============================================================================================================================
    subroutine allocateGlobalArrays
    
        implicit none
        
        maxn = 2*ntotal+nvirt
        maxinter = 262*maxn
        
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
