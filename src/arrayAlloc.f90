subroutine allocateGlobalArrays( )

use param
implicit none

allocate( parts(ntotal+nvirt) )
			
maxinter = 11*(ntotal+nvirt)
allocate( pairs(maxinter) )
			
end

subroutine deallocateGlobalArrays( )

use param
implicit none

deallocate( parts )
deallocate( pairs )

end