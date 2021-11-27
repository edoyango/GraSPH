subroutine allocateGlobalArrays( )

use param
use globvar
implicit none

allocate( parts(ntotal+nvirt) )
			
maxinter = 11*(ntotal+nvirt)
allocate( pairs(maxinter) )
			
end

subroutine deallocateGlobalArrays( )

use param
use globvar
implicit none

deallocate( parts )
deallocate( pairs )

end