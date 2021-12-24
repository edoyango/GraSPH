!-----------------------------------------------------------------------------------------!
!    mmmmmmmmm 	      mmmmmmmmmmmm	  mmmm	   mmmm	   mmmmm          mmmmmmmmm           !
!   mmmmmmmmmm        mmmmmmmmmmmmm	  mmmm	   mmmm	   mmmmmmm        mmmmmmmmmm          !
!   mmmmmmmmm         mmmm	   mmmm	  mmmm	   mmmm	        mmmm      mmmm   mmmm         !
!     mmmmmmmm        mmmm	   mmmm	  mmmm	   mmmm	        mmmm      mmmm    mmmm        !
!       mmmmmmmm      mmmmmmmmmmmm	  mmmmmmmmmmmmm       mmmmm       mmmm    mmmm        !
!         mmmmmmmm    mmmmmmmmm		  mmmmmmmmmmmmm      mmmmm        mmmm    mmmm        !
!         mmmmmmmm    mmmmm			  mmmm	   mmmm	    mmmm          mmmm    mmmm        !
!        mmmmmmmmm    mmmmm			  mmmm	   mmmm	   mmmmm          mmmm   mmmm         !
!      mmmmmmmmm      mmmmm			  mmmm	   mmmm	   mmmmmmmmm      mmmmmmmmmm          !
!   mmmmmmmmm         mmmmm			  mmmm	   mmmm	   mmmmmmmmm      mmmmmmmm            !
!-----------------------------------------------------------------------------------------!
! ATTENTION:
! The copy right reserved. No part of this program may be reproduced, stored in a retrieval
! system, or transmitted in any form or by any means without the prior permission of:     !
!									BUI HONG HA                                           !
!					Monash Computational Geomechanics (MCG) Lab							  !
!				Department of Civil Engineering, Monash University						  !
!			    R132, 23 College Walk, Clayton, Vic3800, Australia						  !
!				  Tel: + 61 3 9905 2599  I  Fax: +61 3 9905 4944						  !
!				  Email: Ha.Bui@monash.edu I Web: www.monash.edu.au						  !
!---------------------------------- JANUARY 2012 -----------------------------------------!

program SPH

	use globvar
	use param

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