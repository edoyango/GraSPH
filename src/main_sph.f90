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

use param
implicit none
real(8):: Hmax

call time_print

! defining kernel parameters
if (skf.eq.1) scale_k = 2d0
if (skf.eq.2) scale_k = 3d0
current_ts = 0
nstart = 0

mp = 50    ! number of real particles in x
np = 50    ! number of real particles in y
op = 150   ! number of virtual particles in x
pp = 80	   ! number of virtual particles in y
ntotal = mp*np 		! total number of real particles
nvirt = 2*op + 2*pp	! total number of virtual particles

call allocateGlobalArrays( )

!Created fluid domain using particles
call input( )
	
!Created virtual boundary particles
call virt_part( )

1 write(*,*)'  ***************************************************' 
  write(*,*)'          Please input the maximal time steps '
  write(*,*)'  ***************************************************'
  read(*,*) maxtimestep
  
call time_integration( )

write(*,*)'  ***************************************************'
write(*,*) 'Are you going to run more time steps ? (0=No, 1=yes)'
write(*,*)'  ***************************************************'
read (*,*) yesorno
if (yesorno.ne.0) go to 1

call time_print
write (*,*)'   Elapsed CPU time = ', cputime

call deallocateGlobalArrays( )

end