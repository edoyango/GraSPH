program SPH

   use datatypes, only: particles, interactions, time_tracking
   use input_m, only: read_input_and_allocate
   use kernel_m, only: kernel_k
   use ORB_m, only: partition_track, ORB
   use param, only: f, skf, dim
   use summary_m, only: preamble, print_summary
   use time_integration_m, only: time_integration
   use output_m, only: output

   implicit none
   integer:: thisImage, numImages, maxtimestep, print_step, save_step, ntotal, nvirt, maxnloc, maxinter, i
   integer, codimension[*]:: ntotal_loc, nhalo_loc, nvirt_loc
   real(f):: bounds_loc(2*dim)
   integer, allocatable:: nexti(:)
   type(particles), allocatable, codimension[:]:: parts(:)
   type(interactions), allocatable:: pairs(:)
   type(time_tracking):: timings

   ! Saving image id and total number of images
   thisImage = this_image()
   numImages = num_images()

   ! Printing preamble to screen
   call preamble(thisImage, numImages, maxtimestep, print_step, save_step)

   ! Reading input hdf5 file and allocating persistent arrays
   call read_input_and_allocate(thisImage, numImages, ntotal, ntotal_loc, nvirt, nvirt_loc, parts, pairs, nexti, &
      maxnloc, maxinter)

   if (thisImage .eq. 1) then
      write (*, '(A,I0,A)') 'Total simulation size of ', ntotal, ' physical particles, and'
      write (*, '(A,I0,A)') '                         ', nvirt, ' virtual particles.'
   end if

   nhalo_loc = 0
   call ORB(0, thisImage, numImages, ntotal, ntotal_loc, nvirt, nvirt_loc, nhalo_loc, parts, timings)
   call output(0, save_step, thisImage, numImages, ntotal_loc, nhalo_loc, nvirt_loc, parts, ntotal)

   ! Entering discretized time-integration loop
   call time_integration(maxtimestep, print_step, save_step, thisImage, numImages, maxnloc, maxinter, timings, &
                         ntotal_loc, nvirt_loc, nhalo_loc, ntotal, nvirt, parts, pairs, nexti)

   !Printing post-amble to terminal
   call print_summary(thisImage, numImages, timings, partition_track)

end program SPH
