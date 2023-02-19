program SPH

   use datatypes, only: particles, interactions, time_tracking
   use input_m, only: return_ntotal, return_nvirt, allocatePersistentArrays, generate_real_part
   use kernel_m, only: kernel_k
   use ORB_m, only: partition_track, ORB
   use output_m, only: output
   use param, only: f, skf
   use summary_m, only: preamble, print_summary
   use time_integration_m, only: time_integration

   implicit none
   integer:: thisImage, numImages, maxtimestep, print_step, save_step, ntotal, nvirt, maxnloc, maxinter, i
   integer, codimension[*]:: ntotal_loc, nhalo_loc, nvirt_loc
   real(f):: scale_k
   integer, allocatable:: nexti(:)
   type(particles), allocatable, codimension[:]:: parts(:)
   type(interactions), allocatable:: pairs(:)
   type(time_tracking):: timings

   ! Saving image id and total number of images
   thisImage = this_image()
   numImages = num_images()

   ! Printing preamble to screen
   call preamble(thisImage, numImages, maxtimestep, print_step, save_step)

   ! Retrieving kernel k parameter for use in the program
   scale_k = kernel_k(skf)

   ! Retrieving how many particles are to be generated.
   ntotal = return_ntotal()
   nvirt = return_nvirt()

   if (thisImage .eq. 1) write (*, '(A,1x,I9,1x,A)') 'Total simulation size of', ntotal, 'physical particles.'

   ! Allocating particle and interaction arrays
   call allocatePersistentArrays(ntotal, nvirt, parts, pairs, nexti, maxnloc, maxinter)

   ! Generating initial geometry, performing initial partition, and assigning virtual particles
   call generate_real_part(thisImage, numImages, ntotal, ntotal_loc, parts)
   call ORB(0, thisImage, numImages, scale_k, ntotal, ntotal_loc, nhalo_loc, nvirt_loc, parts, timings)

   call output(0, save_step, thisImage, numImages, ntotal_loc, nhalo_loc, nvirt_loc, parts, ntotal)

   ! Entering discretized time-integration loop
   call time_integration(maxtimestep, print_step, save_step, thisImage, numImages, maxnloc, maxinter, timings, scale_k, &
      ntotal_loc, nvirt_loc, nhalo_loc, ntotal, parts, pairs, nexti)

   !Printing post-amble to terminal
   call print_summary(thisImage, numImages, timings, partition_track)

end program SPH
