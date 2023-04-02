program SPH

   use datatypes, only: particles, interactions, time_tracking
   use input_m, only: read_input_and_allocate
   use iso_fortran_env, only: output_unit
   use kernel_m, only: kernel_k
   use ORB_m, only: partition_track, ORB
   use param, only: f, skf, dim
   use summary_m, only: print_summary, time_print
   use time_integration_m, only: time_integration
   use output_m, only: output
   use flap, only: command_line_interface

   implicit none
   integer:: thisImage, numImages, maxtimestep, print_step, save_step, ntotal, nvirt, maxnloc, maxinter, i
   integer, codimension[*]:: ntotal_loc, nhalo_loc, nvirt_loc
   real(f):: scale_k, bounds_loc(2*dim)
   integer, allocatable:: nexti(:)
   type(particles), allocatable, codimension[:]:: parts(:)
   type(interactions), allocatable:: pairs(:)
   type(time_tracking):: timings
   type(command_line_interface):: argparse
   integer:: argerror

   ! Saving image id and total number of images
   thisImage = this_image()
   numImages = num_images()

   ! Parsing arguments
   call argparse%init(description='GraSPH is a program to simulate fluids and granular flows using the SPH numerical ' &
                                //'method.', &
                      version = 'v0.1.0', &
                      authors = 'Edward Yang', &
                      license = 'MIT', &
                      help = 'Usage: ', &
                      examples = ['./sph -m 10000 -p 100 -s 100 -o outputdata -i example/dambreak.h5'])

   call argparse%add(switch='--maxtimestep', &
                     switch_ab='-m', &
                     help='Max number of time-steps to run the simulation for.', &
                     required=.true., &
                     act='store', &
                     error=argerror)
   if (argerror /= 0) stop
           
   call argparse%add(switch='--printtimestep', &
                     switch_ab='-p', &
                     help='Interval of time-steps to write to terminal e.g., -p 1000 will print information every ' &
                        //'1000 time-steps.', &
                     required=.true., &
                     act='store', &
                     error=argerror)
   if (argerror /= 0) stop

   call argparse%add(switch='--savetimestep', &
                     switch_ab='-s', &
                     help='Interval of time-steps to save data e.g., -s 1000 will save output data every 1000 ' &
                        //'time-steps. Data is saved to <output-directory>/sph_out*.h5.', &
                     required=.true., &
                     act='store', &
                     error=argerror)
   if (argerror /= 0) stop

   call argparse%get(switch='-m', val=maxtimestep, error=argerror)
   if (argerror /= 0) stop
   call argparse%get(switch='-p', val=print_step, error=argerror)
   if (argerror /= 0) stop
   call argparse%get(switch='-s', val=save_step, error=argerror)
   if (argerror /= 0) stop

   ! printing preamble to screen
   if (thisImage == 1) call time_print

   if (numImages .eq. 1) then
      write (*, '(A)') 'Executing code in serial!'
   else
      if (thisImage .eq. 1) write (*, '(A,I4,A)') 'Executing code in parallel with ', numImages, ' images!'
   end if

   if (thisImage == 1) then
      write (output_unit, '(A,I7,A)') 'Running ', maxtimestep, ' step(s).'
      write (output_unit, '(A,I7,A)') 'Printing summary to screen every ', print_step, ' step(s).'
      write (output_unit, '(A,I7,A)') 'Writing output to disc every ', save_step, ' step(s).'
   end if

   ! Retrieving kernel k parameter for use in the program
   scale_k = kernel_k(skf)

   ! Reading input hdf5 file and allocating persistent arrays
   call read_input_and_allocate(thisImage, numImages, ntotal, ntotal_loc, nvirt, nvirt_loc, parts, pairs, nexti, &
      maxnloc, maxinter)

   if (thisImage .eq. 1) then
      write (*, '(A,I0,A)') 'Total simulation size of ', ntotal, ' physical particles, and'
      write (*, '(A,I0,A)') '                         ', nvirt, ' virtual particles.'
   end if

   nhalo_loc = 0
   call ORB(0, thisImage, numImages, scale_k, ntotal, ntotal_loc, nvirt, nvirt_loc, nhalo_loc, parts, timings)
   call output(0, save_step, thisImage, numImages, ntotal_loc, nhalo_loc, nvirt_loc, parts, ntotal)

   ! Entering discretized time-integration loop
   call time_integration(maxtimestep, print_step, save_step, thisImage, numImages, maxnloc, maxinter, timings, scale_k, &
                         ntotal_loc, nvirt_loc, nhalo_loc, ntotal, nvirt, parts, pairs, nexti)

   !Printing post-amble to terminal
   call print_summary(thisImage, numImages, timings, partition_track)

end program SPH
