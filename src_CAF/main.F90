program SPH

   use datatypes, only: particles, interactions, time_tracking
   use input_m, only: read_input_and_allocate
   use ORB_m, only: partition_track, ORB
#ifdef PARALLEL
   use mpi
#endif
   use param, only: f, skf, dim
   use summary_m, only: preamble, print_summary
   use time_integration_m, only: time_integration
   use output_m, only: output

   implicit none
   integer:: my_rank, num_ranks, maxtimestep, print_step, save_step, maxinter, i, ierr
   real(f):: bounds_loc(2*dim)
   integer, allocatable:: nexti(:)
   type(particles):: parts
   type(interactions), allocatable:: pairs(:)
   type(time_tracking):: timings

   ! Initializing MPI and obtaining this process' rank and number of MPI ranks participating
#ifdef PARALLEL
   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, num_ranks, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
#else
   num_ranks = 1
   my_rank = 0
#endif

   ! Printing preamble to screen
   call preamble(my_rank, num_ranks, maxtimestep, print_step, save_step)

   ! Reading input hdf5 file and allocating persistent arrays
   call read_input_and_allocate(my_rank, num_ranks, parts, pairs, nexti, maxinter)

   if (my_rank == 0) then
      write (*, '(A,I0,A)') 'Total simulation size of ', parts%ntotal, ' physical particles, and'
      write (*, '(A,I0,A)') '                         ', parts%nvirt, ' virtual particles.'
   end if

   parts%nhalo_loc = 0
   call ORB(0, my_rank, num_ranks, parts, timings)
   call output(0, save_step, my_rank, num_ranks, parts)

   ! Entering discretized time-integration loop
   call time_integration(maxtimestep, print_step, save_step, my_rank, num_ranks, maxinter, timings, &
                         parts, pairs, nexti)

   !Printing post-amble to terminal
   ! call print_summary(my_rank, num_ranks, timings, partition_track)

#ifdef PARALLEL
   call MPI_Finalize(ierr)
#endif

end program SPH
