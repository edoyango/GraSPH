program SPH
   
   use datatypes, only: particles, interactions, time_tracking
   use param, only: f, skf
   use param_para, only: MPI_derived_types
   use ORB_m, only: partition_track, ORB
   
   use mpi_f08
   use output_m, only: output
   use input_m, only: return_ntotal, return_nvirt, allocatePersistentArrays, generate_real_part, generate_ghost_part
   use kernel_m, only: kernel_k
   use summary_m, only: preamble, time_print, print_summary
   use time_integration_m, only: time_integration

   implicit none
   integer:: maxtimestep, print_step, save_step, procid, numprocs, ntotal_loc = 0, nhalo_loc = 0, nvirt_loc = 0, nghos_loc = 0, &
   ntotal = 0, nvirt = 0, maxnloc, maxinter, ierr
   real(f):: scale_k
   type(MPI_derived_types):: MPI_types
   type(time_tracking):: timings
   integer,allocatable:: gind(:)
   type(particles),allocatable:: parts(:)
   type(interactions),allocatable:: pairs(:)

   !Initializing MPI
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, procid, ierr) !Retrieving process rank
   call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr) !Retrieving total number of processes

   ! Creating MPI datatypes for simulation
   call MPI_types%CreateMPITypes

   !Printing preamble to screen
   call preamble(procid, numprocs, maxtimestep, print_step, save_step)

   ! retrieving kernel k parameter for use in the program
   scale_k = kernel_k(skf)

   !Retrieving how many particles are to be generated. Hence the generate=.false.
   ntotal = return_ntotal()
   nvirt = return_nvirt()

   if (procid .eq. 0) write (*, '(A24,1x,I9,1x,A19)') 'Total simulation size of', ntotal, 'physical particles.'

   !Allocating particle and interaction arrays
   call allocatePersistentArrays(ntotal,nvirt,parts,pairs,gind,maxnloc,maxinter)

   ! Generating initial geometry, performing initial partition, and assigning virtual and ghost particles.
   call generate_real_part(procid, numprocs, ntotal, ntotal_loc, parts)
   call ORB(0, procid, numprocs, MPI_types, scale_k, ntotal, ntotal_loc, nhalo_loc, nvirt_loc, parts, timings)
   call generate_ghost_part(scale_k, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, parts, gind)
   call output(0, save_step, procid, numprocs, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, parts, ntotal)

   !Entering discretized time-integration loop
   call time_integration(maxtimestep, print_step, save_step, procid, numprocs, maxnloc, maxinter, MPI_types, timings, scale_k, &
   ntotal_loc, nvirt_loc, nhalo_loc, nghos_loc, ntotal, parts, pairs, gind)

   !Printing post-amble to terminal
   if (procid .eq. 0) call time_print
   call print_summary(procid, numprocs, timings, partition_track)

   !Cleaning up and ending MPI
   call MPI_FINALIZE(ierr)

end
