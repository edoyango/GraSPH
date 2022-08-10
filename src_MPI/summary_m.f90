module summary_m

   use datatypes, only: time_tracking
   use param, only: f, tf
   use mpi_f08

contains

   !==============================================================================================================================
   subroutine time_print
      !***************************************************************************************************
      !   TIME_PRINT:      Print out the current date and time.
      !   The standard Fortran 90 routine DATE_AND_TIME is used to get the current
      !   date and time strings.
      !***************************************************************************************************
      implicit none
      integer, parameter :: output = 6

      !Local scalars.
      character(len=8) :: datstr
      character(len=10) :: timstr

      !Get the current date and time.
      call date_and_time(datstr, timstr)

      !Write out the date and time.
      write (output, "(/31x,A)") "Date = "//datstr(7:8)//"/"//datstr(5:6)//"/"//datstr(1:4)
      write (output, "(30x,A)") "Time = "//timstr(1:2)//":"//timstr(3:4)//":"//timstr(5:10)

   end subroutine time_print

   !==============================================================================================================================
   subroutine preamble(procid, numprocs, maxtimestep, print_step, save_step)
      ! Prints preamble information to terminal and checks that maxtimestep, print_step, save_step have been supplied

      implicit none
      integer, intent(in):: procid, numprocs
      integer, intent(out):: maxtimestep, print_step, save_step
      integer:: ierr
      character(len=100):: args(3)

      if (command_argument_count() .lt. 3) then
         if (procid .eq. 0) then
            write (*, '(A79)') '!!!ERROR!!! --- 3 commandline arguments must be provided: maxtimestep,'
            write (*, '(A56)') '                print_step, save_step. Program ending...'
            stop
         end if
         call MPI_ABORT(MPI_COMM_WORLD, 100, ierr)
      else
         call get_command_argument(1, args(1))
         call get_command_argument(2, args(2))
         call get_command_argument(3, args(3))
         read (args, *) maxtimestep, print_step, save_step
      end if

      if (procid .eq. 0) call time_print

      if (numprocs .eq. 1) then
         write (*, '(A25)') 'Executing code in serial!'
      else
         if (procid .eq. 0) write (*, '(A31,1x,I4,1x,A10)') 'Executing code in parallel with', numprocs, 'processes!'
      end if

      if (procid .eq. 0) then
         write (*, '(A8,I7,A9)') 'Running ', maxtimestep, ' step(s).'
         write (*, '(A33,I7,A9)') 'Printing summary to screen every ', print_step, ' step(s).'
         write (*, '(A29,I7,A9)') 'Writing output to disc every ', save_step, ' step(s).'
      end if

   end subroutine preamble

   !==============================================================================================================================
   subroutine print_summary(procid, numprocs, timings, partition_track)
      ! Obtains and prints final MPI summary data e.g. number of partitions, cut axes reorientations, and average wall-times (broken
      ! down)

      use param_para, only: partition_tracking

      implicit none
      integer, intent(in):: procid, numprocs
      type(time_tracking), intent(in):: timings
      type(partition_tracking), intent(in):: partition_track
      integer:: ierr
      double precision:: t_wall_avg, t_ORB_avg, t_dist_avg, t_output_avg

      if (procid .eq. 0) then
         call MPI_REDUCE(timings%t_wall, t_wall_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(timings%t_ORB, t_ORB_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(timings%t_dist, t_dist_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(timings%t_output, t_output_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         t_wall_avg = t_wall_avg/numprocs
         t_ORB_avg = t_ORB_avg/numprocs
         t_dist_avg = t_dist_avg/numprocs
         t_output_avg = t_output_avg/numprocs

         write (*, '(A79)') '================================= TIME SUMMARY ================================'
         write (*, '(A29,F14.7)') 'Average Wall time (s)      = ', t_wall_avg
         write (*, '(A29,F14.7)') 'Average Partition time (s) = ', t_ORB_avg
         write (*, '(A29,F14.7)') 'Average Send/recv time (s) = ', t_dist_avg
         write (*, '(A29,F14.7)') 'Average Output time (s)    = ', t_output_avg
         write (*, '(A79)') '============================== PARTITION SUMMARY =============================='
         write (*, '(A35,I7)') '            Number of Partitions = ', partition_track%n_parts
         if (partition_track%n_parts .eq. 1) then
            write (*, '(A42)') '    Avg Timesteps B/N Partitions =     N/A'
            write (*, '(A42)') '    Max Timesteps B/N Partitions =     N/A'
            write (*, '(A42)') '    Min Timesteps B/N Partitions =     N/A'
         else
            write (*, '(A35,I7)') '    Avg Timesteps B/N Partitions = ', (partition_track%prev_part_tstep - 1)/ &
               (partition_track%n_parts - 1)
            write (*, '(A35,I7)') '    Max Timesteps B/N Partitions = ', partition_track%maxtstep_bn_part
            write (*, '(A35,I7)') '    Min Timesteps B/N Partitions = ', partition_track%mintstep_bn_part
         end if

         write (*, *)
         write (*, '(A35,I7)') 'Number of Cut Axis Reorientation = ', partition_track%n_reorients

         if (partition_track%n_reorients .eq. 1) then
            write (*, '(A42)') 'Avg Timesteps B/N Reorientations =     N/A'
            write (*, '(A42)') 'Max Timesteps B/N Reorientations =     N/A'
            write (*, '(A42)') 'Min Timesteps B/N Reorientations =     N/A'
         else
            write (*, '(A35,I7)') 'Avg Timesteps B/N Reorientations = ', (partition_track%prev_reorient_tstep - 1)/ &
               (partition_track%n_reorients - 1)
            write (*, '(A35,I7)') 'Max Timesteps B/N Reorientations = ', partition_track%maxtstep_bn_reorient
            write (*, '(A35,I7)') 'Min Timesteps B/N Reorientations = ', partition_track%mintstep_bn_reorient
         end if
      else
         call MPI_REDUCE(timings%t_wall, t_wall_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(timings%t_ORB, t_ORB_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(timings%t_dist, t_dist_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(timings%t_output, t_output_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      end if

   end subroutine print_summary

   !==================================================================================================================================
   subroutine print_loadbalance(procid, numprocs, timings, ntotal_loc, nhalo_loc, nvirt_loc, niac, itimestep, time, maxtimestep)
      ! Prints load balance statistics. Called occasionally as determined by print_step

      implicit none
      integer, intent(in):: procid, numprocs, ntotal_loc, nhalo_loc, nvirt_loc, niac, itimestep, maxtimestep
      type(time_tracking), intent(in):: timings
      real(f), intent(in):: time
      integer:: i, n_glob(numprocs, 4), min_n(4), max_n(4), mean_n(4), stdev_n(4), ierr
      type(MPI_Request):: request(4)
      type(MPI_Status):: status(MPI_STATUS_SIZE*4)

      call MPI_IGATHER(ntotal_loc, 1, MPI_INTEGER, n_glob(:, 1), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, request(1), ierr)
      call MPI_IGATHER(nhalo_loc, 1, MPI_INTEGER, n_glob(:, 2), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, request(2), ierr)
      call MPI_IGATHER(nvirt_loc, 1, MPI_INTEGER, n_glob(:, 3), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, request(3), ierr)
      call MPI_IGATHER(niac, 1, MPI_INTEGER, n_glob(:, 4), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, request(4), ierr)

      call MPI_WAITALL(4, request, status, ierr)

      if (procid .eq. 0) then

         !calculating summary statistics for load balance
         do i = 1, 4
            min_n(i) = MINVAL(n_glob(:, i))
            max_n(i) = MAXVAL(n_glob(:, i))
            mean_n(i) = SUM(n_glob(:, i))/numprocs
         end do

         stdev_n(:) = 0_f
         do i = 1, numprocs
            stdev_n(:) = stdev_n(:) + (n_glob(i, :) - mean_n(:))**2
         end do
         stdev_n(:) = int(SQRT(ABS(DBLE(stdev_n(:)))/numprocs))

         write (*, '(A79)') '_______________________________________________________________________________'
         write (*, '(A22,I7,A1,I7,9x,A19,F14.7)') '  current time step = ', itimestep, '/', maxtimestep, 'current time (s) = ', time
         write (*, '(A65,F14.7)') '                                                  Walltime (s) = ', timings%t_wall
         write (*, '(A79)') '_______________________________________________________________________________'
         write (*, 9999) 'ntotal_loc statistics: mean = ', mean_n(1), ' stdev = ', stdev_n(1)
         write (*, 9999) '                       min  = ', min_n(1), ' max   = ', max_n(1)
         write (*, 9999) ' nhalo_loc statistics: mean = ', mean_n(2), ' stdev = ', stdev_n(2)
         write (*, 9999) '                       min  = ', min_n(2), ' max   = ', max_n(2)
         write (*, 9999) ' nvirt_loc statistics: mean = ', mean_n(3), ' stdev = ', stdev_n(3)
         write (*, 9999) '                       min  = ', min_n(3), ' max   = ', max_n(3)
         write (*, 9999) '      niac statistics: mean = ', mean_n(4), ' stdev = ', stdev_n(4)
         write (*, 9999) '                       min  = ', min_n(4), ' max   = ', max_n(4)

9999     format(A30, I10, A9, I10)

      end if

   end subroutine print_loadbalance

end module summary_m
