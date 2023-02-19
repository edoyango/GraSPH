module summary_m

   use iso_fortran_env, only: output_unit, error_unit, event_type
   use datatypes, only: time_tracking
   use param, only: f, tf

   private
   public:: preamble, print_summary, print_loadbalance

contains

   !====================================================================================================================
   subroutine time_print
      !***************************************************************************************************
      !   TIME_PRINT:      Print out the current date and time.
      !   The standard Fortran 90 routine DATE_AND_TIME is used to get the current
      !   date and time strings.
      !***************************************************************************************************
      implicit none

      !Local scalars.
      character(len=8) :: datstr
      character(len=10) :: timstr

      !Get the current date and time.
      call date_and_time(datstr, timstr)

      !Write out the date and time.
      write (output_unit, "(/31x,A)") "Date = "//datstr(7:8)//"/"//datstr(5:6)//"/"//datstr(1:4)
      write (output_unit, "(30x,A)") "Time = "//timstr(1:2)//":"//timstr(3:4)//":"//timstr(5:10)

   end subroutine time_print

   !====================================================================================================================
   subroutine preamble(thisImage, numImages, maxtimestep, print_step, save_step)
      ! Prints preamble information to terminal and checks that maxtimestep, print_step, save_step have been supplied

      implicit none
      integer, intent(in):: thisImage, numImages
      integer, intent(out):: maxtimestep, print_step, save_step
      character(len=100):: args(3)

      if (command_argument_count() .lt. 3) then
         if (thisImage .eq. 1) then
            write (error_unit, '(A79)') '!!!ERROR!!! --- 3 commandline arguments must be provided: maxtimestep,'
            write (error_unit, '(A56)') '                print_step, save_step. Program ending...'
            stop
         end if
         error stop
      else
         call get_command_argument(1, args(1))
         call get_command_argument(2, args(2))
         call get_command_argument(3, args(3))
         read (args, *) maxtimestep, print_step, save_step
      end if

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

   end subroutine preamble

   !====================================================================================================================
   subroutine print_summary(thisImage, numImages, timings, partition_track)

      ! Obtains and prints final MPI summary data e.g. number of partitions, cut axes reorientations, and average 
      ! wall-times (broken down)

      use param_para, only: partition_tracking

      implicit none
      integer, intent(in):: thisImage, numImages
      type(time_tracking), intent(in):: timings
      type(partition_tracking), intent(in):: partition_track
      real(tf):: t_wall_avg, t_ORB_avg, t_dist_avg, t_output_avg

      t_wall_avg = timings%t_wall
      t_ORB_avg = timings%t_ORB
      t_dist_avg = timings%t_dist
      t_output_avg = timings%t_output

      call co_sum(t_wall_avg, result_image=1)
      call co_sum(t_ORB_avg, result_image=1)
      call co_sum(t_dist_avg, result_image=1)
      call co_sum(t_output_avg, result_image=1)

      if (thisImage==1) then

         call time_print
         t_wall_avg = t_wall_avg/numImages
         t_ORB_avg = t_ORB_avg/numImages
         t_dist_avg = t_dist_avg/numImages
         t_output_avg = t_output_avg/numImages

         write(*, '(A)') '================================= TIME SUMMARY ================================'
         write (*, '(A,F14.7)') 'Average Wall time (s)      = ', t_wall_avg
         write (*, '(A,F14.7)') 'Average Partition time (s) = ', t_ORB_avg
         write (*, '(A,F14.7)') 'Average Send/recv time (s) = ', t_dist_avg
         write (*, '(A,F14.7)') 'Average Output time (s)    = ', t_output_avg
         write (*, '(A)') '============================== PARTITION SUMMARY =============================='
         write (*, '(A,I0)') '            Number of Partitions = ', partition_track%n_parts
         if (partition_track%n_parts .eq. 1) then
            write (*, '(A)') '    Avg Timesteps B/N Partitions = N/A'
            write (*, '(A)') '    Max Timesteps B/N Partitions = N/A'
            write (*, '(A)') '    Min Timesteps B/N Partitions = N/A'
         else
            write (*, '(A,I0)') '    Avg Timesteps B/N Partitions = ', (partition_track%prev_part_tstep - 1)/ &
               (partition_track%n_parts - 1)
            write (*, '(A,I0)') '    Max Timesteps B/N Partitions = ', partition_track%maxtstep_bn_part
            write (*, '(A,I0)') '    Min Timesteps B/N Partitions = ', partition_track%mintstep_bn_part
         end if

         write (*, *)
         write (*, '(A,I0)') 'Number of Cut Axis Reorientation = ', partition_track%n_reorients

         if (partition_track%n_reorients .eq. 1) then
            write (*, '(A)') 'Avg Timesteps B/N Reorientations = N/A'
            write (*, '(A)') 'Max Timesteps B/N Reorientations = N/A'
            write (*, '(A)') 'Min Timesteps B/N Reorientations = N/A'
         else
            write (*, '(A,I0)') 'Avg Timesteps B/N Reorientations = ', (partition_track%prev_reorient_tstep - 1)/ &
               (partition_track%n_reorients - 1)
            write (*, '(A,I0)') 'Max Timesteps B/N Reorientations = ', partition_track%maxtstep_bn_reorient
            write (*, '(A,I0)') 'Min Timesteps B/N Reorientations = ', partition_track%mintstep_bn_reorient
         end if

      end if

   end subroutine print_summary

   !==================================================================================================================================
   subroutine print_loadbalance(thisImage, numImages, wtime, ntotal_loc, nhalo_loc, nvirt_loc, niac, itimestep, &
         time, maxtimestep)
      ! Prints load balance statistics. Called occasionally as determined by print_step

      implicit none
      integer, intent(in):: thisImage, numImages, ntotal_loc, nhalo_loc, nvirt_loc, niac, itimestep, maxtimestep
      real(tf), intent(in):: wtime
      real(f), intent(in):: time
      integer:: i, min_n(4), max_n(4), mean_n(4), stdev_n(4)
      integer, codimension[:], allocatable:: n_glob(:, :)
      type(event_type), codimension[*], save:: putevents

      allocate(n_glob(numImages, 4)[*])

      n_glob(thisImage, 1)[1] = ntotal_loc
      n_glob(thisImage, 2)[1] = nhalo_loc
      n_glob(thisImage, 3)[1] = nvirt_loc
      n_glob(thisImage, 4)[1] = niac

      event post (putevents[1])

      if (thisImage==1) then
         event wait (putevents, until_count=numImages)

         !calculating summary statistics for load balance
         do i = 1, 4
            min_n(i) = MINVAL(n_glob(:, i))
            max_n(i) = MAXVAL(n_glob(:, i))
            mean_n(i) = SUM(n_glob(:, i))/numImages
         end do

         stdev_n(:) = 0._f
         do i = 1, numImages
            stdev_n(:) = stdev_n(:) + (n_glob(i, :) - mean_n(:))**2
         end do
         stdev_n(:) = int(SQRT(ABS(real(stdev_n(:), kind=f))/numImages))

         write (*, '(A)') '_______________________________________________________________________________'
         write (*, '(A,I7,A,I7,9x,A,F14.7)') '  current time step = ', itimestep, '/', maxtimestep, &
            'current time (s) = ', time
         write (*, '(A,F14.7)') '                                                  Walltime (s) = ', wtime
         write (*, '(A)') '_______________________________________________________________________________'
         write (*, 9999) 'ntotal_loc statistics: mean = ', mean_n(1), ' stdev = ', stdev_n(1)
         write (*, 9999) '                       min  = ', min_n(1), ' max   = ', max_n(1)
         write (*, 9999) ' nhalo_loc statistics: mean = ', mean_n(2), ' stdev = ', stdev_n(2)
         write (*, 9999) '                       min  = ', min_n(2), ' max   = ', max_n(2)
         write (*, 9999) ' nvirt_loc statistics: mean = ', mean_n(3), ' stdev = ', stdev_n(3)
         write (*, 9999) '                       min  = ', min_n(3), ' max   = ', max_n(3)
         write (*, 9999) '      niac statistics: mean = ', mean_n(4), ' stdev = ', stdev_n(4)
         write (*, 9999) '                       min  = ', min_n(4), ' max   = ', max_n(4)

9999     format(A, I10, A, I10)

      end if

   end subroutine print_loadbalance

end module summary_m
