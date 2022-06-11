module summary_m

   use globvar, only: ntotal_loc, nhalo_loc, nvirt_loc, niac, maxtimestep, print_step, save_step, itimestep, time, cputime, &
                      t_graph, t_dist, output_time
   use globvar_para, only: procid, numprocs, mintstep_bn_part, mintstep_bn_reorient, maxtstep_bn_part, maxtstep_bn_reorient, &
                           prev_part_tstep, prev_reorient_tstep, n_parts, n_reorients, ierr
   use param, only: f
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
      write (output, "(/A)") "                  Date = "//datstr(7:8)//"/"// &
         datstr(5:6)//"/"// &
         datstr(1:4)
      write (output, "(A)") "                  Time = "//timstr(1:2)//":"// &
         timstr(3:4)//":"// &
         timstr(5:10)
      write (output, *)

   end subroutine time_print

   !==============================================================================================================================
   subroutine preamble
      ! Prints preamble information to terminal and checks that maxtimestep, print_step, save_step have been supplied

      implicit none
      character(len=100):: args(3)

      if (command_argument_count() .lt. 3) then
         if (procid .eq. 0) then
            write (*, '(A92)') '!!!ERROR!!! --- 3 commandline arguments must be provided: maxtimestep, print_step, save_step'
            write (*, '(A33)') '                Program ending...'
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
   subroutine print_summary
      ! Obtains and prints final MPI summary data e.g. number of partitions, cut axes reorientations, and average wall-times (broken
      ! down)

      implicit none
      double precision:: cputime_total, t_graph_total, t_dist_total, output_time_total

      if (procid .eq. 0) then
         call MPI_REDUCE(cputime, cputime_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(t_graph, t_graph_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(t_dist, t_dist_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(output_time, output_time_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         cputime = cputime/numprocs
         t_graph = t_graph/numprocs
         t_dist = t_dist/numprocs
         output_time = output_time/numprocs

         write (*, *) '================================= TIME SUMMARY ================================'
         write (*, *) 'Average Total CPU time = ', cputime
         write (*, *) 'Average Partition time = ', t_graph
         write (*, *) 'Average Send/recv time = ', t_dist
         write (*, *) 'Average Output time =    ', output_time
         write (*, *) '============================== PARTITION SUMMARY =============================='
         write (*, *) '            Number of Partitions = ', n_parts
         if (n_parts .eq. 1) then
            write (*, *) '    Avg Timesteps B/N Partitions =          N/A'
            write (*, *) '    Max Timesteps B/N Partitions =          N/A'
            write (*, *) '    Min Timesteps B/N Partitions =          N/A'
         else
            write (*, *) '    Avg Timesteps B/N Partitions = ', (prev_part_tstep - 1)/(n_parts - 1)
            write (*, *) '    Max Timesteps B/N Partitions = ', maxtstep_bn_part
            write (*, *) '    Min Timesteps B/N Partitions = ', mintstep_bn_part
         end if

         write (*, *) ' '
         write (*, *) 'Number of Cut Axis Reorientation = ', n_reorients

         if (n_reorients .eq. 1) then
            write (*, *) 'Avg Timesteps B/N Reorientations =          N/A'
            write (*, *) 'Max Timesteps B/N Reorientations =          N/A'
            write (*, *) 'Min Timesteps B/N Reorientations =          N/A'
         else
            write (*, *) 'Avg Timesteps B/N Reorientations = ', (prev_reorient_tstep - 1)/(n_reorients - 1)
            write (*, *) 'Max Timesteps B/N Reorientations = ', maxtstep_bn_reorient
            write (*, *) 'Min Timesteps B/N Reorientations = ', mintstep_bn_reorient
         end if
      else
         call MPI_REDUCE(cputime, cputime_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(t_graph, t_graph_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(t_dist, t_dist_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(output_time, output_time_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      end if

   end subroutine print_summary

   !==================================================================================================================================
   subroutine print_loadbalance
      ! Prints load balance statistics. Called occasionally as determined by print_step

      implicit none
      integer:: i, n_glob(numprocs, 4), min_n(4), max_n(4), mean_n(4), stdev_n(4)
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

         write (*, *) '_______________________________________________________________________________'
         write (*, *) '  current number of time step =', itimestep, '    current time=', real(time)
         write (*, *) '                                                 Walltime    =', &
            real(cputime + t_dist + t_graph + output_time)
         write (*, *) '_______________________________________________________________________________'
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
