module summary_m

   use datatypes, only: particles, interactions, time_tracking
   use param, only: f

   public:: time_print, preamble, print_summary, print_update

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

   !==================================================================================================================================
   subroutine preamble(maxtimestep, print_step, save_step)

      implicit none
      integer, intent(out):: maxtimestep, print_step, save_step
      character(len=100):: args(3)

      if (command_argument_count() .lt. 3) then
         write (*, '(A79)') '!!!ERROR!!! --- 3 commandline arguments must be provided: maxtimestep,'
         write (*, '(A56)') '                print_step, save_step. Program ending...'
         stop
      else
         call get_command_argument(1, args(1))
         call get_command_argument(2, args(2))
         call get_command_argument(3, args(3))
         read (args, *) maxtimestep, print_step, save_step
      end if

      call time_print

      write (*, '(A36)') 'Serial code - running on one core...'

      write (*, '(A8,I7,A9)') 'Running ', maxtimestep, ' step(s).'
      write (*, '(A33,I7,A9)') 'Printing summary to screen every ', print_step, ' step(s).'
      write (*, '(A29,I7,A9)') 'Writing output to disc every ', save_step, ' step(s).'

   end subroutine preamble

   !==================================================================================================================================
   subroutine print_summary(timings)

      implicit none
      type(time_tracking), intent(in):: timings

      write (*, '(A79)') '================================= TIME SUMMARY ================================'
      write (*, '(A29,F14.7)') 'Average Wall time (s)   = ', timings%walltime()
      write (*, '(A29,F14.7)') 'Average Output time (s) = ', timings%t_output

   end subroutine print_summary

   !==================================================================================================================================
   subroutine print_update(itimestep, maxtimestep, ntotal, nvirt, nghos, niac, parts, pairs, nexti, time, timings)

      implicit none
      integer, intent(in):: ntotal, nvirt, nghos, niac, itimestep, maxtimestep, nexti(:)
      type(particles), intent(in):: parts(:)
      type(interactions), intent(in):: pairs(:)
      type(time_tracking), intent(in):: timings
      real(f), intent(in):: time
      integer, allocatable:: ns(:)
      integer:: maxp, minp, sumiac, maxiac, miniac, noiac, i, j, k

      write (*, '(A79)') '_______________________________________________________________________________'
      write (*, '(A22,I7,A1,I7,9x,A19,F14.7)') '  current time step = ', itimestep, '/', maxtimestep, &
         'current time (s) = ', time
      write (*, '(A65,F14.7)') '                                                  Walltime (s) = ', timings%walltime()
      write (*, '(A79)') '_______________________________________________________________________________'

      !Statistics for the interaction, Print information to screen
      allocate (ns(ntotal + nvirt + nghos))
      ns(:) = 0
      do i = 1, ntotal + nvirt + nghos
         do k = nexti(i), nexti(i + 1) - 1
            j = pairs(k)%j
            ns(parts(i)%ind) = ns(parts(i)%ind) + 1
            ns(parts(j)%ind) = ns(parts(j)%ind) + 1
         end do
      end do

      sumiac = 0
      maxiac = -huge(1)
      noiac = 0
      miniac = huge(1)
      maxp = 1
      do i = 1, ntotal + nvirt + nghos
         sumiac = sumiac + ns(parts(i)%ind)
         if (ns(parts(i)%ind) > maxiac) then
            maxiac = ns(parts(i)%ind)
            maxp = i
         end if
         if (ns(parts(i)%ind) < miniac) then
            miniac = ns(parts(i)%ind)
            minp = i
         end if
         if (ns(parts(i)%ind) .eq. 0) noiac = noiac + 1
      end do

      write (*, '(A40,I7)') ' >> Statistics: no. of ghost particles: ', nghos
      write (*, '(A42)') ' >> Statistics: interactions per particle:'
      write (*, '(A14,I7,A23,I4)') '    Particle: ', maxp, ' maximal interactions: ', maxiac
      write (*, '(A14,I7,A23,I4)') '    Particle: ', minp, ' minimal interactions: ', miniac
      write (*, '(A14,F14.7)') '    Average : ', real(sumiac)/real(ntotal + nvirt + nghos)
      write (*, '(A18,I9)') '    Total pairs : ', sumiac
      write (*, '(A36,I7)') '    Particles with no interactions: ', noiac

   end subroutine print_update

end module summary_m
