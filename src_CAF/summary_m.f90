module summary_m

   use iso_fortran_env, only: output_unit, error_unit

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

end module summary_m