module error_msg_m

   public:: error_msg

contains

   !==============================================================================================================================
   subroutine error_msg(err_type, p_ind)
      ! subroutine to output pre-prepared messages.
      ! err_type = 1: When generating physical particles, array bounds would be exceeded
      !          = 2: When receiving physical particles from other processes, total number of particles would exceed array bounds
      !          = 3: When receiving halo particles from other processes, total number of particles would exceed array bounds
      !          = 4: When allocating virt particles, total number of particles would exceed array bounds
      !          = 5: a particle flagged for sending could not be sent to a process (adjacent comm)
      !          = 6: a particle flagged for sending could not be sent to a process (global comm)
      !          = 8: Too many particle interactions

      use globvar, only: itimestep
      use globvar_para, only: procid, ierr
      use mpi_f08

      use summary_m, only: print_summary

      implicit none
      integer, intent(in):: err_type, p_ind

      select case (err_type)
      case (1)
         write (*, '(A)') '>>>ERROR<<< Particle array is not large enough!'
         write (*, '(A)') '            Error invoked in inputMPI'
      case (2)
         write (*, '(A24,1x,I8,A9,1x,I4,A37)') '>>>ERROR<<< At time-step', itimestep, ', process', procid, &
            ': Particle array is not large enough!'
         write (*, '(A)') '            Error invoked in ORB_sendrecv_phys'
      case (3)
         write (*, '(A24,1x,I8,A9,1x,I4,A37)') '>>>ERROR<<< At time-step', itimestep, ', process', procid, &
            ': Particle array is not large enough!'
         write (*, '(A)') '            Error invoked in ORB_sendrecv_halo'
      case (4)
         write (*, '(A24,1x,I8,A9,1x,I4,A37)') '>>>ERROR<<< At time-step', itimestep, ', process', procid, &
            ': Particle array is not large enough!'
         write (*, '(A)') '            Error invoked in virt_part'
      case (5)
         write (*, '(A24,1x,I8,A9,1x,I4,1x,A32,1x,I8,A)') '>>>ERROR<<< At time-step', itimestep, ', process', procid, &
            ': Did not find home for particle', p_ind
         write (*, '(A)') '            Error invoked in ORB_sendrecv_phys'
      case (6)
         write (*, '(A24,1x,I8,A9,1x,I4,1x,A32,1x,I8,A)') '>>>ERROR<<< At time-step', itimestep, ', process', procid, &
            ': Did not find home for particle', p_ind
         write (*, '(A)') '            Error invoked in ORB_sendrecv_phys_total'
      case (7)
         write (*, '(A24,1x,I8,A9,1x,I4,1x,A32,1x,I8,A)') '>>>ERROR<<< At time-step', itimestep, ', process', procid, &
            ': Mismatched neighbours'
         write (*, '(A)') '            Error invoked in subdomain_actual_neighbour'
      case (8)
         write (*, '(A24,1x,I8,A9,1x,I4,1x,A32,1x,I8,A)') '>>>ERROR<<< At time-step', itimestep, ', process', procid, &
            ': Too many interactions'
      end select

      if (err_type .ne. 0) call print_summary

      call MPI_ABORT(MPI_COMM_WORLD, err_type, ierr)

   end subroutine error_msg

end module error_msg_m
