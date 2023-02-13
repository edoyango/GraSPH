module time_integration_m

   use datatypes, only: particles, interactions, time_tracking
   use mpi_f08
   use param, only: f, tf, dim, dt, rh0, c, gamma, g
   use param_para, only: MPI_derived_types

   use input_m, only: update_virt_part
   use flink_list_m, only: flink_list
   use ORB_m, only: ORB, neighbours, n_process_neighbour
   use ORB_sr_m, only: ORB_sendrecv_haloupdate
   use output_m, only: output
   use single_step_m, only: single_step
   use summary_m, only: print_loadbalance

   private
   public:: time_integration

contains

   !==============================================================================================================================
   subroutine time_integration(maxtimestep, print_step, save_step, procid, numprocs, maxnloc, maxinter, MPI_types, timings, &
                               scale_k, ntotal_loc, nvirt_loc, nhalo_loc, nghos_loc, ntotal, parts, pairs, nexti, gind)
      ! Subroutine responsible for the main time-integration loop

      implicit none
      integer, intent(in):: procid, numprocs, maxtimestep, print_step, save_step, maxnloc, maxinter, ntotal
      type(MPI_derived_types), intent(in):: MPI_types
      real(f), intent(in):: scale_k
      type(time_tracking), intent(inout):: timings
      integer, intent(inout):: ntotal_loc, nvirt_loc, nghos_loc, nhalo_loc, nexti(:), gind(:)
      type(particles), intent(inout):: parts(maxnloc)
      type(interactions), intent(out):: pairs(maxinter)
      integer:: i, ki, itimestep, niac
      real(f):: time = 0._f
      double precision:: tmptime ! used to record wall time measurements
      real(f), allocatable:: v_min(:, :), rho_min(:), dvxdt(:, :), drho(:), vw(:)

      allocate (v_min(dim, maxnloc), rho_min(maxnloc), dvxdt(dim, maxnloc), drho(maxnloc), vw(maxnloc))

      drho(:) = 0._f
      dvxdt(1:dim-1,:) = 0._f
      dvxdt(dim, :) = -g
      
      timings%t_wall = timings%t_wall - MPI_WTIME()

      ! Time-integration (Leap-Frog)
      do itimestep = 1, maxtimestep

         ! save properties at start of step, update properties to mid-step.
         do i = 1,ntotal_loc
            parts(i)%v_min(:) = parts(i)%vx(:)
            parts(i)%vx(:) = parts(i)%vx(:) + 0.5_f*dt*dvxdt(:,i)
            parts(i)%rho_min = parts(i)%rho
            parts(i)%rho = parts(i)%rho + 0.5_f*dt*drho(i)
         end do

         ! distributing particles
         call ORB(itimestep, procid, numprocs, MPI_types, scale_k, ntotal, ntotal_loc, nhalo_loc, nvirt_loc, parts, timings)

         !Interaction parameters, calculating neighboring particles
         call flink_list(maxinter, scale_k, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, niac, parts, pairs, nexti)

            ! update halo particles after first increment
            ! if (ki > 1) then
            !    timings%t_dist = timings%t_dist - MPI_WTIME()
            !    call ORB_sendrecv_haloupdate(ki, MPI_types%haloupdatetype, n_process_neighbour, neighbours, ntotal_loc, parts)
            !    timings%t_dist = timings%t_dist + MPI_WTIME()
            ! end if

            call update_virt_part(ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, parts, niac, pairs, nexti, vw)

            ! update pressure of newly updated real and halo particles
            do i = 1, ntotal_loc + nhalo_loc + nvirt_loc
               parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1_f)/gamma
            end do

            ! calculating forces
            call single_step(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, dvxdt, drho, nexti)

            ! call RK4_update(ki, ntotal_loc, v_min, rho_min, dvxdt, drho, parts)

            do i = 1,ntotal_loc
               parts(i)%vx(:) = parts(i)%v_min(:) + dt*dvxdt(:,i)
               parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)
               parts(i)%rho = parts(i)%rho_min + dt*drho(i)
            end do

         time = time + dt

         timings%t_output = timings%t_output - MPI_WTIME()

         ! write output data
         if (mod(itimestep, save_step) .eq. 0) then
            call output(itimestep, save_step, procid, numprocs, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, parts, ntotal)
         end if

         if (mod(itimestep, print_step) .eq. 0) then
            tmptime = MPI_WTIME()
            timings%t_wall = timings%t_wall + tmptime
            call print_loadbalance(procid, numprocs, timings, ntotal_loc, nhalo_loc, nvirt_loc, niac, itimestep, time, maxtimestep)
            timings%t_wall = timings%t_wall - tmptime
         end if
         
         timings%t_output = timings%t_output + MPI_WTIME()

      end do
      
      timings%t_wall = timings%t_wall + MPI_WTIME()

   end subroutine time_integration

   !================================================================================================================================
   pure subroutine RK4_update(ki, ntotal_loc, v_min, rho_min, dvxdt, drhodt, parts)

      implicit none
      integer, intent(in):: ki, ntotal_loc
      real(f), intent(in):: v_min(:, :), rho_min(:), dvxdt(:, :, :), drhodt(:, :)
      type(particles), intent(inout):: parts(:)
      integer:: i

      select case (ki)
      case default
         do i = 1, ntotal_loc
            parts(i)%vx(:) = v_min(:, i) + 0.5_f*dt*dvxdt(:, i, ki)
            parts(i)%rho = rho_min(i) + 0.5_f*dt*drhodt(i, ki)
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1._f)/gamma
         end do
      case (3)
         do i = 1, ntotal_loc
            parts(i)%vx(:) = v_min(:, i) + dt*dvxdt(:, i, 3)
            parts(i)%rho = rho_min(i) + dt*drhodt(i, 3)
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1._f)/gamma
         end do
      case (4)
         do i = 1, ntotal_loc
            parts(i)%vx(:) = v_min(:, i) + &
                             dt/6_f*(dvxdt(:, i, 1) + 2._f*dvxdt(:, i, 2) + 2._f*dvxdt(:, i, 3) + dvxdt(:, i, 4))
            parts(i)%rho = rho_min(i) + dt/6_f*(drhodt(i, 1) + 2._f*drhodt(i, 2) + 2._f*drhodt(i, 3) + drhodt(i, 4))
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1._f)/gamma
            parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)
         end do
      end select

   end subroutine RK4_update

end module time_integration_m
