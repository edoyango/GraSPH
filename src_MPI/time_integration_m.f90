module time_integration_m

   use datatypes, only: particles, interactions
   use globvar, only: parts, ntotal_loc, time, cputime, output_time, t_graph, t_dist, test_time, itimestep, maxtimestep, &
                      print_step, save_step, maxnloc, nvirt_loc, nghos_loc, nhalo_loc, niac, pairs, maxinter, scale_k
   use mpi_f08
   use param, only: f, dim, dt, rh0, c, gamma
   use param_para, only: MPI_derived_types

   use input_m, only: gind, generate_ghost_part, update_ghost_part
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
   subroutine time_integration(procid, numprocs, MPI_types)
      ! Subroutine responsible for the main time-integration loop

      implicit none
      integer, intent(in):: procid, numprocs
      type(MPI_derived_types), intent(in):: MPI_types
      integer:: i, ki
      real(f), allocatable:: v_min(:, :), rho_min(:), dvxdt(:, :, :), drho(:, :)

      allocate (v_min(dim, maxnloc), rho_min(maxnloc), dvxdt(dim, maxnloc, 4), drho(maxnloc, 4))
      allocate (gind(maxnloc))

      ! Time-integration (Leap-Frog)
      do itimestep = 1, maxtimestep

         cputime = cputime - MPI_WTIME()

         ! distributing particles
         call ORB(procid, numprocs, MPI_types)

         do i = 1, ntotal_loc + nhalo_loc
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1_f)/gamma
         end do

         ! generating ghost particles
         call generate_ghost_part(ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, parts, gind)

         ! Storing velocity and density at initial time-step
         do i = 1, ntotal_loc
            v_min(:, i) = parts(i)%vx(:)
            rho_min(i) = parts(i)%rho
         end do

         !Interaction parameters, calculating neighboring particles
         call flink_list(maxinter, scale_k, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, niac, parts, pairs)

         do ki = 1, 4

            ! update halo particles after first increment
            if (ki > 1) then
               t_dist = t_dist - MPI_WTIME()
               call ORB_sendrecv_haloupdate(ki, MPI_types%haloupdatetype, n_process_neighbour, neighbours)
               t_dist = t_dist + MPI_WTIME()
            end if

            ! update pressure of newly updated real and halo particles
            do i = 1, ntotal_loc + nhalo_loc
               parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1_f)/gamma
            end do

            ! update ghost particles
            if (ki > 1) call update_ghost_part(ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, parts)

            ! calculating forces
            call single_step(ki, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, parts, niac, pairs, dvxdt(:, :, ki), drho(:, ki))

            call RK4_update(ki, ntotal_loc, v_min, rho_min, dvxdt, drho, parts)

         end do

         time = time + dt

         cputime = cputime + MPI_WTIME()

         output_time = output_time - MPI_WTIME()

         ! write output data
         if (mod(itimestep, save_step) .eq. 0) then
            call output(procid, numprocs)
         end if

         output_time = output_time + MPI_WTIME()

         if (mod(itimestep, print_step) .eq. 0) then
            call print_loadbalance(procid, numprocs)
         end if

      end do

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
