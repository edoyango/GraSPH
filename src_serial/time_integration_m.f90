module time_integration_m

   use datatypes, only: particles, interactions
   use param, only: f, dim, rh0, gamma, c, dt

   use flink_list_m, only: flink_list
   use input_m, only: generate_ghost_part, update_ghost_part, update_virt_part
   use output_m, only: output
   use single_step_m, only: single_step
   use summary_m, only: print_update

   public:: time_integration

contains

   !==============================================================================================================================
   subroutine time_integration(time, cputime, output_time, test_time, ntotal, nvirt, nghos, parts, print_step, save_step, &
                               maxtimestep, niac, pairs, scale_k, gind)
      ! Subroutine responsible for the main time-integration loop

      implicit none
      integer, intent(in):: ntotal, nvirt, print_step, save_step, maxtimestep
      real(f), intent(in):: scale_k
      integer, intent(inout):: niac, nghos, gind(:)
      real(f), intent(inout):: time, cputime, output_time, test_time
      type(particles), intent(inout):: parts(:)
      type(interactions), intent(inout):: pairs(:)
      integer:: i, itimestep
      real(f):: t1, t2, t3, t4
      real(f), allocatable:: v_min(:, :), rho_min(:), dvxdt(:, :, :), drho(:, :)

      allocate (v_min(dim, ntotal), rho_min(ntotal), dvxdt(dim, ntotal, 4), drho(ntotal, 4))

      call CPU_TIME(t1)

      ! Time-integration (Leap-Frog)
      do itimestep = 1, maxtimestep

         ! Update density and velocity half a time step (not at first time-step)
         do i = 1, ntotal
            v_min(:, i) = parts(i)%vx(:)
            rho_min(i) = parts(i)%rho
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1_f)/gamma
         end do

         call generate_ghost_part(scale_k, ntotal, nvirt, nghos, parts, gind)

         !Interaction parameters, calculating neighboring particles
         call flink_list(scale_k, ntotal, nvirt, nghos, parts, niac, pairs)

!~                call update_virt_part

         ! calculating forces (k1)
         call single_step(1, ntotal, nvirt, nghos, niac, pairs, parts, dvxdt(:, :, 1), drho(:, 1))

         ! updating data for mid-timestep base on k1
         do i = 1, ntotal
            parts(i)%vx(:) = v_min(:, i) + 0.5_f*dt*dvxdt(:, i, 1)
            parts(i)%rho = rho_min(i) + 0.5_f*dt*drho(i, 1)
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1_f)/gamma
         end do

         call update_ghost_part(ntotal, nvirt, nghos, gind, parts)

!~                call update_virt_part

         ! calculating forces (k2)
         call single_step(2, ntotal, nvirt, nghos, niac, pairs, parts, dvxdt(:, :, 2), drho(:, 2))

         ! updating data for mid-timestep base on k2
         do i = 1, ntotal
            parts(i)%vx(:) = v_min(:, i) + 0.5_f*dt*dvxdt(:, i, 2)
            parts(i)%rho = rho_min(i) + 0.5_f*dt*drho(i, 2)
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1_f)/gamma
         end do

         call update_ghost_part(ntotal, nvirt, nghos, gind, parts)

!~                call update_virt_part

         ! calculating forces (k3)
         call single_step(3, ntotal, nvirt, nghos, niac, pairs, parts, dvxdt(:, :, 3), drho(:, 3))

         ! updating data for mid-timestep base on k3
         do i = 1, ntotal
            parts(i)%vx(:) = v_min(:, i) + dt*dvxdt(:, i, 3)
            parts(i)%rho = rho_min(i) + dt*drho(i, 3)
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1_f)/gamma
         end do

         call update_ghost_part(ntotal, nvirt, nghos, gind, parts)

!~                call update_virt_part

         call single_step(4, ntotal, nvirt, nghos, niac, pairs, parts, dvxdt(:, :, 4), drho(:, 4))

         ! updating data for mid-timestep base on k1, k2, k3, k4
         do i = 1, ntotal
            parts(i)%vx(:) = v_min(:, i) + dt/6_f*(dvxdt(:, i, 1) + 2_f*dvxdt(:, i, 2) + 2_f*dvxdt(:, i, 3) + dvxdt(:, i, 4))

            parts(i)%rho = rho_min(i) + dt/6_f*(drho(i, 1) + 2_f*drho(i, 2) + 2_f*drho(i, 3) + drho(i, 4))

            parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)

            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1_f)/gamma
         end do

         call update_ghost_part(ntotal, nvirt, nghos, gind, parts)

!~                call update_virt_part

         time = time + dt

         call CPU_TIME(t3)

         ! write output data
         if (mod(itimestep, save_step) .eq. 0) then
            call output(itimestep, save_step, ntotal, nvirt, nghos, parts)
         end if

         call CPU_TIME(t4)
         output_time = output_time + t4 - t3

         if (mod(itimestep, print_step) .eq. 0) then
            call CPU_TIME(t2)
            cputime = t2 - t1
            call print_update(itimestep, maxtimestep, ntotal, nvirt, nghos, niac, parts, pairs, time, cputime, output_time)
         end if

      end do

      call CPU_TIME(t2)
      cputime = t2 - t1

   end subroutine time_integration

end module time_integration_m
