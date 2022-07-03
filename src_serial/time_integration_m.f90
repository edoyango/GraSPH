module time_integration_m
   
   use datatypes, only: particles, interactions, time_tracking, system_clock_timer
   use param, only: f, dim, rh0, gamma, c, dt

   use flink_list_m, only: flink_list
   use input_m, only: generate_ghost_part, update_ghost_part
   use output_m, only: output
   use single_step_m, only: single_step
   use summary_m, only: print_update

   public:: time_integration

contains

   !==============================================================================================================================
   subroutine time_integration(time, timings, ntotal, nvirt, nghos, parts, print_step, save_step, maxtimestep, niac, pairs, nexti, &
                               scale_k, gind)
      ! Subroutine responsible for the main time-integration loop

      implicit none
      integer, intent(in):: ntotal, nvirt, print_step, save_step, maxtimestep
      real(f), intent(in):: scale_k
      integer, intent(inout):: niac, nghos, nexti(:), gind(:)
      real(f), intent(inout):: time
      type(time_tracking), intent(inout):: timings
      type(particles), intent(inout):: parts(:)
      type(interactions), intent(inout):: pairs(:)
      integer:: i, itimestep, ki, maxn
      real(f), allocatable:: v_min(:, :), rho_min(:), dvxdt(:, :, :), drhodt(:, :)
      
      maxn = size(parts)
      
      allocate (v_min(dim, ntotal), rho_min(ntotal), dvxdt(dim, maxn, 4), drhodt(maxn, 4))

      ! Time-integration (Leap-Frog)
      do itimestep = 1, maxtimestep
         
         timings%t_compute = timings%t_compute - system_clock_timer()
         
         ! Update density and velocity half a time step (not at first time-step)
         do i = 1, ntotal
            v_min(:, i) = parts(i)%vx(:)
            rho_min(i) = parts(i)%rho
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1_f)/gamma
         end do

         call generate_ghost_part(scale_k, ntotal, nvirt, nghos, parts, gind)

         !Interaction parameters, calculating neighboring particles
         call flink_list(scale_k, ntotal, nvirt, nghos, parts, niac, pairs, nexti)

         do ki = 1, 4

            ! calculating rate of change of speed and density on particles
            call single_step(ki, ntotal, nvirt, nghos, niac, pairs, parts, dvxdt(:, :, ki), drhodt(:, ki), nexti)

            ! applying update to particles
            call RK4_update(ki, ntotal, v_min, rho_min, dvxdt, drhodt, parts)

            ! updating ghost particles to reflect real particles
            call update_ghost_part(ntotal, nvirt, nghos, gind, parts)

         end do

         time = time + dt
         
         timings%t_compute = timings%t_compute + system_clock_timer()
         timings%t_output = timings%t_output - system_clock_timer()

         ! write output data
         if (mod(itimestep, save_step) .eq. 0) then
            call output(itimestep, save_step, ntotal, nvirt, nghos, parts)
         end if

         timings%t_output = timings%t_output + system_clock_timer()
         
         if (mod(itimestep, print_step) .eq. 0) then
            call print_update(itimestep, maxtimestep, ntotal, nvirt, nghos, niac, parts, pairs, nexti, time, timings)
         end if

      end do

   end subroutine time_integration

   pure subroutine RK4_update(ki, ntotal, v_min, rho_min, dvxdt, drhodt, parts)

      implicit none
      integer, intent(in):: ki, ntotal
      real(f), intent(in):: v_min(:, :), rho_min(:), dvxdt(:, :, :), drhodt(:, :)
      type(particles), intent(inout):: parts(:)
      integer:: i

      select case (ki)
      case default
         do i = 1, ntotal
            parts(i)%vx(:) = v_min(:, i) + 0.5_f*dt*dvxdt(:, i, ki)
            parts(i)%rho = rho_min(i) + 0.5_f*dt*drhodt(i, ki)
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1._f)/gamma
         end do
      case (3)
         do i = 1, ntotal
            parts(i)%vx(:) = v_min(:, i) + dt*dvxdt(:, i, 3)
            parts(i)%rho = rho_min(i) + dt*drhodt(i, 3)
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1._f)/gamma
         end do
      case (4)
         do i = 1, ntotal
            parts(i)%vx(:) = v_min(:, i) + &
                             dt/6_f*(dvxdt(:, i, 1) + 2._f*dvxdt(:, i, 2) + 2._f*dvxdt(:, i, 3) + dvxdt(:, i, 4))
            parts(i)%rho = rho_min(i) + dt/6_f*(drhodt(i, 1) + 2._f*drhodt(i, 2) + 2._f*drhodt(i, 3) + drhodt(i, 4))
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1._f)/gamma
            parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)
         end do
      end select

   end subroutine RK4_update

end module time_integration_m
