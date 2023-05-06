module time_integration_m

   use datatypes, only: particles, interactions, time_tracking, system_clock_timer
   use flink_list_m, only: flink_list
   ! use input_m, only: update_virt_part
   use ORB_m, only: ORB, neighbours, n_process_neighbour
!    use ORB_sr_m, only: ORB_sendrecv_haloupdate
   use param, only: f, tf, dim, dt, rh0, c, gamma, g
   use single_step_m, only: single_step
   use summary_m, only: print_loadbalance
   use output_m, only: output

   private
   public:: time_integration

contains

   subroutine time_integration(maxtimestep, print_step, save_step, my_rank, num_ranks, maxinter, timings, parts, &
                               pairs, nexti)

      integer, intent(in):: maxtimestep, print_step, save_step, my_rank, num_ranks, maxinter
      type(time_tracking), intent(inout):: timings
      integer, intent(inout):: nexti(:)
      type(particles), intent(inout):: parts
      type(interactions), intent(inout):: pairs(maxinter)
      integer:: i, j, k, d, n, itimestep, niac
      real(f):: time
      real(tf):: tmptime
      real(f), allocatable:: dvxdt(:, :), drhodt(:), vw(:)

      allocate (dvxdt(dim, parts%maxn), drhodt(parts%maxn), vw(parts%maxn))

      ! initializing
      time = 0._f
      drhodt(:) = 0._f
      dvxdt(1:dim - 1, :) = 0._f
      dvxdt(dim, :) = -g

      timings%t_wall = timings%t_wall - system_clock_timer()

      ! Time-integration (Leap-Frog)
      do itimestep = 1, maxtimestep

         ! save properties at start of step, update properties to mid-step.
         do i = 1, parts%ntotal_loc+parts%nvirt_loc
            if (parts%itype(i)==1) then
               do d = 1, dim
                  parts%v_min(d, i) = parts%vx(d, i)
                  parts%vx(d, i) = parts%vx(d, i) + 0.5_f*dt*dvxdt(d, i)
               end do
               parts%rho_min(i) = parts%rho(i)
               parts%rho(i) = parts%rho(i) + 0.5_f*dt*drhodt(i)
            end if
         end do

         ! distributing particles
         call ORB(itimestep, my_rank, num_ranks, parts, timings)

         ! Finding neighbours within kh
         call flink_list(maxinter, niac, parts, pairs, nexti)

         ! call update_virt_part(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, nexti, vw)

         ! update pressure of newly updated real and halo particles
         do i = 1, parts%ntotal_loc + parts%nhalo_loc + parts%nvirt_loc
            parts%p(i) = rh0*c**2*((parts%rho(i)/rh0)**gamma - 1._f)/gamma
         end do

         ! calculating forces
         call single_step(parts, niac, pairs, dvxdt, drhodt, nexti)

         ! updating positions and velocity to full timestep
         do i = 1, parts%ntotal_loc+parts%nvirt_loc
            if (parts%itype(i)==1) then
               parts%rho(i) = parts%rho_min(i) + dt*drhodt(i)
               do d = 1, dim
                  parts%vx(d, i) = parts%v_min(d, i) + dt*dvxdt(d, i)
                  parts%x(d, i) = parts%x(d, i) + dt*parts%vx(d, i)
               end do
            end if
         end do

         time = time + dt

         timings%t_output = timings%t_output - system_clock_timer()

         ! write output data
         if (mod(itimestep, save_step) .eq. 0) then
            call output(itimestep, save_step, my_rank, num_ranks, parts)
         end if

         if (mod(itimestep, print_step) .eq. 0) then
            tmptime = timings%t_wall + system_clock_timer()
            call print_loadbalance(my_rank, num_ranks, tmptime, parts%ntotal_loc, parts%nhalo_loc, parts%nvirt_loc, &
               niac, itimestep, time, maxtimestep)
         end if

         timings%t_output = timings%t_output + system_clock_timer()

      end do

      timings%t_wall = timings%t_wall + system_clock_timer()

   end subroutine time_integration

end module time_integration_m
