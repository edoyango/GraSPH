module time_integration_m

   use datatypes, only: particles, interactions, time_tracking, system_clock_timer
   use flink_list_m, only: flink_list
   use input_m, only: update_virt_part
   use ORB_m, only: ORB, neighbours, n_process_neighbour
!    use ORB_sr_m, only: ORB_sendrecv_haloupdate
   use param, only: f, tf, dim, dt, rh0, c, gamma, g
   use single_step_m, only: single_step
   use summary_m, only: print_loadbalance
#ifdef PARALLEL
   use output_m, only: output_parallel
#else
   use output_m, only: output_serial
#endif

   private
   public:: time_integration

contains

   subroutine time_integration(maxtimestep, print_step, save_step, thisImage, numImages, maxnloc, maxinter, timings, &
                               scale_k, ntotal_loc, nvirt_loc, nhalo_loc, ntotal, nvirt, parts, pairs, nexti)

      integer, intent(in):: maxtimestep, print_step, save_step, thisImage, numImages, maxnloc, maxinter, ntotal, nvirt
      real(f), intent(in):: scale_k
      type(time_tracking), intent(inout):: timings
      integer, intent(inout):: nexti(:)
      integer, codimension[*], intent(inout):: ntotal_loc, nvirt_loc, nhalo_loc
      type(particles), codimension[*], intent(inout):: parts(maxnloc)
      type(interactions), intent(inout):: pairs(maxinter)
      integer:: i, j, k, d, n, itimestep, niac
      real(f):: time
      real(tf):: tmptime
      real(f), allocatable:: dvxdt(:, :), drhodt(:), vw(:)

      allocate (dvxdt(dim, maxnloc), drhodt(maxnloc), vw(maxnloc))

      ! initializing
      time = 0._f
      drhodt(:) = 0._f
      dvxdt(1:dim - 1, :) = 0._f
      dvxdt(dim, :) = -g

      timings%t_wall = timings%t_wall - system_clock_timer()

      ! Time-integration (Leap-Frog)
      do itimestep = 1, maxtimestep

         ! save properties at start of step, update properties to mid-step.
         do i = 1, ntotal_loc+nvirt_loc
            if (parts(i)%itype==1) then
               parts(i)%v_min(:) = parts(i)%vx(:)
               parts(i)%vx(:) = parts(i)%vx(:) + 0.5_f*dt*dvxdt(:, i)
               parts(i)%rho_min = parts(i)%rho
               parts(i)%rho = parts(i)%rho + 0.5_f*dt*drhodt(i)
            end if
         end do

         ! distributing particles
         call ORB(itimestep, thisImage, numImages, scale_k, ntotal, ntotal_loc, nvirt, nvirt_loc, nhalo_loc, parts, &
            timings)

         ! Finding neighbours within kh
         call flink_list(maxinter, scale_k, ntotal_loc, nhalo_loc, nvirt_loc, niac, parts, pairs, nexti)

         call update_virt_part(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, nexti, vw)

         ! update pressure of newly updated real and halo particles
         do i = 1, ntotal_loc + nhalo_loc + nvirt_loc
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1._f)/gamma
         end do

         ! calculating forces
         call single_step(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, dvxdt, drhodt, nexti)

         ! updating positions and velocity to full timestep
         do i = 1, ntotal_loc+nvirt_loc
            if (parts(i)%itype==1) then
               parts(i)%rho = parts(i)%rho_min + dt*drhodt(i)
               parts(i)%vx(:) = parts(i)%v_min(:) + dt*dvxdt(:, i)
               parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)
            end if
         end do

         time = time + dt

         timings%t_output = timings%t_output - system_clock_timer()

         ! write output data
         if (mod(itimestep, save_step) .eq. 0) then
#ifdef PARALLEL
            call output_parallel(itimestep, save_step, thisImage, numImages, ntotal_loc, nhalo_loc, nvirt_loc, parts, ntotal)
#else
            call output_serial(itimestep, save_step, ntotal, nvirt, parts)
#endif
         end if

         if (mod(itimestep, print_step) .eq. 0) then
            tmptime = timings%t_wall + system_clock_timer()
            call print_loadbalance(thisImage, numImages, tmptime, ntotal_loc, nhalo_loc, nvirt_loc, niac, itimestep, &
                                   time, maxtimestep)
         end if

         timings%t_output = timings%t_output + system_clock_timer()

      end do

      timings%t_wall = timings%t_wall + system_clock_timer()

   end subroutine time_integration

end module time_integration_m
