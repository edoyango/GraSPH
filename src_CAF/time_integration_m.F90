module time_integration_m

   use datatypes, only: particles, interactions, time_tracking, system_clock_timer
   use flink_list_m, only: flink_list
   use input_m, only: update_virt_part
   use iso_fortran_env
   use ORB_m, only: ORB, neighbours, n_process_neighbour
!    use ORB_sr_m, only: ORB_sendrecv_haloupdate
   use param, only: f, tf, dt, rh0, c, gamma, g, hsml
   use single_step_m, only: single_step
   use summary_m, only: print_loadbalance
   use omp_lib
   use output_m, only: output

   private
   public:: time_integration

contains

   subroutine time_integration(maxtimestep, print_step, save_step, thisImage, numImages, maxnloc, maxinter, timings, &
                               scale_k, ntotal_loc, nvirt_loc, nhalo_loc, ntotal, nvirt, parts, pairs, nexti)

      use param, only: dim, pi

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

      ! linked list variabls
      real(f):: mingridx(dim), maxgridx(dim), dcell, dx(dim), r, tdwdx(dim), tw, factor, q
      integer:: ngridx(dim), icell, jcell, kcell, inttmp, jth, xi, yi, zi
      integer, allocatable:: pincell(:, :, :), cells(:, :, :, :), gridind(:, :)
      integer, parameter:: maxpcell = 125
      integer, parameter:: sweep(3, 13) = reshape((/-1, -1, -1, &
                                                    -1, -1, 0, &
                                                    -1, -1, 1, &
                                                    -1, 0, -1, &
                                                    -1, 0, 0, &
                                                    -1, 0, 1, &
                                                    -1, 1, -1, &
                                                    -1, 1, 0, &
                                                    -1, 1, 1, &
                                                    0, -1, -1, &
                                                    0, -1, 0, &
                                                    0, -1, 1, &
                                                    0, 0, -1/), (/3, 13/))

      allocate (dvxdt(dim, maxnloc), drhodt(maxnloc), vw(maxnloc))
      allocate (gridind(dim, maxnloc))

      ! initializing
      time = 0._f
      drhodt(:) = 0._f
      dvxdt(1:dim - 1, :) = 0._f
      dvxdt(dim, :) = -g

      timings%t_wall = timings%t_wall - system_clock_timer()

      call omp_set_num_threads(12)

      !$omp parallel default(shared) private(icell, jcell, kcell, inttmp, jth, dx, r, factor, q, xi, yi, zi)

      ! Time-integration (Leap-Frog)
      do itimestep = 1, maxtimestep

         ! save properties at start of step, update properties to mid-step.
         !$omp do
         do i = 1, ntotal_loc+nvirt_loc
            if (parts(i)%itype==1) then
               parts(i)%v_min(:) = parts(i)%vx(:)
               parts(i)%vx(:) = parts(i)%vx(:) + 0.5_f*dt*dvxdt(:, i)
               parts(i)%rho_min = parts(i)%rho
               parts(i)%rho = parts(i)%rho + 0.5_f*dt*drhodt(i)
            end if
         end do
         !$omp end do nowait

         ! ! distributing particles
         ! call ORB(itimestep, thisImage, numImages, scale_k, ntotal, ntotal_loc, nvirt, nvirt_loc, nhalo_loc, parts, &
         !    timings)

         ! ! Finding neighbours within kh
         ! call flink_list(maxinter, scale_k, ntotal_loc, nhalo_loc, nvirt_loc, niac, parts, pairs, nexti)
         !$omp do reduction(min:mingridx) reduction(max:maxgridx)
         do i = 1, ntotal_loc+nvirt_loc
            do d = 1, dim
               mingridx(d) = min(mingridx(d), parts(i)%x(d))
               maxgridx(d) = max(maxgridx(d), parts(i)%x(d))
            end do
         end do
         !$omp end do nowait

         !$omp single
         dcell = scale_k*hsml
         maxgridx(:) = maxgridx(:) + 2._f*dcell
         mingridx(:) = mingridx(:) - 2._f*dcell
         ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
         maxgridx(:) = mingridx(:) + ngridx(:)*dcell

         allocate (pincell(ngridx(1), ngridx(2), ngridx(3)), &
                cells(maxpcell, ngridx(1), ngridx(2), ngridx(3)))
         !$omp end single

         !$omp workshare
         pincell(:, :, :) = 0
         !$omp end workshare

         !$omp do
         do i = 1, ntotal_loc + nvirt_loc
            icell = int((parts(i)%x(1) - mingridx(1))/dcell) + 1
            jcell = int((parts(i)%x(2) - mingridx(2))/dcell) + 1
            kcell = int((parts(i)%x(3) - mingridx(3))/dcell) + 1
            gridind(1, i) = icell
            gridind(2, i) = jcell
            gridind(3, i) = kcell
            !$omp atomic capture
            pincell(icell, jcell, kcell) = pincell(icell, jcell, kcell) + 1
            inttmp = pincell(icell, jcell, kcell)
            !$omp end atomic
            
            cells(inttmp, icell, jcell, kcell) = i
         end do
         !$omp end do

         !$omp single
         niac = 0
         !$omp end single

         !$omp do
         do i = 1, ntotal_loc + nvirt_loc

            icell = gridind(1, i)
            jcell = gridind(2, i)
            kcell = gridind(3, i)

            do j = 1, pincell(icell, jcell, kcell)
               jth = cells(j, icell, jcell, kcell)
               if (jth > i .and. (parts(i)%itype > 0 .or. parts(jth)%itype > 0)) then
                  dx(:) = parts(i)%x(:) - parts(jth)%x(:)
                  r = sqrt(sum(dx*dx))
                  if (r < scale_k*hsml) then
                     !$omp atomic capture
                     niac = niac + 1
                     inttmp = niac
                     !$omp end atomic
                     if (inttmp < maxinter) then
                        pairs(inttmp)%j = jth
                        factor = 21._f/(256._f*pi*hsml*hsml*hsml)
                        q = r/hsml
                        pairs(inttmp)%w = factor*max(0._f, 2._f-q)**4*(2._f*q + 1._f)
                        pairs(inttmp)%dwdx = -factor*10._f*q*max(0._f, 2._f-q)**3*dx(:)/(r*hsml)
                     else
                        error stop
                     end if
                  end if
               end if
            end do

            ! finding pairs within cells adjacent to i's cell
            do k = 1, 13
               xi = icell + sweep(1, k)
               yi = jcell + sweep(2, k)
               zi = kcell + sweep(3, k)
               do j = 1, pincell(xi, yi, zi)
                  jth = cells(j, xi, yi, zi)
                  if (parts(i)%itype > 0 .or. parts(jth)%itype > 0) then
                     dx(:) = parts(i)%x(:) - parts(jth)%x(:)
                     r = sqrt(sum(dx*dx))
                     if (r < scale_k*hsml) then
                        !$omp atomic capture
                        niac = niac + 1
                        inttmp = niac
                        !$omp end atomic
                        if (inttmp < maxinter) then
                           pairs(inttmp)%j = jth
                           factor = 21._f/(256._f*pi*hsml*hsml*hsml)
                           q = r/hsml
                           pairs(inttmp)%w = factor*max(0._f, 2._f-q)**4*(2._f*q + 1._f)
                           pairs(inttmp)%dwdx = -factor*10._f*q*max(0._f, 2._f-q)**3*dx(:)/(r*hsml)
                        else
                           error stop
                        end if
                     end if
                  end if
               end do
            end do
         end do
         !$omp end do

            

         !$omp single
         write(*,*) pairs(1)%w, pairs(2)%w, pairs(2)%w, pairs(1)%dwdx !sum(pincell), cells(1:pincell(3, 3, 3), 3, 3, 3)
         deallocate(pincell, cells)
         !$omp end single

         ! call update_virt_part(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, nexti, vw)

         ! ! update pressure of newly updated real and halo particles
         ! do i = 1, ntotal_loc + nhalo_loc + nvirt_loc
         !    parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1._f)/gamma
         ! end do

         ! ! calculating forces
         ! call single_step(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, dvxdt, drhodt, nexti)

         ! ! updating positions and velocity to full timestep
         ! do i = 1, ntotal_loc+nvirt_loc
         !    if (parts(i)%itype==1) then
         !       parts(i)%rho = parts(i)%rho_min + dt*drhodt(i)
         !       parts(i)%vx(:) = parts(i)%v_min(:) + dt*dvxdt(:, i)
         !       parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)
         !    end if
         ! end do

         time = time + dt

         timings%t_output = timings%t_output - system_clock_timer()

         !$omp master
         ! write output data
!          if (mod(itimestep, save_step) .eq. 0) then
! ! #ifdef PARALLEL
!             call output(itimestep, save_step, thisImage, numImages, ntotal_loc, nhalo_loc, nvirt_loc, parts, ntotal)
! ! #else
! !             call output_serial(itimestep, save_step, ntotal, nvirt, parts)
! ! #endif
!          end if

         if (mod(itimestep, print_step) .eq. 0) then
            tmptime = timings%t_wall + system_clock_timer()
            call print_loadbalance(thisImage, numImages, tmptime, ntotal_loc, nhalo_loc, nvirt_loc, niac, itimestep, &
                                   time, maxtimestep)
         end if
         !$omp end master

         timings%t_output = timings%t_output + system_clock_timer()

      end do

      !$omp end parallel

      timings%t_wall = timings%t_wall + system_clock_timer()

   end subroutine time_integration

end module time_integration_m
