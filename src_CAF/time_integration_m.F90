module time_integration_m

   use datatypes, only: particles, interactions, time_tracking, system_clock_timer
   use flink_list_m, only: flink_list
   use input_m, only: update_virt_part
   use iso_fortran_env
   use ORB_m, only: ORB, neighbours, n_process_neighbour
!    use ORB_sr_m, only: ORB_sendrecv_haloupdate
   use param, only: f, tf, dt, rh0, c, gamma, g, hsml, irho, mass
   use single_step_m, only: single_step
   use summary_m, only: print_loadbalance
   use omp_lib
   use output_m, only: output
   use kernel_m, only: kernel

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
      real(f):: mingridx(dim), maxgridx(dim), dcell, dx(dim), r, tdwdx(dim), tw, factor, q, flttmp, vxj(dim)
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

      call omp_set_num_threads(2)

      !$omp parallel default(shared) &
      !$omp private(icell, jcell, kcell, inttmp, jth, dx, r, factor, q, xi, yi, zi, i, j, k, flttmp, vxj)

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
                        pairs(inttmp)%i = i
                        pairs(inttmp)%j = jth
                        call kernel(r, dx, hsml, pairs(inttmp)%w, pairs(inttmp)%dwdx)
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
                           pairs(inttmp)%i = i
                           pairs(inttmp)%j = jth
                           call kernel(r, dx, hsml, pairs(inttmp)%w, pairs(inttmp)%dwdx)
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
         deallocate(pincell, cells)
         !$omp end single

         ! call update_virt_part(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, nexti, vw)
         !$omp do
         do i = 1, ntotal_loc+nvirt_loc
            if (parts(i)%itype < 0) then
               vw(i) = 0._f
               parts(i)%rho = 0._f
               parts(i)%vx(:) = 0._f
            end if
         end do
         !$omp end do

         !$omp do
         do k = 1, niac
            i = pairs(k)%i
            j = pairs(k)%j
            if ( parts(i)%itype < 0 .and. parts(j)%itype > 0) then
               flttmp = mass*pairs(k)%w/parts(j)%rho
               !$omp atomic
               vw(i) = vw(i) + flttmp
               !$omp atomic
               parts(i)%rho = parts(i)%rho + mass*pairs(k)%w
               vxj = parts(j)%vx(:)*flttmp
               select case (parts(i)%itype)
               case default
                  !$omp atomic
                  parts(i)%vx(1) = parts(i)%vx(1) - vxj(1)
                  !$omp atomic
                  parts(i)%vx(2) = parts(i)%vx(2) - vxj(2)
                  !$omp atomic
                  parts(i)%vx(3) = parts(i)%vx(3) - vxj(3)
               case (-2)
                  !$omp atomic
                  parts(i)%vx(1) = parts(i)%vx(1) + vxj(1)
                  !$omp atomic
                  parts(i)%vx(2) = parts(i)%vx(2) + vxj(2)
                  !$omp atomic
                  parts(i)%vx(3) = parts(i)%vx(3) - vxj(3)
               case (-3)
                  !$omp atomic
                  parts(i)%vx(1) = parts(i)%vx(1) - vxj(1)
                  !$omp atomic
                  parts(i)%vx(2) = parts(i)%vx(2) + vxj(2)
                  !$omp atomic
                  parts(i)%vx(3) = parts(i)%vx(3) + vxj(3)
               case (-4)
                  !$omp atomic
                  parts(i)%vx(1) = parts(i)%vx(1) + vxj(1)
                  !$omp atomic
                  parts(i)%vx(2) = parts(i)%vx(2) - vxj(2)
                  !$omp atomic
                  parts(i)%vx(3) = parts(i)%vx(3) + vxj(3)
               case (-5)
                  !$omp atomic
                  parts(i)%vx(1) = parts(i)%vx(1) - vxj(1)
                  !$omp atomic
                  parts(i)%vx(2) = parts(i)%vx(2) - vxj(2)
                  !$omp atomic
                  parts(i)%vx(3) = parts(i)%vx(3) + vxj(3)
               end select
            else if (parts(j)%itype < 0 .and. parts(i)%itype > 0) then
               flttmp = mass*pairs(k)%w/parts(i)%rho
               !$omp atomic
               vw(j) = vw(j) + flttmp
               !$omp atomic
               parts(j)%rho = parts(j)%rho + mass*pairs(k)%w
               vxj = parts(i)%vx(:)*flttmp
               select case (parts(j)%itype)
               case default
                  !$omp atomic
                  parts(j)%vx(1) = parts(j)%vx(1) - vxj(1)
                  !$omp atomic
                  parts(j)%vx(2) = parts(j)%vx(2) - vxj(2)
                  !$omp atomic
                  parts(j)%vx(3) = parts(j)%vx(3) - vxj(3)
               case (-2)
                  !$omp atomic
                  parts(j)%vx(1) = parts(j)%vx(1) + vxj(1)
                  !$omp atomic
                  parts(j)%vx(2) = parts(j)%vx(2) + vxj(2)
                  !$omp atomic
                  parts(j)%vx(3) = parts(j)%vx(3) - vxj(3)
               case (-3)
                  !$omp atomic
                  parts(j)%vx(1) = parts(j)%vx(1) - vxj(1)
                  !$omp atomic
                  parts(j)%vx(2) = parts(j)%vx(2) + vxj(2)
                  !$omp atomic
                  parts(j)%vx(3) = parts(j)%vx(3) + vxj(3)
               case (-4)
                  !$omp atomic
                  parts(j)%vx(1) = parts(j)%vx(1) + vxj(1)
                  !$omp atomic
                  parts(j)%vx(2) = parts(j)%vx(2) - vxj(2)
                  !$omp atomic
                  parts(j)%vx(3) = parts(j)%vx(3) + vxj(3)
               case (-5)
                  !$omp atomic
                  parts(j)%vx(1) = parts(j)%vx(1) - vxj(1)
                  !$omp atomic
                  parts(j)%vx(2) = parts(j)%vx(2) - vxj(2)
                  !$omp atomic
                  parts(j)%vx(3) = parts(j)%vx(3) + vxj(3)
               end select
            end if
         end do
         !$omp end do

         !$omp do
         do i = 1, ntotal_loc+nvirt_loc
            if (parts(i)%itype < 0) then
               if (vw(i) > 0._f) then
                  parts(i)%rho = parts(i)%rho/vw(i)
                  parts(i)%vx(:) = parts(i)%vx(:)/vw(i)
               else
                  parts(i)%rho = irho
                  parts(i)%vx(:) = 0._f
               end if
            end if
         end do
         !$omp end do


         ! update pressure of newly updated real and halo particles
         !$omp do
         do i = 1, ntotal_loc + nvirt_loc
            parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma - 1._f)/gamma
         end do
         !$omp end do

         ! ! calculating forces
         ! call single_step(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, dvxdt, drhodt, nexti)

         ! updating positions and velocity to full timestep
         !$omp do
         do i = 1, ntotal_loc+nvirt_loc
            if (parts(i)%itype==1) then
               parts(i)%rho = parts(i)%rho_min + dt*drhodt(i)
               parts(i)%vx(:) = parts(i)%v_min(:) + dt*dvxdt(:, i)
               parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)
            end if
         end do
         !$omp end do

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
