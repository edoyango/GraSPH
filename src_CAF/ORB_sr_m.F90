module ORB_sr_m

   use datatypes, only: particles, system_clock_timer
   use iso_fortran_env, only: lock_type, event_type
#ifdef PARALLEL
   use mpi
#endif
   use param_para, only: neighbour_data
   use param, only: tf, f, dim, hsml, halotype, scale_k

   private
   public:: ORB_sendrecv_diffuse, ORB_sendrecv_halo

contains

   !===================================================================================================================
   subroutine ORB_sendrecv_diffuse(itimestep, my_rank, my_bounds, repartition_mode, n_process_neighbour, &
                                   neighbours, old_ntotal_loc, ntotal_loc, parts)
      ! Recursive function to exchange physical particles. In cases were subdomain boundaries are updated, the possibility of needing
      ! diffusion is considered

#ifdef PARALLEL

      implicit none
      integer, intent(in):: itimestep, my_rank, repartition_mode, n_process_neighbour
      real(f), intent(in):: my_bounds(2*dim)
      integer, intent(inout):: ntotal_loc
      type(particles), intent(inout):: parts
      type(neighbour_data), intent(inout):: neighbours(:)
      integer, intent(out):: old_ntotal_loc
      integer:: i, j, k, d, n, searchrange(2), entrydepth, nphys_send_all, ndiffuse, diff_dest, displ0, displ1, &
                neighbourImageIDs(n_process_neighbour), displ, ierr, request(2*8*n_process_neighbour), &
                status(MPI_STATUS_SIZE, 2*8*n_process_neighbour)
      real(f):: xmin_loc(dim), xmax_loc(dim), xi(dim), xmin_rem(dim), xmax_rem(dim), dr_min, dr
      real(tf):: tmptime
      logical:: diffuse
      integer, allocatable:: removal_list(:)

      ! Initialization
      diffuse = .true.
      entrydepth = 0
      searchrange(1) = 1
      searchrange(2) = ntotal_loc
      xmin_loc(:) = my_bounds(1:dim)
      xmax_loc(:) = my_bounds(dim + 1:2*dim)

      diffuseloop: do while (diffuse)

         ! Searching particles to remove within indices of searchrange(1) and searchrange(2), inclusive.
         ! At node 0, searchrange(1:2) = [1,ntotal_loc)
         allocate (removal_list(searchrange(2) - searchrange(1) + 2))
         nphys_send_all = 0

         ! If there are any particles that do not belong to the host process,
         ! begin searching for neighbouring processes to send the particle to.
         ! If particle is not contained with subdomain boundaries, send the particle to nearest process neighbouring current host.
         do i = 1, n_process_neighbour
            neighbours(i)%nphys_send = 0
            neighbours(i)%PhysPackSend%maxn = ntotal_loc
            call neighbours(i)%PhysPackSend%allocate_particles()
         end do

         tmptime = -system_clock_timer()

         ndiffuse = 0
         loop_through_parts: do i = searchrange(1), searchrange(2)
            xi(:) = parts%x(:, i)
            if (any(xi(:) < xmin_loc(:)) .or. any(xi(:) >= xmax_loc(:))) then
               nphys_send_all = nphys_send_all + 1
               removal_list(nphys_send_all) = i
               do n = 1, n_process_neighbour
                  xmin_rem(:) = neighbours(n)%bounds(1:dim)
                  xmax_rem(:) = neighbours(n)%bounds(dim + 1:2*dim)
                  if (.not. (any(xi(:) < xmin_rem(:)) .or. any(xi(:) >= xmax_rem(:)))) then
                     neighbours(n)%nphys_send = neighbours(n)%nphys_send + 1
                     ! neighbours(n)%PhysPackSend(neighbours(n)%nphys_send) = parts(i)
                     neighbours(n)%PhysPackSend%itype(neighbours(n)%nphys_send) = parts%itype(i)
                     neighbours(n)%PhysPackSend%indglob(neighbours(n)%nphys_send) = parts%indglob(i)
                     neighbours(n)%PhysPackSend%rho(neighbours(n)%nphys_send) = parts%rho(i)
                     neighbours(n)%PhysPackSend%rho_min(neighbours(n)%nphys_send) = parts%rho_min(i)
                     neighbours(n)%PhysPackSend%p(neighbours(n)%nphys_send) = parts%p(i)
                     neighbours(n)%PhysPackSend%x(:, neighbours(n)%nphys_send) = parts%x(:, i)
                     neighbours(n)%PhysPackSend%vx(:, neighbours(n)%nphys_send) = parts%vx(:, i)
                     neighbours(n)%PhysPackSend%v_min(:, neighbours(n)%nphys_send) = parts%v_min(:, i)
                     cycle loop_through_parts
                  end if
               end do

               ! if particle makes it past the do loop above, then it belongs to an image that is not the current
               ! image's neighbour
               ndiffuse = ndiffuse + 1
               dr_min = huge(1._f)
               diff_dest = 0
               do n = 1, n_process_neighbour
                  dr = 0._f
                  do d = 1, dim
                     dr = dr + MAX(0._f, &
                                   neighbours(n)%bounds(d) - xi(d), &
                                   xi(d) - neighbours(n)%bounds(dim + d))**2
                  end do
                  if (dr < dr_min) then
                     diff_dest = n
                     dr_min = dr
                  end if
               end do

               neighbours(diff_dest)%nphys_send = neighbours(diff_dest)%nphys_send + 1
               ! neighbours(diff_dest)%PhysPackSend(neighbours(diff_dest)%nphys_send) = parts(i)
               neighbours(diff_dest)%PhysPackSend%itype(neighbours(diff_dest)%nphys_send) = parts%itype(i)
               neighbours(diff_dest)%PhysPackSend%indglob(neighbours(diff_dest)%nphys_send) = parts%indglob(i)
               neighbours(diff_dest)%PhysPackSend%rho(neighbours(diff_dest)%nphys_send) = parts%rho(i)
               neighbours(diff_dest)%PhysPackSend%rho_min(neighbours(diff_dest)%nphys_send) = parts%rho_min(i)
               neighbours(diff_dest)%PhysPackSend%p(neighbours(diff_dest)%nphys_send) = parts%p(i)
               neighbours(diff_dest)%PhysPackSend%x(:, neighbours(diff_dest)%nphys_send) = parts%x(:, i)
               neighbours(diff_dest)%PhysPackSend%vx(:, neighbours(diff_dest)%nphys_send) = parts%vx(:, i)
               neighbours(diff_dest)%PhysPackSend%v_min(:, neighbours(diff_dest)%nphys_send) = parts%v_min(:, i)
            end if

         end do loop_through_parts

         do i = 1, n_process_neighbour
            call MPI_Irecv(neighbours(i)%nphys_recv, 1, MPI_INTEGER, neighbours(i)%rank, 0, MPI_COMM_WORLD, &
               request(2*i-1), ierr)
            call MPI_Isend(neighbours(i)%nphys_send, 1, MPI_INTEGER, neighbours(i)%rank, 0, MPI_COMM_WORLD, &
               request(2*i), ierr)
         end do

         removal_list(nphys_send_all+1) = 0

         ! Shifting information to remove sent particles
         n = 0
         if (nphys_send_all > 0 .and. nphys_send_all < ntotal_loc) then
            do i = removal_list(1), ntotal_loc
               if (i == removal_list(n + 1)) then
                  n = n + 1
               else
                  ! parts(i - n) = parts(i)
                  parts%itype(i-n) = parts%itype(i)
                  parts%indglob(i-n) = parts%indglob(i)
                  parts%rho(i-n) = parts%rho(i)
                  parts%rho_min(i-n) = parts%rho_min(i)
                  parts%p(i-n) = parts%p(i)
                  parts%x(:, i-n) = parts%x(:, i)
                  parts%vx(:, i-n) = parts%vx(:, i)
                  parts%v_min(:, i-n) = parts%v_min(:, i)
               end if
            end do
         end if

         ntotal_loc = ntotal_loc - nphys_send_all
         old_ntotal_loc = ntotal_loc
         call MPI_Waitall(2*n_process_neighbour, request, status, ierr)

         displ = ntotal_loc
         do n = 1, n_process_neighbour

            if (neighbours(n)%nphys_recv > 0) then
               call MPI_Irecv(parts%itype(displ+1), neighbours(n)%nphys_recv, MPI_INTEGER, neighbours(n)%rank, 0, &
                  MPI_COMM_WORLD, request(2*8*(n-1)+1), ierr)
               call MPI_Irecv(parts%indglob(displ+1), neighbours(n)%nphys_recv, MPI_INTEGER, neighbours(n)%rank, 0, &
                  MPI_COMM_WORLD, request(2*8*(n-1)+2), ierr)
               call MPI_Irecv(parts%rho(displ+1), neighbours(n)%nphys_recv, MPI_DOUBLE_PRECISION, neighbours(n)%rank, &
                  0, MPI_COMM_WORLD, request(2*8*(n-1)+3), ierr)
               call MPI_Irecv(parts%rho_min(displ+1), neighbours(n)%nphys_recv, MPI_DOUBLE_PRECISION, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+4), ierr)
               call MPI_Irecv(parts%p(displ+1), neighbours(n)%nphys_recv, MPI_DOUBLE_PRECISION, neighbours(n)%rank, &
                  0, MPI_COMM_WORLD, request(2*8*(n-1)+5), ierr)
               call MPI_Irecv(parts%x(1, displ+1), dim*neighbours(n)%nphys_recv, MPI_DOUBLE_PRECISION, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+6), ierr)
               call MPI_Irecv(parts%vx(1, displ+1), dim*neighbours(n)%nphys_recv, MPI_DOUBLE_PRECISION, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+7), ierr)
               call MPI_Irecv(parts%v_min(1, displ+1), dim*neighbours(n)%nphys_recv, MPI_DOUBLE_PRECISION, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+8), ierr)
               displ = displ + neighbours(n)%nphys_recv
            else
               request(2*8*(n-1)+1:2*8*(n-1)+8) = MPI_REQUEST_NULL
            end if
            
            if (neighbours(n)%nphys_send > 0) then
               call MPI_Isend(neighbours(n)%PhysPackSend%itype(1), neighbours(n)%nphys_send, MPI_INTEGER, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+9), ierr)
               call MPI_Isend(neighbours(n)%PhysPackSend%indglob(1), neighbours(n)%nphys_send, MPI_INTEGER, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+10), ierr)
               call MPI_Isend(neighbours(n)%PhysPackSend%rho(1), neighbours(n)%nphys_send, MPI_DOUBLE_PRECISION, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+11), ierr)
               call MPI_Isend(neighbours(n)%PhysPackSend%rho_min(1), neighbours(n)%nphys_send, MPI_DOUBLE_PRECISION, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+12), ierr)
               call MPI_Isend(neighbours(n)%PhysPackSend%p(1), neighbours(n)%nphys_send, MPI_DOUBLE_PRECISION, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+13), ierr)
               call MPI_Isend(neighbours(n)%PhysPackSend%x(1, 1), dim*neighbours(n)%nphys_send, MPI_DOUBLE_PRECISION, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+14), ierr)
               call MPI_Isend(neighbours(n)%PhysPackSend%vx(1, 1), dim*neighbours(n)%nphys_send, MPI_DOUBLE_PRECISION, &
                  neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+15), ierr)
               call MPI_Isend(neighbours(n)%PhysPackSend%v_min(1, 1), dim*neighbours(n)%nphys_send, &
                  MPI_DOUBLE_PRECISION, neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*8*(n-1)+16), ierr)
            else
               request(2*8*(n-1)+9:2*8*n) = MPI_REQUEST_NULL
            end if
         end do
         ntotal_loc = displ

         call MPI_Waitall(2*8*n_process_neighbour, request, status, ierr)

         ! if subdomain boundary update has occurred, check if diffusion is necessary.
         ! Perform if necessary
         if (repartition_mode >= 2) then

            if (my_rank == 1) write (*, '(6x, A)') 'Checking whether diffusion is needed...'

            ! call co_sum(ndiffuse)
            call MPI_Allreduce(MPI_IN_PLACE, ndiffuse, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

            if (ndiffuse == 0) then
               diffuse = .false.
               if (my_rank == 0) then
                  write (*, '(6x, A)') 'No diffusion needed. Continuing...'
                  write (*, '(A)') "_______________________________________________________________________________"
               end if
            else
               if (my_rank .eq. 0) then
                  write (*, '(6x, A,I0)') 'Diffusion occuring... Current timestep:         ', itimestep
                  write (*, '(6x, A,I0)') '                      Current depth:            ', entrydepth
                  write (*, '(6x, A,I0)') '                      Particles to be diffused: ', ndiffuse
               end if

               entrydepth = entrydepth + 1

               searchrange(1) = old_ntotal_loc + 1
               searchrange(2) = ntotal_loc

            end if
         else
            diffuse = .false.
         end if

         deallocate (removal_list)

         do i = 1, n_process_neighbour
            call neighbours(i)%PhysPackSend%deallocate_particles()
         end do

      end do diffuseloop

#endif

   end subroutine ORB_sendrecv_diffuse

   !====================================================================================================================
   subroutine ORB_sendrecv_halo(my_rank, my_bounds, n_process_neighbour, neighbours, old_ntotal_loc, ntotal_loc, parts)

      !subroutine responsible for sending sending halo particle information between processes, given
      !predetermined subdomain boundaires.
      !Note: subdomain boundaries are used as inputs (bounds_glob).

#ifdef PARALLEL

      implicit none
      integer, intent(in):: my_rank, n_process_neighbour, old_ntotal_loc, ntotal_loc
      real(f), intent(in):: my_bounds(2*dim)
      type(neighbour_data), intent(inout):: neighbours(n_process_neighbour)
      type(particles), intent(inout):: parts
      integer:: i, j, k, d, n, searchrange(2), displ, ierr, &
         request(2*5*n_process_neighbour), status(MPI_STATUS_SIZE, 2*5*n_process_neighbour)
      real(f):: xmin_loc(dim), xmax_loc(dim), xi(dim), tmptime

      ! initialization
      parts%nhalo_loc = 0

      do i = 1, n_process_neighbour
         neighbours(i)%nhalo_send = 0
         allocate( neighbours(i)%halo_pindex(ntotal_loc))
         neighbours(i)%HaloPackSend%maxn = ntotal_loc
         call neighbours(i)%HaloPackSend%allocate_particles()
      end do

      ! halo particle send location determination
      ! first loop loops over currently held particles
      searchrange(1) = 1
      searchrange(2) = old_ntotal_loc
      tmptime = -system_clock_timer()
      ! begin search
      xmin_loc = my_bounds(1:dim) + 2._f*scale_k*hsml
      xmax_loc = my_bounds(dim + 1:2*dim) - 2._f*scale_k*hsml
      do k = 1, 1
         !do i = searchrange(1), searchrange(2)
         do i = 1, ntotal_loc
            xi(:) = parts%x(:, i)
            if (any([xi(:) <= xmin_loc(:), xi(:) >= xmax_loc(:)])) then
               do j = 1, n_process_neighbour
                  if (all([xi(:) >= neighbours(j)%bounds(1:dim) - 2._f*scale_k*hsml, &
                           xi(:) <= neighbours(j)%bounds(dim + 1:2*dim) + 2._f*scale_k*hsml])) then
                     neighbours(j)%nhalo_send = neighbours(j)%nhalo_send + 1
                     neighbours(j)%halo_pindex(neighbours(j)%nhalo_send) = i
                     ! neighbours(j)%HaloPackSend(neighbours(j)%nhalo_send) = parts(i)
                     ! neighbours(j)%HaloPackSend(neighbours(j)%nhalo_send)%itype = &
                        ! neighbours(j)%HaloPackSend(neighbours(j)%nhalo_send)%itype + sign(halotype, parts(i)%itype)
                     neighbours(j)%HaloPackSend%itype(neighbours(j)%nhalo_send) = parts%itype(i) + &
                        sign(halotype, parts%itype(i))
                     neighbours(j)%HaloPackSend%indglob(neighbours(j)%nhalo_send) = parts%indglob(i)
                     neighbours(j)%HaloPackSend%rho(neighbours(j)%nhalo_send) = parts%rho(i)
                     neighbours(j)%HaloPackSend%x(:, neighbours(j)%nhalo_send) = parts%x(:, i)
                     neighbours(j)%HaloPackSend%vx(:, neighbours(j)%nhalo_send) = parts%vx(:, i)
                     
                  end if
               end do
            end if
         end do

      end do

      do n = 1, n_process_neighbour
         call MPI_Irecv(neighbours(n)%nhalo_recv, 1, MPI_INTEGER, neighbours(n)%rank, 0, MPI_COMM_WORLD, &
            request(2*n-1), ierr)
         call MPI_Isend(neighbours(n)%nhalo_send, 1, MPI_INTEGER, neighbours(n)%rank, 0, MPI_COMM_WORLD, &
            request(2*n), ierr)
      end do

      call MPI_Waitall(2*n_process_neighbour, request, status, ierr)

      displ = ntotal_loc
      do n = 1, n_process_neighbour
         if (neighbours(n)%nhalo_recv > 0) then
            call MPI_Irecv(parts%itype(displ+1), neighbours(n)%nhalo_recv, &
               MPI_INTEGER, neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*5*(n-1)+1), ierr)
            call MPI_Irecv(parts%indglob(displ+1), neighbours(n)%nhalo_recv, &
               MPI_INTEGER, neighbours(n)%rank, 1, MPI_COMM_WORLD, request(2*5*(n-1)+2), ierr)
            call MPI_Irecv(parts%rho(displ+1), neighbours(n)%nhalo_recv, &
               MPI_DOUBLE_PRECISION, neighbours(n)%rank, 2, MPI_COMM_WORLD, request(2*5*(n-1)+3), ierr)
            call MPI_Irecv(parts%x(1, displ+1), dim*neighbours(n)%nhalo_recv, &
               MPI_DOUBLE_PRECISION, neighbours(n)%rank, 3, MPI_COMM_WORLD, request(2*5*(n-1)+4), ierr)
            call MPI_Irecv(parts%vx(1, displ+1), dim*neighbours(n)%nhalo_recv, &
               MPI_DOUBLE_PRECISION, neighbours(n)%rank, 4, MPI_COMM_WORLD, request(2*5*(n-1)+5), ierr)

            displ = displ + neighbours(n)%nhalo_recv
         else
            request(2*5*(n-1)+1:2*5*(n-1)+5) = MPI_REQUEST_NULL
         end if

         if (neighbours(n)%nhalo_send > 0) then
            call MPI_Isend(neighbours(n)%HaloPackSend%itype(1), neighbours(n)%nhalo_send, MPI_INTEGER, &
               neighbours(n)%rank, 0, MPI_COMM_WORLD, request(2*5*(n-1)+6), ierr)
            call MPI_Isend(neighbours(n)%HaloPackSend%indglob(1), neighbours(n)%nhalo_send, MPI_INTEGER, &
               neighbours(n)%rank, 1, MPI_COMM_WORLD, request(2*5*(n-1)+7), ierr)
            call MPI_Isend(neighbours(n)%HaloPackSend%rho(1), neighbours(n)%nhalo_send, MPI_DOUBLE_PRECISION, &
               neighbours(n)%rank, 2, MPI_COMM_WORLD, request(2*5*(n-1)+8), ierr)
            call MPI_Isend(neighbours(n)%HaloPackSend%x(1, 1), dim*neighbours(n)%nhalo_send, MPI_DOUBLE_PRECISION, &
               neighbours(n)%rank, 3, MPI_COMM_WORLD, request(2*5*(n-1)+9), ierr)
            call MPI_Isend(neighbours(n)%HaloPackSend%vx(1, 1), dim*neighbours(n)%nhalo_send, MPI_DOUBLE_PRECISION, &
               neighbours(n)%rank, 4, MPI_COMM_WORLD, request(2*5*(n-1)+10), ierr)
         else
            request(2*5*(n-1)+6:2*5*(n-1)+10) = MPI_REQUEST_NULL
         end if
      end do

      parts%nhalo_loc = displ - ntotal_loc

      call MPI_Waitall(2*5*n_process_neighbour, request, status, ierr)

      do n = 1, n_process_neighbour
         deallocate(neighbours(n)%halo_pindex)
         call neighbours(n)%HaloPackSend%deallocate_particles
      end do

      tmptime = tmptime + system_clock_timer()

#endif

   end subroutine ORB_sendrecv_halo

end module ORB_sr_m
