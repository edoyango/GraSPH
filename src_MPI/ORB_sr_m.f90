module ORB_sr_m

   use globvar, only: ntotal_loc, nhalo_loc, parts, maxnloc, itimestep, scale_k
   use globvar_para, only: n_process_neighbour, parttype, neighbour_data
   use mpi_f08
   use param, only: f, dim, hsml
   !use error_msg_m, only: error_msg
   !use output_m, only: output

   private
   public:: ORB_sendrecv_diffuse, ORB_sendrecv_halo, ORB_sendrecv_haloupdate

contains

   !==============================================================================================================================
   subroutine ORB_sendrecv_diffuse(procid, bounds_loc, repartition_mode, neighbours, nrequest, request, n_recv_all)
      ! Recursive function to exchange physical particles. In cases were subdomain boundaries are updated, the possibility of needing
      ! diffusion is considered

      implicit none
      integer, intent(in):: procid, repartition_mode
      real(f), intent(in):: bounds_loc(2*dim)
      type(neighbour_data), intent(inout):: neighbours(:)
      integer, intent(out):: nrequest, n_recv_all
      type(MPI_Request), intent(out):: request(:)
      integer:: d, i, j, n, pos_recv, ierr
      integer:: nphys_send_all, diff_dest, ndiffuse_loc, ndiffuse_all, searchrange(2), entrydepth
      real(f):: xmin_loc(dim), xmax_loc(dim), xmin_rem(dim), xmax_rem(dim), xi(dim), dr, dr_min
      type(MPI_STATUS):: status(2*n_process_neighbour + 1)
      logical:: diffuse
      integer, allocatable:: removal_list(:)

      ! Initialization
      searchrange(:) = [1, ntotal_loc]
      nrequest = 0
      diffuse = .true.
      entrydepth = 0

      xmin_loc(:) = bounds_loc(1:dim)
      xmax_loc(:) = bounds_loc(dim+1:2*dim)

      do while (diffuse)

         ! Searching particles to remove within indices of searchrange(1) and searchrange(2), inclusive.
         ! At node 0, searchrange(1:2) = [1,ntotal_loc)
         allocate (removal_list(ntotal_loc + 1))
         nphys_send_all = 0
         do i = searchrange(1), searchrange(2)
            xi(:) = parts(i)%x(:)
            if (any([xi(:) .lt. xmin_loc(:), xi(:) .ge. xmax_loc(:)])) then
               nphys_send_all = nphys_send_all + 1
               removal_list(nphys_send_all) = i
            end if
         end do
         removal_list(nphys_send_all + 1) = 0 ! This is needed due to a quirk in the loop below

         ! If there are any particles that do not belong to the host process,
         ! begin searching for neighbouring processes to send the particle to.
         ! If particle is not contained with subdomain boundaries, send the particle to nearest process neighbouring current host.
         !allocate (PhysPackSend(nphys_send_all, n_process_neighbour))
         do i = 1, n_process_neighbour
            neighbours(i)%nphys_send = 0
            allocate (neighbours(i)%PhysPackSend(nphys_send_all))
         end do
         ndiffuse_loc = 0
         if (nphys_send_all .gt. 0) then
            loop_through_parts: do j = 1, nphys_send_all
               i = removal_list(j)
               xi(:) = parts(i)%x(:) !placeholder variable of particle position for code readibility
               do n = 1, n_process_neighbour
                  xmin_rem(:) = neighbours(n)%bounds(1:dim)
                  xmax_rem(:) = neighbours(n)%bounds(dim + 1:2*dim)
                  if (.not. any([xi(:) .lt. xmin_rem(:), xi(:) .ge. xmax_rem(:)])) then

                     call neighbours(n)%Pack_PhysPart(parts(i))

                     cycle loop_through_parts ! stop searching for a subdomain neighbour and move to next particle
                  end if
               end do

               ! particle belongs to non-neighbouring process.
               ! Evaluates closeness by considering minimum distance between particle and a processes edge,face, or vertice
               ndiffuse_loc = ndiffuse_loc + 1
               dr_min = huge(1._f)
               diff_dest = 0
               do n = 1, n_process_neighbour
                  dr = 0_f
                  do d = 1, dim
                     dr = dr + MAX(0d0, neighbours(n)%bounds(d) - parts(i)%x(d), parts(i)%x(d) - neighbours(n)%bounds(dim + d))**2
                  end do
                  dr = SQRT(dr)
                  if (dr .lt. dr_min) then
                     diff_dest = n
                     dr_min = dr
                  end if
               end do
               call neighbours(diff_dest)%Pack_PhysPart(parts(i))

               !call error_msg(itimestep,procid,4,ind(i))
            end do loop_through_parts
         end if

         ! Posting non-blocking send/recv to exchange info with neighers of no. particles being sent
         do n = 1, n_process_neighbour
            call MPI_IRECV(neighbours(n)%nphys_recv, 1, MPI_INTEGER, neighbours(n)%pid, 0, MPI_COMM_WORLD, request(2*n - 1), ierr)
            call MPI_ISEND(neighbours(n)%nphys_send, 1, MPI_INTEGER, neighbours(n)%pid, 0, MPI_COMM_WORLD, request(2*n), ierr)
         end do

         ! Shifting information to remove sent particles
         n = 0
         if ((nphys_send_all .gt. 0) .and. (nphys_send_all .lt. ntotal_loc)) then
            do i = removal_list(1), ntotal_loc
               if (i .eq. removal_list(n + 1)) then
                  n = n + 1
               else
                  parts(i - n) = parts(i) ! moving particle to new local position
               end if
            end do
         end if

         ntotal_loc = ntotal_loc - nphys_send_all

         ! Wait for non-blocking send to complete
         call MPI_WAITALL(2*n_process_neighbour, request, status, ierr)

         ! Calculating total number of physical particles to be received. Check if this exceeds particle array boundaries
         n_recv_all = SUM(neighbours(1:n_process_neighbour)%nphys_recv)

         !if (ntotal_loc + n_recv_all .gt. maxnloc) call error_msg(2, parts(ntotal_loc)%indloc)

         ! Non-blocking sends to exchange physical particles that have moved processes ---------------------------------------------------
         pos_recv = ntotal_loc + 1
         nrequest = 0
         do n = 1, n_process_neighbour

            if (neighbours(n)%nphys_send .gt. 0) then

               nrequest = nrequest + 1
               call MPI_ISEND(neighbours(n)%PhysPackSend(1), neighbours(n)%nphys_send, parttype, neighbours(n)%pid, 0, &
                              MPI_COMM_WORLD, request(nrequest), ierr)
            end if

            if (neighbours(n)%nphys_recv .gt. 0) then

               nrequest = nrequest + 1
               call MPI_IRECV(parts(pos_recv), neighbours(n)%nphys_recv, parttype, neighbours(n)%pid, 0, MPI_COMM_WORLD, &
                              request(nrequest), ierr)
               pos_recv = pos_recv + neighbours(n)%nphys_recv

            end if

         end do

         ! if subdomain boundary update has occurred, check if diffusion is necessary.
         ! Perform if necessary
         if (repartition_mode .ge. 2) then

            ! Checking if any process has particles that need to be diffused
            call MPI_ALLREDUCE(ndiffuse_loc, ndiffuse_all, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

            ! If any particles need to be diffused, repeat process on newly received particles
            if (ndiffuse_all .eq. 0) then
               diffuse = .false.
            else
               if (procid .eq. 0) then
                  write (*, '(A48,I7)') 'Diffusion occuring... Current timestep:         ', itimestep
                  write (*, '(A48,I7)') '                      Current depth:            ', entrydepth
                  write (*, '(A48,I7)') '                      Particles to be diffused: ', ndiffuse_all
               end if

               ! wait for particle exchange to complete
               call MPI_WAITALL(nrequest, request, status, ierr)

               entrydepth = entrydepth + 1
               searchrange = [ntotal_loc + 1, ntotal_loc + n_recv_all]
               ! Update ntotal_loc as required.
               ntotal_loc = ntotal_loc + n_recv_all

               ! cleaning up to prepare for next iteration
               deallocate (removal_list)
               do i = 1, n_process_neighbour
                  deallocate (neighbours(i)%PhysPackSend)
               end do
            end if
         else
            diffuse = .false.
         end if

      end do

   end subroutine ORB_sendrecv_diffuse

   !==============================================================================================================================
   subroutine ORB_sendrecv_halo(procid, bounds_loc, neighbours, request_in, request_out, nphys_recv_all, nrequest)

      !subroutine responsible for sending sending halo particle information between processes, given
      !predetermined subdomain boundaires.
      !Note: subdomain boundaries are used as inputs (bounds_glob).

      use globvar_para, only: halotype

      implicit none
      integer, intent(in):: procid
      real(f), intent(in):: bounds_loc(2*dim)
      type(neighbour_data), intent(inout):: neighbours(:)
      type(MPI_Request), intent(inout):: request_in(:)
      integer, intent(inout):: nphys_recv_all, nrequest
      type(MPI_Request), intent(out):: request_out(:)
      type(MPI_Status):: status(2*n_process_neighbour)
      integer:: i, j, k, n, pos0_recv, pos1_recv, pos0, pos1, maxloop, ierr
      real(f):: xmin_rem(dim), xmax_rem(dim), xi(dim), xmin_loc(dim), xmax_loc(dim)
      logical:: wait_for_phys

      ! Initialization
      do i = 1, n_process_neighbour
         neighbours(i)%nhalo_send = 0
         if (allocated(neighbours(i)%halo_pindex)) deallocate (neighbours(i)%halo_pindex)
         allocate (neighbours(i)%halo_pindex(ntotal_loc + nphys_recv_all))
      end do

      if (nphys_recv_all .eq. 0) then
         wait_for_phys = .false.
         maxloop = 1
      else
         wait_for_phys = .true.
         maxloop = 2
      end if

      ! Halo particle send location determination
      ! First loop loops over currently held particles
      pos0 = 1
      pos1 = ntotal_loc

      ! Begin search
      xmin_loc = bounds_loc(1:dim) + scale_k*hsml
      xmax_loc = bounds_loc(dim+1:2*dim) - scale_k*hsml
      do k = 1, maxloop
         do i = pos0, pos1
            xi(:) = parts(i)%x(:)
            if (any([xi(:) .le. xmin_loc(:), xi(:) .ge. xmax_loc(:)])) then ! if particle is potentially neighbour's halo
               do j = 1, n_process_neighbour
                  xmin_rem(:) = neighbours(j)%bounds(1:dim) - scale_k*hsml
                  xmax_rem(:) = neighbours(j)%bounds(dim + 1:2*dim) + scale_k*hsml
                  if (all([xi(:) .ge. xmin_rem(:), xi(:) .le. xmax_rem(:)])) then

                     neighbours(j)%nhalo_send = neighbours(j)%nhalo_send + 1
                     neighbours(j)%halo_pindex(neighbours(j)%nhalo_send) = i

                  end if

               end do
            end if
         end do

         ! Waiting for physical particles to complete exchange if needed
         if (wait_for_phys) then
            pos0 = ntotal_loc + 1
            pos1 = ntotal_loc + nphys_recv_all
            wait_for_phys = .false.
            call MPI_WAITALL(nrequest, request_in, status, ierr) !wait for new physical particles to arrive
         end if
      end do

      ! Posting non-blocking send for nhalo_send exchange
      do n = 1, n_process_neighbour
         call MPI_IRECV(neighbours(n)%nhalo_recv, 1, MPI_INTEGER, neighbours(n)%pid, 0, MPI_COMM_WORLD, request_out(2*n - 1), ierr)
         call MPI_ISEND(neighbours(n)%nhalo_send, 1, MPI_INTEGER, neighbours(n)%pid, 0, MPI_COMM_WORLD, request_out(2*n), ierr)
      end do

      ! Creating indexed derived types to send halo particles to each neighbouring process
      do n = 1, n_process_neighbour
         if (neighbours(n)%nhalo_send > 0) then
            call neighbours(n)%create_indexed_halotypes
         end if
      end do

      ntotal_loc = ntotal_loc + nphys_recv_all
      parts(1:ntotal_loc)%indloc = (/(i, i=1, ntotal_loc)/) ! updating all the loca indices of new physical particles

      do i = 1, n_process_neighbour
         deallocate (neighbours(i)%PhysPackSend)
      end do

      ! Wait for non-blocking send to complete
      call MPI_WAITALL(2*n_process_neighbour, request_out, status, ierr)

      ! stopping program if array bounds are exceeded
      nhalo_loc = SUM(neighbours(1:n_process_neighbour)%nhalo_recv)
      !if (ntotal_loc + nhalo_loc .gt. maxnloc) call error_msg(3, parts(i)%indloc)

      ! Non-blocking sends to exchange physical particles that have moved processes
      nrequest = 0
      pos1_recv = ntotal_loc
      do n = 1, n_process_neighbour

         if (neighbours(n)%nhalo_recv .gt. 0) then

            pos0_recv = pos1_recv + 1
            pos1_recv = pos0_recv - 1 + neighbours(n)%nhalo_recv

            nrequest = nrequest + 1

            call MPI_IRECV(parts(pos0_recv), neighbours(n)%nhalo_recv, halotype, neighbours(n)%pid, 0, MPI_COMM_WORLD, &
                           request_out(nrequest), ierr)

         end if

         if (neighbours(n)%nhalo_send .gt. 0) then

            nrequest = nrequest + 1

            call MPI_ISEND(parts(1), 1, neighbours(n)%halotype_indexed, neighbours(n)%pid, 0, MPI_COMM_WORLD, &
                           request_out(nrequest), ierr)

         end if

      end do

   end subroutine ORB_sendrecv_halo

   !==============================================================================================================================
   subroutine ORB_sendrecv_haloupdate(ki, neighbours)
      ! Reduced version of ORB_sendrecv_halo where no searching occurs (only the exchange)

      use globvar_para, only: haloupdatetype

      implicit none
      integer, intent(in):: ki
      type(neighbour_data), intent(inout):: neighbours(:)
      integer:: n, i, pos0_recv, pos1_recv, n_request, ierr
      type(MPI_Request):: request(2*n_process_neighbour)
      type(MPI_Status):: status(2*n_process_neighbour)

      !3. halo particle send/receive ----------------------------------------------------------------------------------------------------
      n_request = 0
      pos1_recv = ntotal_loc
      do n = 1, n_process_neighbour

         if (neighbours(n)%nhalo_recv .gt. 0) then

            ! Indices where receiving data will be placed
            pos0_recv = pos1_recv + 1
            pos1_recv = pos0_recv - 1 + neighbours(n)%nhalo_recv

            n_request = n_request + 1
            call MPI_IRECV(parts(pos0_recv), neighbours(n)%nhalo_recv, haloupdatetype, neighbours(n)%pid, 0, MPI_COMM_WORLD, &
                           request(n_request), ierr)

         end if

         if (neighbours(n)%nhalo_send .gt. 0) then

            n_request = n_request + 1
            call MPI_ISEND(parts(1), 1, neighbours(n)%haloupdatetype_indexed, neighbours(n)%pid, 0, MPI_COMM_WORLD, &
                           request(n_request), ierr)

         end if

      end do

      ! Freeing up halotype used just for first exchange
      if (ki .eq. 2) then
         do i = 1, n_process_neighbour
            if (neighbours(i)%nhalo_send > 0) call MPI_TYPE_FREE(neighbours(i)%halotype_indexed, ierr)
         end do
      end if

      call MPI_WAITALL(n_request, request, status, ierr)

      ! Freeing up halo type used for 2nd - 3rd exchanges
      if (ki .eq. 4) then
         do i = 1, n_process_neighbour
            if (neighbours(i)%nhalo_send > 0) call MPI_TYPE_FREE(neighbours(i)%haloupdatetype_indexed, ierr)
         end do
      end if

   end subroutine ORB_sendrecv_haloupdate

end module ORB_sr_m
