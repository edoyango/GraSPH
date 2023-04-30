module ORB_m

   use datatypes, only: particles, time_tracking, system_clock_timer
#ifdef PARALLEL
   use mpi_f08
#endif
   ! use ORB_sr_m, only: ORB_sendrecv_diffuse, ORB_sendrecv_halo
   use param, only: f, dim, hsml, scale_k
   use param_para, only: neighbour_data, partition_tracking, dcell_ORB, box_ratio_threshold

   private
   real(f):: box_ratio_previous(dim, dim) = TINY(1._f)
   type(neighbour_data), allocatable:: neighbours(:)
   type(partition_tracking):: partition_track
   integer:: prev_load, n_process_neighbour, repartition_mode
   integer, allocatable:: node_cax(:), pincell_ORB(:, :, :), cellmins(:, :), cellmaxs(:, :), &
      neighbour_ranks(:)
   real(f), allocatable:: my_bounds(:)
   real(f), parameter:: dcell = hsml*dcell_ORB

   public:: ORB, partition_track, neighbours, n_process_neighbour

contains

   subroutine ORB(itimestep, my_rank, num_ranks, parts, timings)

      use param_para, only: ORBcheck1, ORBcheck2

      implicit none
      integer, intent(in):: my_rank, num_ranks, itimestep
      type(particles), intent(inout):: parts
      type(time_tracking), intent(inout):: timings
      integer:: tree_layers, maxnode, stepsSincePrevious, ngridx(3), i, j, k, d, imageRange_ini(2), gridind_ini(3, 2), &
                old_ntv_loc, ngridx_real(3)
      integer:: ntv_loc
      real(f):: mingridx_ini(3), maxgridx_ini(3), current_to_previous(dim, dim), box_ratio_current(dim, dim)
      logical:: free_face(2*dim)

! preprocessor directive ensures subroutine exits straightaway if not in parallel mode
#ifdef PARALLEL

      !allocating partitioning arrays and initialising diagnostic variables ------------------------------------------
      timings%t_ORB = timings%t_ORB - system_clock_timer() ! commence timing of ORB algoirthm
      if (itimestep == 0) then
         tree_layers = ceiling(log(real(num_ranks, kind=f))/log(2.d0))
         maxnode = 2*2**tree_layers - 1
         allocate (node_cax(maxnode), neighbours(num_ranks))
         allocate (my_bounds(2*dim), cellmins(3, num_ranks), cellmaxs(3, num_ranks), neighbour_ranks(num_ranks))
      end if

      ! Boundary Determiniation Algorithm ----------------------------------------------------------------------------
      repartition_mode = 1 !initially assumes no partition
      ! only checks if boundary needs updating every 50-100 time-steps
      stepsSincePrevious = itimestep - partition_track%prev_part_tstep
      if ((itimestep == 0) .or. &
          ((stepsSincePrevious >= ORBcheck1) .and. &
           (mod(stepsSincePrevious, ORBcheck2) == 0))) then

         ! checking if change in particles on current image is > 5%
         if (parts%ntotal_loc > prev_load + &
            0.05_f*real(parts%ntotal, kind=f)/real(my_rank, kind=f)) then
            repartition_mode = 2
         else if (itimestep == 0) then 
            repartition_mode = 3
         end if

         call MPI_Allreduce(MPI_IN_PLACE, repartition_mode, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD)

         if (repartition_mode > 1) then

            partition_track%n_parts = partition_track%n_parts + 1
            if (itimestep /= 1) then
               partition_track%mintstep_bn_part = min(partition_track%mintstep_bn_part, stepsSincePrevious)
               partition_track%maxtstep_bn_part = max(partition_track%maxtstep_bn_part, stepsSincePrevious)
            end if
            partition_track%prev_part_tstep = itimestep

            ! Creating particle-in-cell grid
            call particle_grid(my_rank, num_ranks, parts, ngridx, mingridx_ini, maxgridx_ini, ngridx_real)

            do d = 1, dim
               box_ratio_current(:, d) = real(ngridx_real(d), kind=f)/real(ngridx_real(:), kind=f)
            end do

            current_to_previous(:, :) = box_ratio_current(:, :)/box_ratio_previous(:, :)
            if (any(current_to_previous(:, :) > 1._f + box_ratio_threshold)) repartition_mode = 3

            ! partition summary info
            if (repartition_mode == 3) then
               box_ratio_previous(:, :) = box_ratio_current(:, :)
               if (itimestep /= 1) then
                  partition_track%maxtstep_bn_part = &
                     max(partition_track%maxtstep_bn_reorient, itimestep - partition_track%prev_reorient_tstep)
                  partition_track%mintstep_bn_reorient = &
                     min(partition_track%mintstep_bn_reorient, itimestep - partition_track%prev_reorient_tstep)
               end if
               partition_track%prev_reorient_tstep = itimestep
               partition_track%n_reorients = partition_track%n_reorients + 1
            end if

            if (my_rank == 0) then
               write(*, '(A)') "_______________________________________________________________________________"
               write(*, '(A)') "INFO: Image subdomain boundaries being updated"
               write(*, '(6x, A, I0)') "Current timestep: ", itimestep
               if (repartition_mode==2) write(*, '(6x, A)') "Repartition mode: Updating cut locations only"
               if (repartition_mode==3) write(*, '(6x, A)') "Repartition mode: Updating cut orientations and locations"
            end if

            ! determine subdomain boundaries using particle distribution
            gridind_ini(:, 1) = 1
            gridind_ini(:, 2) = ngridx(:)
            imageRange_ini(1) = 0
            imageRange_ini(2) = num_ranks-1
            free_face(2*dim) = .true.
            n_process_neighbour = 0

            call ORB_bounds(my_rank, num_ranks, gridind_ini, num_ranks, 1, imagerange_ini, parts%ntotal, &
                            dcell, mingridx_ini, maxgridx_ini, free_face)

            deallocate (pincell_ORB)

         end if

      end if

      timings%t_ORB = timings%t_ORB + system_clock_timer() ! conclude timing of ORB algorithm

      ! Particle distribution (physical, halo) -------------------------------------------------------------------------
      timings%t_dist = timings%t_dist - system_clock_timer() ! commence timing of particle distribution

      ! physical particle distribution
      ! ntv_loc = ntotal_loc + nvirt_loc
      ! call ORB_sendrecv_diffuse(itimestep, my_rank, my_bounds, repartition_mode, n_process_neighbour, neighbours, &
      !                           old_ntv_loc, ntv_loc, parts)

      ! ! halo particle interaction
      ! nhalo_loc = 0
      ! call ORB_sendrecv_halo(my_rank, my_bounds, n_process_neighbour, neighbours, old_ntv_loc, &
      !                        ntv_loc, nhalo_loc, parts)

      ! ntotal_loc = 0
      ! nvirt_loc = 0
      ! do i = 1, ntv_loc
      !    if (parts(i)%itype>0) ntotal_loc = ntotal_loc + 1
      !    if (parts(i)%itype<0) nvirt_loc = nvirt_loc + 1
      ! end do

      ! update virtual particles
      ! call generate_virt_part(my_rank, num_ranks, ntotal, ntotal_loc, nvirt, nvirt_loc, parts)

      timings%t_dist = timings%t_dist + system_clock_timer()

#endif

   end subroutine ORB

   !====================================================================================================================
   subroutine particle_grid(my_rank, num_ranks, parts, ngridx, mingridx, maxgridx, ngridx_real)

      implicit none
      integer, intent(in):: my_rank, num_ranks
      type(particles), intent(in):: parts
      integer, intent(out):: ngridx(3), ngridx_real(3)
      real(f), intent(out):: mingridx(3), maxgridx(3)
      real(f):: minx(3), maxx(3), minx_real(3), maxx_real(3), mingridx_real(3), maxgridx_real(3)
      integer:: i, j, k, d, cellrange(3), icell(3), n_nonzerocells, n, rem_rank, cellmins_loc(3), cellmaxs_loc(3)
      type(MPI_Request):: request(4)
      type(MPI_Status):: status(4)

      n = parts%ntotal_loc+parts%nvirt_loc

      ! Local max, min, in each direction ------------------------------------------------------------------------------
      minx(1:dim) = huge(1._f)
      maxx(1:dim) = -huge(1._f)
      minx_real(1:dim) = huge(1._f)
      maxx_real(1:dim) = -huge(1._f)
      do i = 1, n
         do d = 1, dim
            maxx(d) = max(maxx(d), parts%x(d, i))
            minx(d) = min(minx(d), parts%x(d, i))
            if (parts%itype(i) > 0) then
               maxx_real(d) = max(maxx_real(d), parts%x(d, i))
               minx_real(d) = min(minx_real(d), parts%x(d, i))
            end if
         end do
      end do

      if (dim==2) then
         minx(3) = 0._f
         maxx(3) = 0._f
         maxx_real(3) = 0._f
         minx_real(3) = 0._f
      end if

      ! Global max, min, in each direction -----------------------------------------------------------------------------
      mingridx(:) = minx(:)
      maxgridx(:) = maxx(:)
      mingridx_real(:) = minx_real(:)
      maxgridx_real(:) = maxx_real(:)
      call MPI_Iallreduce(MPI_IN_PLACE, maxgridx, dim, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, request(1))
      call MPI_Iallreduce(MPI_IN_PLACE, mingridx, dim, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, request(2))
      call MPI_Iallreduce(MPI_IN_PLACE, maxgridx_real, dim, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, request(3))
      call MPI_Iallreduce(MPI_IN_PLACE, mingridx_real, dim, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, request(4))
      call MPI_Waitall(4, request, status)

      ! Number of grid cells and adjusting max extent in each direction ------------------------------------------------
      mingridx(:) = mingridx(:) - 0.5_f*dcell
      maxgridx(:) = maxgridx(:) + 0.5_f*dcell
      ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
      maxgridx(:) = mingridx(:) + dcell*ngridx(:)

      mingridx_real(:) = mingridx_real(:) - 0.5_f*dcell
      maxgridx_real(:) = maxgridx_real(:) + 0.5_f*dcell
      ngridx_real(:) = int((maxgridx_real(:) - mingridx_real(:))/dcell) + 1
      maxgridx_real(:) = mingridx_real(:) + dcell*ngridx_real(:)

      ! Reducing search area by locating indices in which local particles are contained within
      cellmins_loc(:) = int((minx(:) - mingridx(:))/dcell) + 1
      cellmaxs_loc(:) = int((maxx(:) - mingridx(:))/dcell) + 1
      cellrange(:) = cellmaxs_loc(:) - cellmins_loc(:) + 1

      call MPI_Iallgather(cellmins_loc, 3, MPI_INTEGER, cellmins, 3, MPI_INTEGER, MPI_COMM_WORLD, request(1))
      call MPI_Iallgather(cellmaxs_loc, 3, MPI_INTEGER, cellmaxs, 3, MPI_INTEGER, MPI_COMM_WORLD, request(2))

      ! Creating list of non-zero grid cells ---------------------------------------------------------------------------
      allocate (pincell_ORB(cellrange(1), cellrange(2), cellrange(3)))

      pincell_ORB(:, :, :) = 0

      ! Counting local particles in each cell, and tracking non-zero cells
      n_nonzerocells = 0
      do i = 1, n
         if (parts%itype(i) > 0) then
            icell(1:dim) = int((parts%x(:, i) - mingridx(:))/dcell) + 1
            if (dim == 2) icell(3) = 1
            icell(:) = icell(:) - cellmins_loc(:) + 1
            pincell_ORB(icell(1), icell(2), icell(3)) = pincell_ORB(icell(1), icell(2), icell(3)) + 1
         end if
      end do

      ! Ensuring all MPI processes leave with correct view of cellmins/maxs
      call MPI_Waitall(2, request, status)

   end subroutine particle_grid

   !====================================================================================================================
   recursive subroutine ORB_bounds(my_rank, num_ranks, gridind_in, num_ranks_in, node_in, procrange_in, &
                                   ntotal_in, dcell, mingridx_in, maxgridx_in, free_face)

      ! Recursive function that performs the 'bisection' part of the ORB algorithm
      ! At the end when each image reaches its leaf node, all images sync and then inspect other images' boundaries
      ! for neighbourness

      use param_para, only: bound_extend

      implicit none
      integer, intent(in):: my_rank, num_ranks, gridind_in(3, 2), num_ranks_in, node_in, &
                            procrange_in(2), ntotal_in
      real(f), intent(in):: mingridx_in(3), maxgridx_in(3), dcell
      logical, intent(inout):: free_face(6) ! array to track which faces are cuts, and which faces are free
      integer:: i, ii, j, k, d, ngridx_trim(3), A(3), pincol, np_per_node, my_limits(2, 3), rem_rank, &
                n, procrange_out(2), gridind_out(3, 2), ntotal_out, num_ranks_out, node_out, procrange_lo(2), &
                procrange_hi(2)
      integer:: cax, cut_loc, n_p
      integer, allocatable:: gridsums_loc(:)
      real(f):: rem_bounds(6), all_bounds(6, num_ranks_in)

      !determining cut axis. 1 = x, 2 = y ------------------------------------------------------------------------------
      if (repartition_mode == 3) then
         call P_trim(my_rank, num_ranks, gridind_in, ngridx_trim)
         A(1) = ngridx_trim(2)*ngridx_trim(3)
         A(2) = ngridx_trim(1)*ngridx_trim(3)
         A(3) = ngridx_trim(1)*ngridx_trim(2)
         cax = minloc(A, 1)
         node_cax(node_in) = cax
      else
         cax = node_cax(node_in)
      end if

      allocate(gridsums_loc(gridind_in(cax,2)-gridind_in(cax,1)+1))

      do i = gridind_in(cax,1), gridind_in(cax,2)
         ii = i - gridind_in(cax, 1) + 1
         if (i >= cellmins(cax, my_rank+1) .and. i <= cellmaxs(cax, my_rank+1)) then
            do d = 1, 3
               my_limits(1, d) = max(1, gridind_in(d, 1) - cellmins(d, my_rank+1) + 1)
               my_limits(2, d) = min(cellmaxs(d, my_rank+1), gridind_in(d, 2)) - cellmins(d, my_rank+1) + 1
            end do
            my_limits(1:2, cax) = i - cellmins(cax, my_rank+1) + 1
         
            gridsums_loc(ii) = sum(pincell_ORB(my_limits(1, 1):my_limits(2, 1), &
                                               my_limits(1, 2):my_limits(2, 2), &
                                               my_limits(1, 3):my_limits(2, 3)))
         else
            gridsums_loc(ii) = 0
         end if
      end do

      ! call co_sum(gridsums_loc)
      call MPI_Allreduce(MPI_IN_PLACE, gridsums_loc, gridind_in(cax,2)-gridind_in(cax,1)+1, MPI_INTEGER, MPI_SUM, &
         MPI_COMM_WORLD)

      !cut location ----------------------------------------------------------------------------------------------------
      cut_loc = gridind_in(cax, 1) - 1
      n_p = 0
      np_per_node = int(ceiling(real(num_ranks_in)/2)/real(num_ranks_in)*ntotal_in)
      do while (n_p < np_per_node)
         cut_loc = cut_loc + 1
         n_p = n_p + gridsums_loc(cut_loc-gridind_in(cax,1)+1)
      end do

      ! extra step to nudge cut location backwards if that's a better position
      if ((np_per_node - (n_p - gridsums_loc(cut_loc-gridind_in(cax,1)+1)) < n_p - np_per_node)) then
         cut_loc = cut_loc - 1
         n_p = n_p - gridsums_loc(cut_loc-gridind_in(cax,1)+1)
      end if

      deallocate(gridsums_loc)

      !saving output information ---------------------------------------------------------------------------------------
      procrange_lo(1) = procrange_in(1)
      procrange_lo(2) = procrange_in(1) + ceiling(real(num_ranks_in)/2) - 1
      procrange_hi(1) = procrange_lo(2) + 1
      procrange_hi(2) = procrange_lo(1) + num_ranks_in - 1
      gridind_out(:, 1) = gridind_in(:, 1)
      gridind_out(:, 2) = gridind_in(:, 2)

      node_out = 2*node_in
      procrange_out(:) = procrange_lo(:)
      ntotal_out = n_p
      gridind_out(:, 1) = gridind_in(:, 1)
      gridind_out(:, 2) = gridind_in(:, 2)
      gridind_out(cax, 1) = gridind_in(cax, 1)
      gridind_out(cax, 2) = cut_loc
      free_face(3 + cax) = .false. ! "low" processes have cut on "upper" face
      num_ranks_out = procrange_out(2) - procrange_out(1) + 1
      if (num_ranks_out > 1) then
         call ORB_bounds(my_rank, num_ranks, gridind_out, num_ranks_out, node_out, procrange_out, &
                         ntotal_out, dcell, mingridx_in, maxgridx_in, free_face)
      else if (my_rank==procrange_out(1)) then
         my_bounds(1:3) = mingridx_in(:) + (gridind_out(:, 1) - 1)*dcell
         my_bounds(4:6) = mingridx_in(:) + gridind_out(:, 2)*dcell
      end if

      node_out = 2*node_in + 1
      procrange_out(:) = procrange_hi(:)
      ntotal_out = ntotal_in - n_p
      gridind_out(:, 1) = gridind_in(:, 1)
      gridind_out(:, 2) = gridind_in(:, 2)
      gridind_out(cax, 1) = cut_loc + 1
      gridind_out(cax, 2) = gridind_in(cax, 2)
      free_face(cax) = .false. ! "upper" processes have cut on "lower" face
      num_ranks_out = procrange_out(2) - procrange_out(1) + 1
      if (num_ranks_out > 1) then
         call ORB_bounds(my_rank, num_ranks, gridind_out, num_ranks_out, node_out, procrange_out, &
                         ntotal_out, dcell, mingridx_in, maxgridx_in, free_face)
      else if (my_rank==procrange_out(1)) then
         my_bounds(1:3) = mingridx_in(:) + (gridind_out(:, 1) - 1)*dcell
         my_bounds(4:6) = mingridx_in(:) + gridind_out(:, 2)*dcell
      end if

      if (node_in==1) then

         call MPI_Allgather(my_bounds, 6, MPI_DOUBLE_PRECISION, all_bounds, 6, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD)

         do rem_rank = 0, num_ranks - 1
            rem_bounds(:) = all_bounds(:, rem_rank+1)
            if (.not. (any(my_bounds(1:3) - 0.5_f*scale_k*hsml > rem_bounds(4:6) + 0.5_f*scale_k*hsml) .or. &
                       any(my_bounds(4:6) + 0.5_f*scale_k*hsml < rem_bounds(1:3) - 0.5_f*scale_k*hsml))) then
               n_process_neighbour = n_process_neighbour + 1
               neighbours(n_process_neighbour)%rank = rem_rank
               neighbours(n_process_neighbour)%bounds(:) = rem_bounds(:)
               write(*, *) my_rank, rem_rank, rem_bounds
            end if
         end do

      end if

   end subroutine ORB_bounds

   !====================================================================================================================
   subroutine P_trim(my_rank, num_ranks, gridind_in, ngridx_trim)
      ! Trims particle-in-cell grid so as to obtain minimal bounding boxes to obtain accurate cut axis orientations

      implicit none
      integer, intent(in):: my_rank, num_ranks, gridind_in(3, 2)
      integer, intent(out):: ngridx_trim(3)
      integer:: i, j, k, n, gridind_new(3, 2)
      type(MPI_Request):: request(2)
      type(MPI_Status):: status(2)

      !Trimming input grid ------------------------------------------------------------------------------------------------------
      !reducing x-length of grid
      !finding new start x-index
      gridind_new(1, 1) = trim_helper(my_rank, num_ranks, 1, gridind_in(1, 1), gridind_in(1, 2), 2, &
                                      gridind_in(2, 1), gridind_in(2, 2), 3, gridind_in(3, 1), gridind_in(3, 2))

      !finding new end x-index
      gridind_new(1, 2) = trim_helper(my_rank, num_ranks, 1, gridind_in(1, 2), gridind_in(1, 1), 2, &
                                      gridind_in(2, 1), gridind_in(2, 2), 3, gridind_in(3, 1), gridind_in(3, 2))

      !reducing y-length of grid
      !finding new start y-index
      gridind_new(2, 1) = trim_helper(my_rank, num_ranks, 2, gridind_in(2, 1), gridind_in(2, 2), 1, &
                                      gridind_in(1, 1), gridind_in(1, 2), 3, gridind_in(3, 1), gridind_in(3, 2))

      !finding new end y-index
      gridind_new(2, 2) = trim_helper(my_rank, num_ranks, 2, gridind_in(2, 2), gridind_in(2, 1), 1, &
                                      gridind_in(1, 1), gridind_in(1, 2), 3, gridind_in(3, 1), gridind_in(3, 2))

      !reducing z-length of grid
      !finding new start z-index
      gridind_new(3, 1) = trim_helper(my_rank, num_ranks, 3, gridind_in(3, 1), gridind_in(3, 2), 1, &
                                      gridind_in(1, 1), gridind_in(1, 2), 2, gridind_in(2, 1), gridind_in(2, 2))

      !finding new end z-index
      gridind_new(3, 2) = trim_helper(my_rank, num_ranks, 3, gridind_in(3, 2), gridind_in(3, 1), 1, &
                                      gridind_in(1, 1), gridind_in(1, 2), 2, gridind_in(2, 1), gridind_in(2, 2))

      ! call co_min(gridind_new(:, 1))
      ! call co_max(gridind_new(:, 2))
      call MPI_Iallreduce(MPI_IN_PLACE, gridind_new(:, 1), 3, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, request(1))
      call MPI_Iallreduce(MPI_IN_PLACE, gridind_new(:, 2), 3, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, request(2))

      call MPI_Waitall(2, request, status)

      ngridx_trim(:) = gridind_new(:, 2) - gridind_new(:, 1) + 1

   end subroutine P_trim

   !=================================================================================
   function trim_helper(my_rank, num_ranks, main_axis, main_start, main_end, off_axis1, off_limits1_start, &
                        off_limits1_end, off_axis2, off_limits2_start, off_limits2_end) result(newindex)
      ! function that performs the scanning to see if there are any non-zero cells

      integer, intent(in):: my_rank, num_ranks, main_axis, main_start, main_end, off_axis1, off_axis2, &
         off_limits1_start, off_limits1_end, off_limits2_start, off_limits2_end
      integer:: newindex, step, a_main, n, my_limits(2, 3), my_rank_maxincell, cellmins_loc(3), cellmaxs_loc(3)

      newindex = main_start

      step = sign(1, main_end - main_start)

      do a_main = main_start, main_end, step
         if (.not. (a_main > cellmaxs(main_axis, my_rank+1) .or. &
                    off_limits1_start > cellmaxs(off_axis1, my_rank+1) .or. &
                    off_limits2_start > cellmaxs(off_axis2, my_rank+1) .or. &
                    cellmins(main_axis, my_rank+1) > a_main .or. &
                    cellmins(off_axis1, my_rank+1) > off_limits1_end .or. &
                    cellmins(off_axis2, my_rank+1) > off_limits2_end)) then
            my_limits(1:2, main_axis) = a_main - cellmins(main_axis, my_rank+1) + 1
            my_limits(1, off_axis1) = max(cellmins(off_axis1, my_rank+1), off_limits1_start) - &
                                              cellmins(off_axis1, my_rank+1) + 1
            my_limits(2, off_axis1) = min(cellmaxs(off_axis1, my_rank+1), off_limits1_end) - &
                                              cellmins(off_axis1, my_rank+1) + 1
            my_limits(1, off_axis2) = max(cellmins(off_axis2, my_rank+1), off_limits2_start) - &
                                              cellmins(off_axis2, my_rank+1) + 1
            my_limits(2, off_axis2) = min(cellmaxs(off_axis2, my_rank+1), off_limits2_end) - &
                                              cellmins(off_axis2, my_rank+1) + 1
            my_rank_maxincell = maxval(pincell_ORB(my_limits(1, 1):my_limits(2, 1), &
                                                      my_limits(1, 2):my_limits(2, 2), &
                                                      my_limits(1, 3):my_limits(2, 3)))
            if (my_rank_maxincell > 0) Then
               newindex = a_main
               return
            end if
         end if
      end do

      newindex = main_end

   end function trim_helper

end module ORB_m

