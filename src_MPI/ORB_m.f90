module ORB_m

   use datatypes, only: particles, time_tracking
   use mpi_f08
   use param, only: f, dim, hsml
   use param_para, only: MPI_derived_types, partition_tracking, neighbour_data

   private
   !Partition frequency variables
   real(f):: box_ratio_previous(dim, dim) = TINY(1.), bounds_loc(2*dim)
   integer:: prev_load, node_cut(2*dim), n_process_neighbour
   type(partition_tracking):: partition_track
   type(neighbour_data), allocatable:: neighbours(:)
   integer, allocatable:: node_cax(:)

   public:: ORB, partition_track, neighbours, n_process_neighbour

contains
   !==============================================================================================================================
   subroutine ORB(itimestep, procid, numprocs, MPI_types, scale_k, ntotal, ntotal_loc, nhalo_loc, nvirt_loc, parts, timings)
      ! Container subroutine for the bulk of the ORB algorithm, including the initial exchange of physical and halo particles

      use param_para, only: dcell_ORB, ORBcheck1, ORBcheck2, box_ratio_threshold
      use ORB_sr_m, only: ORB_sendrecv_diffuse, ORB_sendrecv_halo
      use input_m, only: generate_virt_part

      implicit none
      integer, intent(in):: procid, numprocs, itimestep, ntotal
      real(f), intent(in):: scale_k
      type(MPI_derived_types), intent(in):: MPI_types
      integer, intent(inout):: ntotal_loc, nhalo_loc, nvirt_loc
      type(particles), intent(inout):: parts(:)
      type(time_tracking), intent(inout):: timings
      real(f), parameter:: dcell = hsml*dcell_ORB
      integer:: d, i, ngridx(dim), nphys_recv_all, n_request, procrange_ini(2), tree_layers, maxnode, &
                gridind_ini(dim, 2), repartition_mode_loc, ierr, repartition_mode
      real(f):: bounds_glob(2*dim, numprocs), mingridx_ini(dim), maxgridx_ini(dim), current_to_previous(dim, dim), &
                box_ratio_current(dim, dim)
      type(MPI_Status):: status(4*numprocs)
      type(MPI_Request):: request_phys(2*numprocs), request_halo(2*numprocs)
      integer, allocatable:: pincell_ORB(:, :, :)

      !allocating partitioning arrays and initialising diagnostic variables -------------------------------------------------------------
      timings%t_ORB = timings%t_ORB - MPI_WTIME() ! commence timing of ORB algoirthm
      if (itimestep .eq. 0) then
         tree_layers = CEILING(LOG(DBLE(numprocs))/LOG(2d0))
         maxnode = 2*2**tree_layers - 1
         allocate (node_cax(maxnode), neighbours(numprocs))
      end if

      ! Boundary Determiniation Algorithm ---------------------------------------------------------------------------------------
      repartition_mode = 1 !initially assumes no partition
      ! only checks if boundary needs updating every 50-100 time-steps
      if ((itimestep .eq. 0) .or. ((itimestep - partition_track%prev_part_tstep .ge. ORBcheck1) .and. &
                                   (mod(itimestep - partition_track%prev_part_tstep, ORBcheck2) .eq. 0))) then

         ! checking if change in partilces on current process > 5%
         if (ntotal_loc .gt. prev_load + 0.05_f*DBLE(ntotal)/DBLE(numprocs)) then
            repartition_mode_loc = 2
         end if

         call MPI_ALLREDUCE(repartition_mode_loc, repartition_mode, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

         if ((repartition_mode .gt. 1) .or. (itimestep .eq. 1)) then

            partition_track%n_parts = partition_track%n_parts + 1
            if (itimestep .ne. 1) then
               partition_track%mintstep_bn_part = min(partition_track%mintstep_bn_part, itimestep - partition_track%prev_part_tstep)
               partition_track%maxtstep_bn_part = max(partition_track%maxtstep_bn_part, itimestep - partition_track%prev_part_tstep)
            end if
            partition_track%prev_part_tstep = itimestep

            ! Creating particle-in-cell grid
            call particle_grid(numprocs, MPI_types%ftype, ntotal_loc, parts(1:ntotal_loc), ngridx, dcell, mingridx_ini, &
                               maxgridx_ini, pincell_ORB)

            ! Calculating current aspect ratio.
            do d = 1, dim
               box_ratio_current(:, d) = DBLE(ngridx(d))/DBLE(ngridx(:))
            end do

            current_to_previous(:, :) = box_ratio_current(:, :)/box_ratio_previous(:, :)
            if (any(current_to_previous(:, :) > 1_f + box_ratio_threshold)) repartition_mode = 3

            !partition summary info
            if (repartition_mode .eq. 3) then
               box_ratio_previous(:, :) = box_ratio_current(:, :)
               if (itimestep .ne. 1) then
                  partition_track%maxtstep_bn_reorient = &
                     max(partition_track%maxtstep_bn_reorient, itimestep - partition_track%prev_reorient_tstep)
                  partition_track%mintstep_bn_reorient = &
                     min(partition_track%mintstep_bn_reorient, itimestep - partition_track%prev_reorient_tstep)
               end if
               partition_track%prev_reorient_tstep = itimestep
               partition_track%n_reorients = partition_track%n_reorients + 1
            end if

            ! determine subdomain boundaries using particle distribution
            gridind_ini(:, 1) = 1
            gridind_ini(:, 2) = ngridx(:)
            procrange_ini(1) = 0
            procrange_ini(2) = numprocs - 1
            bounds_glob = ORB_bounds(procid, numprocs, MPI_types%ftype, scale_k, repartition_mode, gridind_ini, numprocs, 1, &
                                     procrange_ini, &
                                     ntotal, pincell_ORB, dcell, mingridx_ini, maxgridx_ini)

            call subdomain_neighbour(procid, numprocs, bounds_glob, scale_k, n_process_neighbour)

            bounds_loc(:) = bounds_glob(:, procid + 1)

            deallocate (pincell_ORB)

         end if

      end if
      timings%t_ORB = timings%t_ORB + MPI_WTIME() ! conclude timing of ORB algorithm

      ! Particle distribution (physical, halo) ----------------------------------------------------------------------------------
      timings%t_dist = timings%t_dist - MPI_WTIME() ! commence timing of particle distribution

      ! physical particle distribution
      call ORB_sendrecv_diffuse(itimestep, procid, bounds_loc, MPI_types%parttype, repartition_mode, n_process_neighbour, &
                                neighbours, n_request, request_phys, nphys_recv_all, ntotal_loc, parts)

      ! halo particle distribution
      call ORB_sendrecv_halo(procid, bounds_loc, scale_k, MPI_types%halotype, MPI_types%haloupdatetype, n_process_neighbour, &
                             neighbours, request_phys, request_halo, nphys_recv_all, n_request, ntotal_loc, nhalo_loc, parts)

      ! update virtual particles
      call generate_virt_part(procid, bounds_loc, scale_k, ntotal, ntotal_loc, nhalo_loc, nvirt_loc, parts)

      if (repartition_mode .gt. 1) prev_load = ntotal_loc

      do i = ntotal_loc + 1, ntotal_loc + nhalo_loc
         parts(i)%indloc = i
      end do

      ! wait for halo particle distribution to complete
      call MPI_WAITALL(n_request, request_halo, status, ierr)

      parts(ntotal_loc + 1:ntotal_loc + nhalo_loc)%itype = 2
      timings%t_dist = timings%t_dist + MPI_WTIME()

   end subroutine ORB

   !==============================================================================================================================
   subroutine particle_grid(numprocs, MPI_ftype, ntotal_loc, parts, ngridx, dcell, mingridx, maxgridx, pincell_ORB)
      ! Subroutine to create a uniform rectangular grid with square cells, and counting the number of particles contained within each
      ! cell. Each MPI process holds a global copy of the entire grid, in preperation for ORB

      implicit none
      integer, intent(in):: numprocs, ntotal_loc
      real(f), intent(in):: dcell
      type(MPI_datatype), intent(in):: MPI_ftype
      type(particles), intent(in):: parts(:)
      integer, intent(out):: ngridx(:)
      real(f), intent(out):: mingridx(:), maxgridx(:)
      integer, allocatable, intent(out):: pincell_ORB(:, :, :)
      integer:: i, d, icell, jcell, kcell, n_nonzerocells, n_nonzerocells_perprocess(numprocs), n_nonzerocells_total, pid, &
                displ(numprocs), cellmins(3), cellmaxs(3), cellrange(3), ierr
      real(f):: minx(3), maxx(3)
      integer:: sendcount, recvcount(numprocs), Plist_size
      integer, allocatable:: Plist_loc(:, :), Plist_all(:, :)
      type(MPI_Status):: status(2)
      type(MPI_Request):: request(2)

      ! Local max, min, in each direction ---------------------------------------------------------------------------------------
      minx(1:dim) = parts(1)%x(1:dim)
      maxx(1:dim) = parts(1)%x(1:dim)
      do i = 2, ntotal_loc
         do d = 1, dim
            if (parts(i)%x(d) .gt. maxx(d)) maxx(d) = parts(i)%x(d)
            if (parts(i)%x(d) .lt. minx(d)) minx(d) = parts(i)%x(d)
         end do
      end do

      ! Global max, min, in each direction --------------------------------------------------------------------------------------
      call MPI_IALLREDUCE(maxx(:), maxgridx(:), dim, MPI_ftype, MPI_MAX, MPI_COMM_WORLD, request(1), ierr) !finding max over all processes
      call MPI_IALLREDUCE(minx(:), mingridx(:), dim, MPI_ftype, MPI_MIN, MPI_COMM_WORLD, request(2), ierr) !finding min over all processes

      call MPI_WAITALL(2, request, status, ierr)

      mingridx(:) = mingridx(:) - dcell
      maxgridx(:) = maxgridx(:) + dcell

      ! Number of grid cells and adjusting max extent in each direction ---------------------------------------------------------
      ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
      maxgridx(:) = mingridx(:) + dcell*ngridx(:)

      ! Reducing search area by locating indices in which local particles are contained within
      cellmins(:) = int((minx(:) - mingridx(:))/dcell) + 1
      cellmaxs(:) = int((maxx(:) - mingridx(:))/dcell) + 1
      cellrange(:) = cellmaxs(:) - cellmins(:) + 1

      ! Creating list of non-zero grid cells ------------------------------------------------------------------------------------
      Plist_size = MIN(ntotal_loc, PRODUCT(cellrange(:)))
      allocate (pincell_ORB(ngridx(1), ngridx(2), ngridx(3)), &
                Plist_loc(4, Plist_size))

      pincell_ORB(cellmins(1):cellmaxs(1), cellmins(2):cellmaxs(2), cellmins(3):cellmaxs(3)) = 0

      ! Counting particles in each cell
      n_nonzerocells = 0
      do i = 1, ntotal_loc
         icell = int((parts(i)%x(1) - mingridx(1))/dcell) + 1
         jcell = int((parts(i)%x(2) - mingridx(2))/dcell) + 1
         kcell = int((parts(i)%x(3) - mingridx(3))/dcell) + 1
         if (pincell_ORB(icell, jcell, kcell) .eq. 0) then
            n_nonzerocells = n_nonzerocells + 1
            Plist_loc(1, n_nonzerocells) = icell
            Plist_loc(2, n_nonzerocells) = jcell
            Plist_loc(3, n_nonzerocells) = kcell
         end if
         pincell_ORB(icell, jcell, kcell) = pincell_ORB(icell, jcell, kcell) + 1
      end do

      ! Exchanging info on how many non-zero cells each process has
      call MPI_IALLGATHER(n_nonzerocells, 1, MPI_INTEGER, n_nonzerocells_perprocess, 1, MPI_INTEGER, MPI_COMM_WORLD, &
                          request(1), ierr)

      ! Populating Plist_loc with non-zero cell entries
      do i = 1, n_nonzerocells
         Plist_loc(4, i) = Pincell_ORB(Plist_loc(1, i), Plist_loc(2, i), Plist_loc(3, i))
      end do

      call MPI_WAIT(request(1), status(1), ierr)

      n_nonzerocells_total = 0
      do pid = 1, numprocs
         displ(pid) = n_nonzerocells_total
         n_nonzerocells_total = n_nonzerocells_total + n_nonzerocells_perprocess(pid)
      end do

      allocate (Plist_all(4, n_nonzerocells_total))

      ! Collecting all process' Plist_loc
      displ = 4*displ
      sendcount = 4*n_nonzerocells
      recvcount = 4*n_nonzerocells_perprocess(:)
      call MPI_IALLGATHERV(Plist_loc, sendcount, MPI_INTEGER, Plist_all, recvcount, displ, MPI_INTEGER, MPI_COMM_WORLD, &
                           request(1), ierr)

      pincell_ORB(:, :, :) = 0

      call MPI_WAIT(request(1), status(1), ierr)

      ! populating the number of particles per cell
      do i = 1, n_nonzerocells_total
         pincell_ORB(Plist_all(1, i), Plist_all(2, i), Plist_all(3, i)) = &
            pincell_ORB(Plist_all(1, i), Plist_all(2, i), Plist_all(3, i)) + Plist_all(4, i)
      end do

      deallocate (Plist_loc, Plist_all)

   end subroutine particle_grid

   !==============================================================================================================================
   recursive function ORB_bounds(procid, numprocs, MPI_ftype, scale_k, repartition_mode, gridind_in, nprocs_in, node_in, &
                                 procrange_in, &
                                 ntotal_in, pincell_ORB, dcell, mingridx_in, maxgridx_in) result(bounds_glob)
      ! Recursive function that performs the 'bisection' part of the ORB algorithm

      use param_para, only: bound_extend

      implicit none
      integer, intent(in):: procid, numprocs, gridind_in(dim, 2), node_in, nprocs_in, procrange_in(2), ntotal_in, &
                            repartition_mode, pincell_ORB(:, :, :)
      type(MPI_datatype), intent(in):: MPI_ftype
      real(f), intent(in):: mingridx_in(dim), maxgridx_in(dim), dcell, scale_k
      integer:: i, node_out, gridind_out(dim, 2), nprocs_out, ntotal_out, procrange_out(2), n_p, cax, np_per_node, pincol, &
                ngridx_trim(dim), A(3), procrange_lo(2), procrange_hi(2), ierr
      real(f):: bounds_out(2*dim), bounds_glob(2*dim, numprocs)

      !determining cut axis. 1 = x, 2 = y ---------------------------------------------------------------------------------------
      if (repartition_mode .eq. 3) then
         call P_trim(gridind_in, ngridx_trim, pincell_ORB)
         A(1) = ngridx_trim(2)*ngridx_trim(3); A(2) = ngridx_trim(1)*ngridx_trim(3); A(3) = ngridx_trim(1)*ngridx_trim(2)
         if ((A(1) .le. A(2)) .and. (A(1) .le. A(3))) cax = 1
         if ((A(2) .le. A(1)) .and. (A(2) .le. A(3))) cax = 2
         if ((A(3) .le. A(1)) .and. (A(3) .le. A(2))) cax = 3
         node_cax(node_in) = cax
      else
         cax = node_cax(node_in)
      end if

      !cut location -------------------------------------------------------------------------------------------------------------
      i = gridind_in(cax, 1) - 1
      n_p = 0
      pincol = 0
      np_per_node = int(ceiling(real(nprocs_in)/2)/real(nprocs_in)*ntotal_in)
      do while (n_p .lt. np_per_node)
         i = i + 1
         if (cax .eq. 1) then
            pincol = SUM(pincell_ORB( &
                         i, &
                         gridind_in(2, 1):gridind_in(2, 2), &
                         gridind_in(3, 1):gridind_in(3, 2)))
         elseif (cax .eq. 2) then
            pincol = SUM(pincell_ORB( &
                         gridind_in(1, 1):gridind_in(1, 2), &
                         i, &
                         gridind_in(3, 1):gridind_in(3, 2)))
         else
            pincol = SUM(pincell_ORB( &
                         gridind_in(1, 1):gridind_in(1, 2), &
                         gridind_in(2, 1):gridind_in(2, 2), &
                         i))
         end if
         n_p = n_p + pincol
      end do

      if ((np_per_node - (n_p - pincol) .lt. n_p - np_per_node)) then
         i = i - 1
         n_p = n_p - pincol
      end if

      !saving output information ------------------------------------------------------------------------------------------------
      procrange_lo(1) = procrange_in(1)
      procrange_lo(2) = procrange_in(1) + ceiling(real(nprocs_in)/2) - 1
      procrange_hi(1) = procrange_lo(2) + 1
      procrange_hi(2) = procrange_lo(1) + nprocs_in - 1
      gridind_out(:, 1) = gridind_in(:, 1)
      gridind_out(:, 2) = gridind_in(:, 2)

      if (procid .le. procrange_lo(2)) then
         node_out = 2*node_in
         procrange_out = procrange_lo
         ntotal_out = n_p
         gridind_out(cax, 1) = gridind_in(cax, 1)
         gridind_out(cax, 2) = i
      else
         node_out = 2*node_in + 1
         procrange_out = procrange_hi
         ntotal_out = ntotal_in - n_p
         gridind_out(cax, 1) = i + 1
         gridind_out(cax, 2) = gridind_in(cax, 2)
      end if

      nprocs_out = procrange_out(2) - procrange_out(1) + 1

      !travelling to child node/saving boundary information ---------------------------------------------------------------------
      if (nprocs_out .ne. 1) then
         bounds_glob = ORB_bounds(procid, numprocs, MPI_ftype, scale_k, repartition_mode, gridind_out, nprocs_out, node_out, &
                                  procrange_out, &
                                  ntotal_out, pincell_ORB, dcell, mingridx_in, maxgridx_in)
      else

         bounds_out(1:dim) = mingridx_in(:) + (gridind_out(:, 1) - 1)*dcell
         bounds_out(dim + 1:2*dim) = mingridx_in(:) + gridind_out(:, 2)*dcell

         if (repartition_mode .eq. 3) call potential_neighbour_process_search(node_out)

         if (node_cut(1) .eq. 0) bounds_out(4) = bounds_out(4) + bound_extend*scale_k*hsml
         if (node_cut(2) .eq. 0) bounds_out(1) = bounds_out(1) - bound_extend*scale_k*hsml
         if (node_cut(3) .eq. 0) bounds_out(5) = bounds_out(5) + bound_extend*scale_k*hsml
         if (node_cut(4) .eq. 0) bounds_out(2) = bounds_out(2) - bound_extend*scale_k*hsml
         if (node_cut(5) .eq. 0) bounds_out(6) = bounds_out(6) + bound_extend*scale_k*hsml
         if (node_cut(6) .eq. 0) bounds_out(3) = bounds_out(3) - bound_extend*scale_k*hsml

         call MPI_ALLGATHER(bounds_out, 2*dim, MPI_ftype, bounds_glob, 2*dim, MPI_ftype, MPI_COMM_WORLD, ierr)

      end if

   end function ORB_bounds

   !==============================================================================================================================
   subroutine P_trim(gridind_in, ngridx_trim, pincell_ORB)
      ! Trims particle-in-cell grid so as to obtain minimal bounding boxes to obtain accurate cut axis orientations

      implicit none
      integer, intent(in):: gridind_in(dim, 2), pincell_ORB(:, :, :)
      integer, intent(out):: ngridx_trim(dim)
      integer:: i, j, k, newi(2), newj(2), newk(2), oldi(2), oldj(2), oldk(2)

      oldi(:) = gridind_in(1, :)
      oldj(:) = gridind_in(2, :)
      oldk(:) = gridind_in(3, :)
      newi(:) = gridind_in(1, :)
      newj(:) = gridind_in(2, :)
      newk(:) = gridind_in(3, :)

      !Trimming input grid ------------------------------------------------------------------------------------------------------
      !reducing x-length of grid
      !finding new start x-index
      newstarti: do i = oldi(1), oldi(2)
         if (any(pincell_ORB(i, oldj(1):oldj(2), oldk(1):oldk(2)) .ne. 0)) then
            newi(1) = i
            exit newstarti
         end if
      end do newstarti

      !finding new end x-index
      newendi: do i = oldi(2), newi(1), -1
         if (any(pincell_ORB(i, oldj(1):oldj(2), oldk(1):oldk(2)) .ne. 0)) then
            newi(2) = i
            exit newendi
         end if
      end do newendi

      !reducing y-length of grid
      !finding new start y-index
      newstartj: do j = oldj(1), oldj(2)
         if (any(pincell_ORB(newi(1):newi(2), j, oldk(1):oldk(2)) .ne. 0)) then
            newj(1) = j
            exit newstartj
         end if
      end do newstartj

      !finding new end y-index
      newendj: do j = oldj(2), newj(1), -1
         if (any(pincell_ORB(newi(1):newi(2), j, oldk(1):oldk(2)) .ne. 0)) then
            newj(2) = j
            exit newendj
         end if
      end do newendj

      !reducing z-length of grid
      !finding new start z-index
      newstartk: do k = oldk(1), oldk(2)
         if (any(pincell_ORB(newi(1):newi(2), oldj(1):oldj(2), k) .ne. 0)) then
            newk(1) = k
            exit newstartk
         end if
      end do newstartk

      !finding new end z-index
      newendk: do k = oldk(2), newk(1), -1
         if (any(pincell_ORB(newi(1):newi(2), oldj(1):oldj(2), k) .ne. 0)) then
            newk(2) = k
            exit newendk
         end if
      end do newendk

      ngridx_trim(1) = newi(2) - newi(1) + 1
      ngridx_trim(2) = newj(2) - newj(1) + 1
      ngridx_trim(3) = newk(2) - newk(1) + 1

   end subroutine P_trim

   !==============================================================================================================================
   subroutine potential_neighbour_process_search(ID_node)
      ! Used to determine which nodes are located at the edge of the domain in each direction.
      ! IGNORE BELOW
      ! Finds any potential neighbours by exploiting the binary tree data structure
      ! First finds nodes in the tree that define the edge of the local process. Then finds processes that also use that node as an
      ! edge but on the opposite side. E.g. if local process has node 2 as an east edge, potential east neighbours of the local
      ! process are the processes that use node 2 as a west edge.

      implicit none
      integer:: ID_node

      !finding nodes that define edge for current process
      node_cut(:) = 0
      do while (ID_node .ne. 1)

         if (mod(ID_node, 2) .eq. 0) then
            ID_node = ID_node/2
            if ((node_cax(ID_node) .eq. 1) .and. (node_cut(1) .eq. 0)) node_cut(1) = ID_node !right face of process defined by 'parent'
            if ((node_cax(ID_node) .eq. 2) .and. (node_cut(3) .eq. 0)) node_cut(3) = ID_node !top face of process defined by 'parent'
            if ((node_cax(ID_node) .eq. 3) .and. (node_cut(5) .eq. 0)) node_cut(5) = ID_node !top face of process defined by 'parent'
         else
            ID_node = (ID_node - 1)/2
            if ((node_cax(ID_node) .eq. 1) .and. (node_cut(2) .eq. 0)) node_cut(2) = ID_node !left of process defined by 'parent'
            if ((node_cax(ID_node) .eq. 2) .and. (node_cut(4) .eq. 0)) node_cut(4) = ID_node !bottom of process defined by 'parent'
            if ((node_cax(ID_node) .eq. 3) .and. (node_cut(6) .eq. 0)) node_cut(6) = ID_node !bottom face of process defined by 'parent'
         end if

      end do

   end subroutine potential_neighbour_process_search

   !==============================================================================================================================
   subroutine subdomain_neighbour(procid, numprocs, bounds_glob, scale_k, n_process_neighbour)
      !creates list of adjacent subdomains for the local subdomain by searching potential neighbours and seeing if they overlap

      implicit none
      integer, intent(in):: procid, numprocs
      real(f), intent(in):: bounds_glob(2*dim, numprocs), scale_k
      integer, intent(out):: n_process_neighbour
      integer:: pid
      real(f):: bounds_loc_min(dim), bounds_loc_max(dim), bounds_rem_min(dim), bounds_rem_max(dim)

      !Checking if potential neighbours overlap with local process. Overlap = adjacent.
      n_process_neighbour = 0
      bounds_loc_min(1:dim) = bounds_glob(1:dim, procid + 1) - 0.5_f*scale_k*hsml
      bounds_loc_max(1:dim) = bounds_glob(dim + 1:2*dim, procid + 1) + 0.5_f*scale_k*hsml
      neighboursearch: do pid = 1, numprocs
         if (pid .ne. procid + 1) then
            bounds_rem_min(:) = bounds_glob(1:dim, pid) - 0.5_f*scale_k*hsml
            bounds_rem_max(:) = bounds_glob(dim + 1:2*dim, pid) + 0.5_f*scale_k*hsml
            ! if local and remote process' extended boundaries don't overlap, check next process
            if (any([bounds_rem_max(:) < bounds_loc_min(:), bounds_rem_min(:) > bounds_loc_max(:)])) cycle neighboursearch
            n_process_neighbour = n_process_neighbour + 1
            neighbours(n_process_neighbour)%pid = pid - 1
            neighbours(n_process_neighbour)%bounds(:) = bounds_glob(:, pid)
         end if
      end do neighboursearch

   end subroutine subdomain_neighbour

end module ORB_m
