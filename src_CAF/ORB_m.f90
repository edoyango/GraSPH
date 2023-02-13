module ORB_m

   use datatypes, only: particles, time_tracking, system_clock_timer
   use param, only: f, dim, hsml
   use param_para, only: neighbour_data, partition_tracking, dcell_ORB

   private
   type(neighbour_data), allocatable:: neighbours(:)
   type(partition_tracking):: partition_track
   integer:: prev_load
   integer, allocatable:: node_cax(:)
   real(f), parameter:: dcell = hsml*dcell_ORB

   public:: ORB

contains

   subroutine ORB(itimestep, thisImage, numImages, scale_k, ntotal, ntotal_loc, nhalo_loc, nvirt_loc, parts, timings)

      use param_para, only: ORBcheck1, ORBcheck2

      implicit none
      integer, intent(in):: thisImage, numImages, itimestep, ntotal
      real(f), intent(in):: scale_k
      type(particles), intent(inout):: parts(:)
      type(time_tracking), intent(inout):: timings
      integer, intent(inout):: ntotal_loc, nhalo_loc, nvirt_loc
      integer:: tree_layers, maxnode, repartition_mode, stepsSincePrevious, ngridx(dim)
      real(f):: mingridx_ini(dim), maxgridx_ini(dim)
      integer, allocatable:: pincell_ORB(:, :, :)[:]

      !allocating partitioning arrays and initialising diagnostic variables ------------------------------------------
      timings%t_ORB = timings%t_ORB - system_clock_timer() ! commence timing of ORB algoirthm
      if (itimestep == 0) then
         tree_layers = ceiling(log(dble(numImages))/log(2.d0))
         maxnode = 2*2**tree_layers - 1
         allocate (node_cax(maxnode), neighbours(numImages))
      end if

      ! Boundary Determiniation Algorithm ----------------------------------------------------------------------------
      repartition_mode = 1 !initially assumes no partition
      ! only checks if boundary needs updating every 50-100 time-steps
      stepsSincePrevious = itimestep - partition_track%prev_part_tstep
      if ((itimestep == 0) .or. &
          ((stepsSincePrevious >= ORBcheck1) .and. &
           (mod(stepsSincePrevious, ORBcheck2) == 0))) then

         ! checking if change in particles on current image is > 5%
         if (ntotal_loc > prev_load + 0.05_f*DBLE(ntotal)/DBLE(thisImage)) then
            repartition_mode = 2
         end if

         call co_max(repartition_mode)

         sync all

         if ((repartition_mode > 1) .or. (itimestep .eq. 1) ) then

            partition_track%n_parts = partition_track%n_parts + 1
            if (itimestep /= 1) then
               partition_track%mintstep_bn_part = min(partition_track%mintstep_bn_part, stepsSincePrevious)
               partition_track%maxtstep_bn_part = max(partition_track%maxtstep_bn_part, stepsSincePrevious)
            end if
            partition_track%prev_part_tstep = itimestep

            ! Creating particle-in-cell grid
            call particle_grid(thisImage, numImages, ntotal_loc, parts(1:ntotal_loc), ngridx, mingridx_ini, maxgridx_ini, &
                               pincell_ORB)

         end if

      end if

   end subroutine ORB

   !====================================================================================================================
   subroutine particle_grid(thisImage, numImages, ntotal_loc, parts, ngridx, mingridx, maxgridx, pincell_ORB)

      implicit none
      integer, intent(in):: thisImage, numImages, ntotal_loc
      type(particles), intent(in):: parts(:)
      integer, intent(out):: ngridx(:)
      real(f), intent(out):: mingridx(:), maxgridx(:)
      integer, allocatable, intent(inout):: pincell_ORB(:, :, :)[:]
      real(f):: minx(3), maxx(3)
      integer:: i, j, k, d, cellmins(3), cellmaxs(3), cellrange(3), Plist_size, icell(3), n_nonzerocells, n, otherImage
      integer, allocatable:: Plist_loc(:, :)

      ! Local max, min, in each direction ------------------------------------------------------------------------------
      minx(1:dim) = parts(1)%x(1:dim)
      maxx(1:dim) = parts(1)%x(1:dim)
      do i = 2, ntotal_loc
         do d = 1, dim
            if (parts(i)%x(d) .gt. maxx(d)) maxx(d) = parts(i)%x(d)
            if (parts(i)%x(d) .lt. minx(d)) minx(d) = parts(i)%x(d)
         end do
      end do

      ! Global max, min, in each direction -----------------------------------------------------------------------------
      mingridx(:) = minx(:)
      maxgridx(:) = maxx(:)
      call co_max(maxgridx)
      call co_min(mingridx)

      sync all
      
      mingridx(:) = mingridx(:) - dcell
      maxgridx(:) = maxgridx(:) + dcell

      ! Number of grid cells and adjusting max extent in each direction ------------------------------------------------
      ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
      maxgridx(:) = mingridx(:) + dcell*ngridx(:)

      ! Reducing search area by locating indices in which local particles are contained within
      cellmins(:) = int((minx(:) - mingridx(:))/dcell) + 1
      cellmaxs(:) = int((maxx(:) - mingridx(:))/dcell) + 1
      cellrange(:) = cellmaxs(:) - cellmins(:) + 1

      ! Creating list of non-zero grid cells ------------------------------------------------------------------------------------
      Plist_size = MIN(ntotal_loc, PRODUCT(cellrange(:)))
      allocate( pincell_ORB(ngridx(1), ngridx(2), ngridx(3))[*], &
                Plist_loc(4, Plist_size))

      pincell_ORB(:, :, :) = 0

      ! Counting local particles in each cell, and tracking non-zero cells
      n_nonzerocells = 0
      do i = 1,ntotal_loc
         icell(1:dim) = int((parts(i)%x(:) - mingridx(:))/dcell) + 1
         if (dim==2) icell(3) = 1
         if (pincell_ORB(icell(1), icell(2), icell(3)) == 0) then
            n_nonzerocells = n_nonzerocells + 1
            Plist_loc(1:3, n_nonzerocells) = icell(:)
         end if
         pincell_ORB(icell(1), icell(2), icell(3)) = pincell_ORB(icell(1), icell(2), icell(3)) + 1
      end do

      ! populating local list of non-zero cells with values
      do i = 1, n_nonzerocells
         icell(:) = Plist_loc(1:3, i)
         Plist_loc(4, i) = pincell_ORB(icell(1), icell(2), icell(3))
      end do

      ! Each image populates other image's pincell
      do n = 1, numImages-1
         otherImage = thisImage + n
         if (otherImage > numImages) otherImage = otherImage - numImages
         sync all
         do i = 1, n_nonzerocells
            icell(:) = Plist_loc(1:3, i)
            pincell_ORB(icell(1), icell(2), icell(3))[otherImage] = &
               pincell_ORB(icell(1), icell(2), icell(3))[otherImage] + Plist_loc(4, i)
         end do
      end do

      ! final sync so images leave subroutine together
      sync all

   end subroutine particle_grid

end module ORB_m

