module ORB_m

   use datatypes, only: particles, time_tracking, system_clock_timer
   use param, only: f, dim, hsml
   use param_para, only: neighbour_data, partition_tracking, dcell_ORB, box_ratio_threshold

   private
   real(f):: box_ratio_previous(dim, dim) = TINY(1._f), bounds_loc(6)
   type(neighbour_data), allocatable:: neighbours(:)
   type(partition_tracking):: partition_track
   integer:: prev_load, n_process_neighbour, repartition_mode
   integer, allocatable:: node_cax(:), pincell_ORB(:, :, :)[:], cellmins(:,:)[:], cellmaxs(:,:)[:]
   real(f), parameter:: dcell = hsml*dcell_ORB

   public:: ORB, partition_track, neighbours, n_process_neighbour

contains

   subroutine ORB(itimestep, thisImage, numImages, scale_k, ntotal, ntotal_loc, nhalo_loc, nvirt_loc, parts, timings)

      use param_para, only: ORBcheck1, ORBcheck2

      implicit none
      integer, intent(in):: thisImage, numImages, itimestep, ntotal
      real(f), intent(in):: scale_k
      type(particles), intent(inout):: parts(:)
      type(time_tracking), intent(inout):: timings
      integer, intent(inout):: ntotal_loc, nhalo_loc, nvirt_loc
      integer:: tree_layers, maxnode, stepsSincePrevious, ngridx(3), i, j, k, d, imageRange_ini(2), gridind_ini(3, 2)
      real(f):: mingridx_ini(3), maxgridx_ini(3), current_to_previous(dim, dim), box_ratio_current(dim, dim)
      real(f), allocatable:: bounds_glob(:, :)[:]

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

         if ((repartition_mode > 1) .or. (itimestep .eq. 1) ) then

            partition_track%n_parts = partition_track%n_parts + 1
            if (itimestep /= 1) then
               partition_track%mintstep_bn_part = min(partition_track%mintstep_bn_part, stepsSincePrevious)
               partition_track%maxtstep_bn_part = max(partition_track%maxtstep_bn_part, stepsSincePrevious)
            end if
            partition_track%prev_part_tstep = itimestep

            ! Creating particle-in-cell grid
            allocate( bounds_glob(2*dim, numImages)[*], cellmins(3, numImages)[*], cellmaxs(3, numImages)[*])
            call particle_grid(thisImage, numImages, ntotal_loc, parts(1:ntotal_loc), ngridx, mingridx_ini, maxgridx_ini, &
                               pincell_ORB, cellmins, cellmaxs)

            do d = 1,dim
               box_ratio_current(:, d) = DBLE(ngridx(d))/DBLE(ngridx(:))
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

            ! determine subdomain boundaries using particle distribution
            gridind_ini(:, 1) = 1
            gridind_ini(:, 2) = ngridx(:)
            imageRange_ini(1) = 1
            imageRange_ini(2) = numImages

            call ORB_bounds(thisImage, numImages, scale_k, gridind_ini, numImages, 1, imagerange_ini, ntotal, &
               dcell, mingridx_ini, maxgridx_ini, bounds_glob)

            deallocate( bounds_glob, cellmins, cellmaxs)

         end if

      end if

   end subroutine ORB

   !====================================================================================================================
   subroutine particle_grid(thisImage, numImages, ntotal_loc, parts, ngridx, mingridx, maxgridx, pincell_ORB, &
      cellmins, cellmaxs)

      implicit none
      integer, intent(in):: thisImage, numImages, ntotal_loc
      type(particles), intent(in):: parts(:)
      integer, intent(out):: ngridx(3), cellmins(3, numImages)[*], cellmaxs(3, numImages)[*]
      real(f), intent(out):: mingridx(3), maxgridx(3)
      integer, allocatable, intent(inout):: pincell_ORB(:, :, :)[:]
      real(f):: minx(3), maxx(3)
      integer:: i, j, k, d, cellrange(3), icell(3), n_nonzerocells, n, otherImage

      ! Local max, min, in each direction ------------------------------------------------------------------------------
      minx(1:dim) = parts(1)%x(1:dim)
      maxx(1:dim) = parts(1)%x(1:dim)
      do i = 2, ntotal_loc
         do d = 1, dim
            if (parts(i)%x(d) > maxx(d)) maxx(d) = parts(i)%x(d)
            if (parts(i)%x(d) < minx(d)) minx(d) = parts(i)%x(d)
         end do
      end do

      ! Global max, min, in each direction -----------------------------------------------------------------------------
      mingridx(:) = minx(:) - 10*dcell
      maxgridx(:) = maxx(:) + 10*dcell
      call co_max(maxgridx)
      call co_min(mingridx)

      ! Number of grid cells and adjusting max extent in each direction ------------------------------------------------
      ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
      maxgridx(:) = mingridx(:) + dcell*ngridx(:)

      ! Reducing search area by locating indices in which local particles are contained within
      cellmins(:, thisImage) = int((minx(:) - mingridx(:))/dcell) + 1
      cellmaxs(:, thisImage) = int((maxx(:) - mingridx(:))/dcell) + 1
      cellrange(:) = cellmaxs(:, thisImage) - cellmins(:, thisImage) + 1

      ! Creating list of non-zero grid cells ------------------------------------------------------------------------------------
      call co_max(cellrange)
      allocate( pincell_ORB(cellrange(1), cellrange(2), cellrange(3))[*] )

      pincell_ORB(:, :, :) = 0

      ! Counting local particles in each cell, and tracking non-zero cells
      n_nonzerocells = 0
      do i = 1,ntotal_loc
         icell(1:dim) = int((parts(i)%x(:) - mingridx(:))/dcell) + 1
         if (dim==2) icell(3) = 1
         icell(:) = icell(:) - cellmins(:, thisImage) + 1
         pincell_ORB(icell(1), icell(2), icell(3)) = pincell_ORB(icell(1), icell(2), icell(3)) + 1
      end do

      ! Getting each image to share it's min-max cell indices to other images
      do i = 1, numImages - 1
         otherImage = thisImage + i
         if (otherImage > numImages) otherImage = otherImage - numImages
         cellmins(:, thisImage)[otherImage] = cellmins(:, thisImage)[thisImage]
         cellmaxs(:, thisImage)[otherImage] = cellmaxs(:, thisImage)[thisImage]
      end do
      
      ! ensuring all images leave with full pincell_ORB
      sync all

   end subroutine particle_grid

   !====================================================================================================================
   recursive subroutine ORB_bounds(thisImage, numImages, scale_k, gridind_in, numImages_in, node_in, imagerange_in, &
      ntotal_in, dcell, mingridx_in, maxgridx_in, bounds_glob)

      ! Recursive function that performs the 'bisection' part of the ORB algorithm

      use param_para, only: bound_extend

      implicit none
      integer, intent(in):: thisImage, numImages, gridind_in(3, 2), numImages_in, node_in, &
         imagerange_in(2), ntotal_in
      real(f), intent(in):: mingridx_in(3), maxgridx_in(3), dcell, scale_k
      real(f), intent(out):: bounds_glob(:, :)[*]
      integer:: i, j, k, d, ngridx_trim(3), A(3), cax, n_p, pincol, np_per_node, otherImage_limits(2, 3), otherImage, &
         n, imagerange_out(2), gridind_out(3, 2), ntotal_out, numImages_out, node_out, imagerange_lo(2), &
         imagerange_hi(2)

      !determining cut axis. 1 = x, 2 = y ------------------------------------------------------------------------------
      if (repartition_mode == 3) then
         call P_trim(thisImage, numImages, gridind_in, ngridx_trim)
         A(1) = ngridx_trim(2)*ngridx_trim(3)
         A(2) = ngridx_trim(1)*ngridx_trim(3)
         A(3) = ngridx_trim(1)*ngridx_trim(2)
         cax = minloc(A, 1)
         node_cax(node_in) = cax
      else
         cax = node_cax(node_in)
      end if
      
      !cut location ----------------------------------------------------------------------------------------------------
      i = gridind_in(cax, 1) - 1
      n_p = 0
      np_per_node = int(ceiling(real(numImages_in)/2)/real(numImages_in)*ntotal_in)
      do while (n_p < np_per_node)
         i = i + 1
         pincol = 0
         do n = 0, numImages-1
            otherImage = mod(thisImage+n, numImages)
            if (otherImage == 0) otherImage = numImages
            if (i >= cellmins(cax, otherImage) .and. i <= cellmaxs(cax, otherImage)) then

               do d = 1, 3
                  otherImage_limits(1, d) = max(1, gridind_in(d, 1) - cellmins(d, otherImage) + 1)
                  otherImage_limits(2, d) = min(cellmaxs(d, otherImage), gridind_in(d, 2)) - cellmins(d, otherImage) + 1
               end do
            
               otherImage_limits(1:2, cax) = i - cellmins(cax, otherImage) + 1
               pincol = pincol + sum(pincell_ORB(otherImage_limits(1,1):otherImage_limits(2, 1), &
                                        otherImage_limits(1,2):otherImage_limits(2, 2), &
                                        otherImage_limits(1,3):otherImage_limits(2, 3))[otherImage])

            end if
         end do
         n_p = n_p + pincol
      end do
      
      ! extra step to nudge cut location backwards if that's a better position
      if ((np_per_node - (n_p - pincol) < n_p - np_per_node)) then
         i = i - 1
         n_p = n_p - pincol
      end if
      
      !saving output information ---------------------------------------------------------------------------------------
      imagerange_lo(1) = imagerange_in(1)
      imagerange_lo(2) = imagerange_in(1) + ceiling(real(numImages_in)/2) - 1
      imagerange_hi(1) = imagerange_lo(2) + 1
      imagerange_hi(2) = imagerange_lo(1) + numImages_in - 1
      gridind_out(:, 1) = gridind_in(:, 1)
      gridind_out(:, 2) = gridind_in(:, 2)

      if (thisImage <= imagerange_lo(2)) then
         node_out = 2*node_in
         imagerange_out(:) = imagerange_lo(:)
         ntotal_out = n_p
         gridind_out(cax, 1) = gridind_in(cax, 1)
         gridind_out(cax, 2) = i
      else
         node_out = 2*node_in + 1
         imagerange_out(:) = imagerange_hi(:)
         ntotal_out = ntotal_in - n_p
         gridind_out(cax, 1) = i + 1
         gridind_out(cax, 2) = gridind_in(cax, 2)
      end if

      numImages_out = imagerange_out(2) - imagerange_out(1) + 1

      if (numImages_out /= 1) then
         call ORB_bounds(thisImage, numImages, scale_k, gridind_out, numImages_out, node_out, imagerange_out, &
            ntotal_out, dcell, mingridx_in, maxgridx_in, bounds_glob)
      else
         sync all
      end if

   end subroutine ORB_bounds

   !====================================================================================================================
   subroutine P_trim(thisImage, numImages, gridind_in, ngridx_trim)
   ! Trims particle-in-cell grid so as to obtain minimal bounding boxes to obtain accurate cut axis orientations

      implicit none
      integer, intent(in):: thisImage, numImages, gridind_in(3, 2)
      integer, intent(out):: ngridx_trim(3)
      integer:: i, j, k, n, newi(2), newj(2), newk(2), oldi(2), oldj(2), oldk(2), inspectrange(2), otherImage

      newi(:) = gridind_in(1, :)
      newj(:) = gridind_in(2, :)
      newk(:) = gridind_in(3, :)

      !Trimming input grid ------------------------------------------------------------------------------------------------------
      !reducing x-length of grid
      !finding new start x-index
      newi(1) = trim_helper(thisImage, numImages, 1, newi(1), newi(2), 2, newj, 3, newk)

      !finding new end x-index
      newi(2) = trim_helper(thisImage, numImages, 1, newi(2), newi(1), 2, newj, 3, newk)

      !reducing y-length of grid
      !finding new start y-index
      newj(1) = trim_helper(thisImage, numImages, 2, newj(1), newj(2), 1, newi, 3, newk)

      !finding new end y-index
      newj(2) = trim_helper(thisImage, numImages, 2, newj(2), newj(1), 1, newi, 3, newk)

      !reducing z-length of grid
      !finding new start z-index
      newk(1) = trim_helper(thisImage, numImages, 3, newk(1), newk(2), 1, newi, 2, newj)

      !finding new end z-index
      newk(2) = trim_helper(thisImage, numImages, 3, newk(2), newk(1), 1, newi, 2, newj)

      ngridx_trim(1) = newi(2) - newi(1) + 1
      ngridx_trim(2) = newj(2) - newj(1) + 1
      ngridx_trim(3) = newk(2) - newk(1) + 1

   end subroutine P_trim

   function trim_helper(thisImage, numImages, main_axis, main_start, main_end, off_axis1, off_limits1, off_axis2, &
         off_limits2) result(newindex)

      integer, intent(in):: thisImage, numImages, main_axis, main_start, main_end, off_axis1, off_limits1(2), &
         off_axis2, off_limits2(2)
      integer:: newindex, step, a_main, a_off1, a_off2, n, icell(3), otherImage, otherImage_limits(2, 3), otherImage_maxincell

      newindex = main_start

      step = sign(1, main_end - main_start)
      do a_main = main_start, main_end, step
         do n = 0, numImages-1
            otherImage = mod(thisImage+n, numImages)
            if (otherImage == 0) otherImage = numImages
            if (.not.(a_main > cellmaxs(main_axis, otherImage) .or. &
                      off_limits1(1) > cellmaxs(off_axis1, otherImage) .or. &
                      off_limits2(1) > cellmaxs(off_axis2, otherImage) .or. &
                      cellmins(main_axis, otherImage) > a_main .or. &
                      cellmins(off_axis1, otherImage) > off_limits1(2) .or. &
                      cellmins(off_axis2, otherImage) > off_limits2(2))) then
               otherImage_limits(1:2, main_axis) = a_main - cellmins(main_axis, otherImage) + 1
               otherImage_limits(1, off_axis1) = max(cellmins(off_axis1, otherImage), off_limits1(1)) - &
                  cellmins(off_axis1, otherImage) + 1
               otherImage_limits(2, off_axis1) = min(cellmaxs(off_axis1, otherImage), off_limits1(2)) - &
                  cellmins(off_axis1, otherImage) + 1
               otherImage_limits(1, off_axis2) = max(cellmins(off_axis2, otherImage), off_limits2(1)) - &
                  cellmins(off_axis2, otherImage) + 1
               otherImage_limits(2, off_axis2) = min(cellmaxs(off_axis2, otherImage), off_limits2(2)) - &
                  cellmins(off_axis2, otherImage) + 1
               otherImage_maxincell = maxval(pincell_ORB(otherImage_limits(1, 1):otherImage_limits(2, 1), &
                                                      otherImage_limits(1, 2):otherImage_limits(2, 2), &
                                                      otherImage_limits(1, 3):otherImage_limits(2, 3))[otherImage])
               if (otherImage_maxincell > 0) Then
                  newindex = a_main
                  return
               end if
            end if
         end do
      end do

   end function trim_helper

end module ORB_m

