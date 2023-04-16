module flink_list_m

   use datatypes, only: particles, interactions
   use iso_fortran_env, only: atomic_int_kind
   use param, only: dim, hsml, f, halotype, index_to_imageindex, imageindex_to_index

   !use error_msg_m, only: error_msg
   use kernel_m, only: kernel

   private
   public:: flink_list

contains

   !==============================================================================================================================
   subroutine flink_list(thisImage, numImages, maxinter, scale_k, maxn, ntotal_loc, nvirt_loc, niac, parts, pairs, &
         nexti)
      ! save as above, but for 3D

      implicit none
      integer, intent(in):: thisImage, numImages, maxn, maxinter, ntotal_loc, nvirt_loc
      real(f), intent(in):: scale_k
      type(particles), intent(in), codimension[*]:: parts(:)
      integer, intent(out):: niac, nexti(:)
      type(interactions), intent(out):: pairs(:)
      integer, parameter:: maxpcell = 125
      integer:: ngridx(3), jth, i, j, k, d, xi, yi, zi, ierr, ncells_loc_avg, ncells, global_cellhash, &
         grid_image_coords(2), global_i, local_i(2), remote_j(2), icell(3), icellhash, jcellhash
      real(f):: mingridx(3), maxgridx(3), dcell
      integer(ATOMIC_INT_KIND):: tmp
      integer(ATOMIC_INT_KIND), allocatable:: pincell(:)[:]
      integer, allocatable:: gridcellhash(:), cells(:, :)[:]
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

      !Determining bounding box extents
      mingridx(:) = parts(1)%x(:)
      maxgridx(:) = parts(1)%x(:)
      do i = 2, ntotal_loc + nvirt_loc
         do d = 1, dim
            mingridx(d) = MIN(mingridx(d), parts(i)%x(d))
            maxgridx(d) = MAX(maxgridx(d), parts(i)%x(d))
         end do
      end do

      call co_max(maxgridx)
      call co_min(mingridx)

      !Determining number of grid cells in each direction
      dcell = scale_k*hsml
      maxgridx(:) = maxgridx(:) + 2._f*dcell
      mingridx(:) = mingridx(:) - 2._f*dcell
      ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
      maxgridx(:) = mingridx(:) + ngridx(:)*dcell

      ncells = product(ngridx)
      ncells_loc_avg = ceiling(real(ncells, f)/numImages)
      allocate(pincell(ncells_loc_avg)[*], &
               cells(maxpcell, ncells_loc_avg)[*], &
               gridcellhash(ntotal_loc+nvirt_loc))

      pincell(:) = 0

      local_i(2) = thisImage
      tmp = -1
      do i = 1, ntotal_loc + nvirt_loc
         icell(:) = int((parts(i)%x(:) - mingridx(:))/dcell) + 1
         global_cellhash = cells2Hash(icell, ngridx)
         gridcellhash(i) = global_cellhash
         grid_image_coords = globalhash2imagehash(global_cellhash, numImages, ncells)
         
         call atomic_fetch_add(pincell(grid_image_coords(1))[grid_image_coords(2)], 1, tmp)
         tmp = tmp + 1 
         
         local_i(1) = i
         global_i = imageindex_to_index(maxn, numImages, local_i)
         cells(tmp, grid_image_coords(1))[grid_image_coords(2)] = global_i
         ! pincell(icell, jcell, kcell) = pincell(icell, jcell, kcell) + 1
         ! cells(pincell(icell, jcell, kcell), icell, jcell, kcell) = i
      end do

      sync all

      niac = 0
      ierr = 0
      nexti(1) = 1
      local_i(2) = thisImage
      do i = 1, ntotal_loc + nvirt_loc
         ! write(*,*) thisImage, i
         local_i(1) = i
         global_i = imageindex_to_index(maxn, numimages, local_i)
         ! Retrieving particle i's grid cell indices
         icellhash = gridcellhash(i)
         grid_image_coords = globalhash2imagehash(icellhash, numImages, ncells)
         ! finding pairs within cell icell
         do j = 1, pincell(grid_image_coords(1))[grid_image_coords(2)]
            jth = cells(j, grid_image_coords(1))[grid_image_coords(2)]            
            ! if (jth > global_i) then
               ! remote_j = index_to_imageindex(maxn, numImages, jth)
               ! call check_if_interact(maxinter, scale_k, i, jth, parts(i), parts(remote_j(1))[remote_j(2)], niac, pairs, ierr)
            ! end if
         end do

      !    ! finding pairs within cells adjacent to i's cell
         do k = 1, 13
            jcellhash = icellhash + dhash(sweep(:, k), ngridx)
            grid_image_coords = globalhash2imagehash(jcellhash, numImages, ncells)
      !       xi = icell + sweep(1, k)
      !       yi = jcell + sweep(2, k)
      !       zi = kcell + sweep(3, k)
      !       do j = 1, pincell(grid_image_coords(1))[grid_image_coords(2)]
      !          jth = cells(j, grid_image_coords(1))[grid_image_coords(2)] 
      !          remote_j = index_to_imageindex(maxn, numImages, jth)
      !          ! call check_if_interact(maxinter, scale_k, i, jth, parts(i), parts(remote_j(1))[remote_j(2)], niac, pairs, ierr)
      ! !          call check_if_interact(maxinter, scale_k, i, jth, parts(i), parts(jth), niac, pairs, ierr)
      !       end do
         end do

         nexti(i + 1) = niac + 1

      end do

      ! if (ierr == 1) then
      !    print *, ' >>> Error <<< : Too many interactions'
      !    stop
      ! end if

   end subroutine flink_list

   !==============================================================================================================================
   pure subroutine check_if_interact(maxinter, scale_k, i, j, p_i, p_j, niac, pairs, ierr)
      ! subroutine to chekc if two particles are interacting and consequently adding to pair list

      implicit none
      integer, intent(in):: maxinter, i, j
      real(f), intent(in):: scale_k
      type(particles), intent(in):: p_i, p_j
      integer, intent(inout):: niac
      type(interactions), intent(inout):: pairs(:)
      integer, intent(inout):: ierr
      real(f):: dxiac(dim), r

      ! only consider interactions between:
      ! real-real
      ! real-virtual
      ! real-halo
      ! halo-virtual - for deterimining boundary
      ! if (.not.((p_i%itype<0 .and. p_j%itype<0) .or. (p_i%itype>halotype .and. p_j%itype>halotype))) then
      if (p_i%itype>0 .or. p_j%itype>0) then
         dxiac(:) = p_i%x(:) - p_j%x(:)
         r = SQRT(SUM(dxiac*dxiac))
         if (r < hsml*scale_k) then
            niac = niac + 1
            if (niac < maxinter) then
               pairs(niac)%j = j
               call kernel(r, dxiac, hsml, pairs(niac)%w, pairs(niac)%dwdx(:))
            else
               ierr = 1
            end if
         end if
      end if

   end subroutine check_if_interact

   pure function cells2Hash(cellCoords, ncells) result(hash)
     
      implicit none
      integer, intent(in):: cellCoords(3), ncells(3)
      integer:: hash
     
      hash = ncells(1)*ncells(2)*(cellCoords(3) - 1) + ncells(1)*(cellCoords(2) - 1) + cellCoords(1)
     
   end function cells2Hash

   pure function hash2Cell(hash, ncells) result(cellCoords)
   
      implicit none
      integer, intent(in):: hash, ncells(3)
      integer:: cellCoords(3), nxy, modhashnxy
     
      nxy = ncells(1)*ncells(2)
      modhashnxy = mod(hash-1,nxy)
     
      cellCoords(3) = int((hash-1)/nxy) + 1
      cellCoords(2) = modhashnxy/ncells(1) + 1
      cellCoords(1) = mod(modhashnxy,ncells(1)) + 1
     
   end function hash2Cell

   integer pure function dhash(dcells, ncells)
      ! h = nx*ny*(k-1) + nx*(j-1) + i = nx*ny*k - nx*ny + nx*j - nx + i
      ! dhdi = 1
      ! dhdj = nx
      ! dhdk = nx*ny
      ! dh = dhdi*di + dhdj*dj + dhdk*dk
      ! dh = 1*di + nx*dj + nx*ny*dk QED

      implicit none
      integer, intent(in):: ncells(3), dcells(3)

      dhash = ncells(1)*ncells(2)*dcells(3) + ncells(1)*dcells(2) + dcells(1)

   end function dhash

   function globalhash2imagehash(hash, numImages, ncells) result(imagehash)

      implicit none
      integer, intent(in):: hash, ncells, numImages
      integer:: imagehash(2), ncells_image_avg

      ncells_image_avg = ceiling(real(ncells, f)/numImages)
      imagehash(2) = ceiling(real(hash, f)/ncells_image_avg)
      imagehash(1) = hash - (imagehash(2)-1)*ncells_image_avg

end function globalhash2imagehash

end module flink_list_m
