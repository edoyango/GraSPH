module flink_list_m

   use datatypes, only: particles, interactions
   use param, only: dim, hsml, f, halotype

   !use error_msg_m, only: error_msg
   use kernel_m, only: kernel

   private
   public:: flink_list

contains

   !==============================================================================================================================
   subroutine flink_list(maxinter, scale_k, ntotal_loc, nhalo_loc, nvirt_loc, niac, parts, pairs, nexti)
      ! save as above, but for 3D

      implicit none
      integer, intent(in):: maxinter, ntotal_loc, nhalo_loc, nvirt_loc
      real(f), intent(in):: scale_k
      type(particles), intent(in):: parts(:)
      integer, intent(out):: niac, nexti(:)
      type(interactions), intent(out):: pairs(:)
      integer, parameter:: maxpcell = 125
      integer:: ngridx(3), jth, i, j, k, d, icell, jcell, kcell, xi, yi, zi, ierr = 0
      real(f):: mingridx(3), maxgridx(3), dcell
      integer, allocatable:: pincell(:, :, :), gridind(:, :), cells(:, :, :, :)
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
      do i = 2, ntotal_loc + nhalo_loc + nvirt_loc
         do d = 1, dim
            mingridx(d) = MIN(mingridx(d), parts(i)%x(d))
            maxgridx(d) = MAX(maxgridx(d), parts(i)%x(d))
         end do
      end do

      !Determining number of grid cells in each direction
      dcell = scale_k*hsml
      maxgridx(:) = maxgridx(:) + 2._f*dcell
      mingridx(:) = mingridx(:) - 2._f*dcell
      ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
      maxgridx(:) = mingridx(:) + ngridx(:)*dcell

      allocate (pincell(ngridx(1), ngridx(2), ngridx(3)), &
                gridind(dim, ntotal_loc + nhalo_loc + nvirt_loc), &
                cells(maxpcell, ngridx(1), ngridx(2), ngridx(3)))

      pincell(:, :, :) = 0

      do i = 1, ntotal_loc + nhalo_loc + nvirt_loc
         gridind(:, i) = int((parts(i)%x(:) - mingridx(:))/dcell) + 1
         icell = gridind(1, i)
         jcell = gridind(2, i)
         kcell = gridind(3, i)
         pincell(icell, jcell, kcell) = pincell(icell, jcell, kcell) + 1
         cells(pincell(icell, jcell, kcell), icell, jcell, kcell) = i
      end do

      niac = 0
      ierr = 0
      nexti(1) = 1
      do i = 1, ntotal_loc + nhalo_loc + nvirt_loc

         ! Retrieving particle i's grid cell indices
         icell = gridind(1, i)
         jcell = gridind(2, i)
         kcell = gridind(3, i)

         ! finding pairs within cell icell,jcell
         do j = 1, pincell(icell, jcell, kcell)
            jth = cells(j, icell, jcell, kcell)
            if (jth > i) then
               call check_if_interact(maxinter, scale_k, i, jth, parts(i), parts(jth), niac, pairs, ierr)
            end if
         end do

         ! finding pairs within cells adjacent to i's cell
         do k = 1, 13
            xi = icell + sweep(1, k)
            yi = jcell + sweep(2, k)
            zi = kcell + sweep(3, k)
            do j = 1, pincell(xi, yi, zi)
               jth = cells(j, xi, yi, zi)
               call check_if_interact(maxinter, scale_k, i, jth, parts(i), parts(jth), niac, pairs, ierr)
            end do
         end do

         nexti(i + 1) = niac + 1

      end do

      if (ierr == 1) then
         print *, ' >>> Error <<< : Too many interactions'
         stop
      end if

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
               pairs(niac)%dx = dxiac
            else
               ierr = 1
            end if
         end if
      end if

   end subroutine check_if_interact

end module flink_list_m
