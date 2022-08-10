module flink_list_m

   use datatypes, only: particles, interactions
   use param, only: dim, f, hsml

   use kernel_m, only: kernel

   public:: flink_list
   private:: check_if_interact

contains

   !==============================================================================================================================
   subroutine flink_list(scale_k, ntotal, nvirt, nghos, parts, niac, pairs, nexti)
      ! save as above, but for 3D

      implicit none
      real(f), intent(in):: scale_k
      integer, intent(in):: ntotal, nvirt, nghos
      type(particles), intent(in):: parts(:)
      integer, intent(out):: niac, nexti(:)
      type(interactions), intent(out):: pairs(:)
      integer, parameter:: maxpcell = 125
      integer:: ngridx(dim), i, j, k, d, icell, jcell, kcell, xi, yi, zi, jth, maxinter, ierr = 0
      real(f):: mingridx(dim), maxgridx(dim), dcell
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
      maxinter = size(pairs)

      !Determining bounding box extents
      mingridx(:) = parts(1)%x(:)
      maxgridx(:) = parts(1)%x(:)
      do i = 2, ntotal + nvirt + nghos
         do d = 1, dim
            mingridx(d) = min(mingridx(d), parts(i)%x(d))
            maxgridx(d) = max(maxgridx(d), parts(i)%x(d))
         end do
      end do

      !Determining number of grid cells in each direction
      dcell = scale_k*hsml
      maxgridx(:) = maxgridx(:) + 2._f*dcell
      mingridx(:) = mingridx(:) - 2._f*dcell
      ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
      maxgridx(:) = mingridx(:) + ngridx(:)*dcell

      allocate (pincell(ngridx(1), ngridx(2), ngridx(3)), &
                gridind(dim, ntotal + nvirt + nghos), &
                cells(maxpcell, ngridx(1), ngridx(2), ngridx(3)))

      pincell(:, :, :) = 0

      do i = 1, ntotal + nvirt + nghos
         gridind(:, i) = int((parts(i)%x(:) - mingridx(:))/dcell) + 1
         icell = gridind(1, i)
         jcell = gridind(2, i)
         kcell = gridind(3, i)
         pincell(icell, jcell, kcell) = pincell(icell, jcell, kcell) + 1
         cells(pincell(icell, jcell, kcell), icell, jcell, kcell) = i
      end do

      niac = 0
      nexti(1) = 1
      do i = 1, ntotal + nvirt + nghos
         icell = gridind(1, i)
         jcell = gridind(2, i)
         kcell = gridind(3, i)
         if (pincell(icell, jcell, kcell) > 1) then

            ! finding pairs within cell icell,jcell
            do j = 1, pincell(icell, jcell, kcell)
               jth = cells(j, icell, jcell, kcell)
               if (i < cells(j, icell, jcell, kcell)) call check_if_interact(maxinter, scale_k, parts(i), parts(jth), niac, pairs, &
                                                                             ierr)
            end do
         end if

         ! finding pairs within cells adjacent to i's cell
         do k = 1, 13
            xi = icell + sweep(1, k)
            yi = jcell + sweep(2, k)
            zi = kcell + sweep(3, k)

            do j = 1, pincell(xi, yi, zi)
               jth = cells(j, xi, yi, zi)
               call check_if_interact(maxinter, scale_k, parts(i), parts(jth), niac, pairs, ierr)
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
   pure subroutine check_if_interact(maxinter, scale_k, p_i, p_j, niac, pairs, ierr)
      ! subroutine to chekc if two particles are interacting and consequently adding to pair list

      implicit none
      integer, intent(in):: maxinter
      real(f), intent(in):: scale_k
      type(particles), intent(in):: p_i, p_j
      integer, intent(inout):: niac, ierr
      type(interactions), intent(inout):: pairs(:)
      real(f):: dxiac(dim), r

      ! only consider interactions when real-real are involved
      if (p_i%itype .eq. 1 .or. p_j%itype .eq. 1) then
         dxiac(:) = p_i%x(:) - p_j%x(:)
         r = SQRT(SUM(dxiac*dxiac))
         if (r < hsml*scale_k) then
            niac = niac + 1
            if (niac < maxinter) then
!~                pairs(niac)%i = p_i%ind
               pairs(niac)%j = p_j%indloc
               call kernel(r, dxiac, hsml, pairs(niac)%w, pairs(niac)%dwdx(:))
            else
               ierr = 1
            end if
         end if
      end if

   end subroutine check_if_interact

end module flink_list_m
