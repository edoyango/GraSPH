module flink_list_m

   use datatypes, only: particles, interactions
   use param, only: dim, hsml, f, halotype, scale_k

   private
   public:: flink_list

contains

   !==============================================================================================================================
   subroutine flink_list(maxinter, niac, parts, pairs, nexti)
      ! save as above, but for 3D

      implicit none
      integer, intent(in):: maxinter
      type(particles), intent(in):: parts
      integer, intent(out):: niac, nexti(:)
      type(interactions), intent(out):: pairs(:)
      type(particles):: p_i
      integer, parameter:: maxpcell = 125
      integer:: ngridx(3), jth, i, j, k, d, icell, jcell, kcell, xi, yi, zi, itypetmp
      real(f):: mingridx(3), maxgridx(3), dcell, xtmp(dim)
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
      mingridx(:) = parts%x(:, 1)
      maxgridx(:) = parts%x(:, 1)
      do i = 2, parts%nsum()
         do d = 1, dim
            mingridx(d) = MIN(mingridx(d), parts%x(d, i))
            maxgridx(d) = MAX(maxgridx(d), parts%x(d, i))
         end do
      end do

      !Determining number of grid cells in each direction
      dcell = scale_k*hsml
      maxgridx(:) = maxgridx(:) + 2._f*dcell
      mingridx(:) = mingridx(:) - 2._f*dcell
      ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
      maxgridx(:) = mingridx(:) + ngridx(:)*dcell

      allocate (pincell(ngridx(1), ngridx(2), ngridx(3)), &
                gridind(dim, parts%nsum()), &
                cells(maxpcell, ngridx(1), ngridx(2), ngridx(3)))

      pincell(:, :, :) = 0

      do i = 1, parts%nsum()
         gridind(:, i) = int((parts%x(:, i) - mingridx(:))/dcell) + 1
         icell = gridind(1, i)
         jcell = gridind(2, i)
         kcell = gridind(3, i)
         pincell(icell, jcell, kcell) = pincell(icell, jcell, kcell) + 1
         cells(pincell(icell, jcell, kcell), icell, jcell, kcell) = i
      end do

      niac = 0
      nexti(1) = 1
      do i = 1, parts%nsum()

         ! Retrieving particle i's grid cell indices
         icell = gridind(1, i)
         jcell = gridind(2, i)
         kcell = gridind(3, i)

         xtmp(:) = parts%x(:, i)
         itypetmp = parts%itype(i)

         ! finding pairs within cell icell,jcell
         do j = 1, pincell(icell, jcell, kcell)
            jth = cells(j, icell, jcell, kcell)
            if (jth > i) then
               call check_if_interact(maxinter, i, jth, itypetmp, xtmp, parts%itype(i), parts%x(:, jth), niac, pairs)
            end if
         end do

         ! finding pairs within cells adjacent to i's cell
         do k = 1, 13
            xi = icell + sweep(1, k)
            yi = jcell + sweep(2, k)
            zi = kcell + sweep(3, k)
            do j = 1, pincell(xi, yi, zi)
               jth = cells(j, xi, yi, zi)
               call check_if_interact(maxinter, i, jth, itypetmp, xtmp, parts%itype(jth), parts%x(:, jth), niac, pairs)
            end do
         end do

         nexti(i + 1) = niac + 1

      end do

   end subroutine flink_list

   !==============================================================================================================================
   pure subroutine check_if_interact(maxinter, i, j, itype_i, x_i, itype_j, x_j, niac, pairs)
      ! subroutine to chekc if two particles are interacting and consequently adding to pair list

      implicit none
      integer, intent(in):: maxinter, i, j, itype_i, itype_j
      real(f), intent(in):: x_i(dim), x_j(dim)
      integer, intent(inout):: niac
      type(interactions), intent(inout):: pairs(:)
      real(f):: dxiac(dim), r

      ! only consider interactions between:
      ! real-real
      ! real-virtual
      ! real-halo
      ! halo-virtual - for deterimining boundary
      ! if (.not.((p_i%itype<0 .and. p_j%itype<0) .or. (p_i%itype>halotype .and. p_j%itype>halotype))) then
      if (itype_i>0 .or. itype_j>0) then
         dxiac(:) = x_i(:) - x_j(:)
         ! r = SQRT(SUM(dxiac*dxiac))
         if (SUM(dxiac*dxiac) < hsml*hsml*scale_k*scale_k) then
            niac = niac + 1
#ifdef DEBUG
            if (niac < maxinter) then
#endif
               pairs(niac)%j = j
               pairs(niac)%dx = dxiac
#ifdef DEBUG
            else
               error stop "Too many interactions!"
            end if
#endif
         end if
      end if

   end subroutine check_if_interact

end module flink_list_m
