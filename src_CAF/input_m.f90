module input_m

   use datatypes, only: particles, interactions
   use param, only: dim, irho, dxo, f, hsml, mp, np, op, pp, qp, rp, nlayer, mass, halotype

   private
   real(f), parameter:: vxmin = 0._f, vymin = 0._f, vzmin = 0._f, &
                        vxmax = vxmin + pp*dxo, vymax = vymin + qp*dxo, vzmax = vzmin + rp*dxo
   real(f), parameter:: rxmin = 0._f, rymin = 0._f, rzmin = 0._f, &
                        rxmax = rxmin + mp*dxo, rymax = rymin + np*dxo, rzmax = rzmin + op*dxo

   public:: return_ntotal, return_nvirt, allocatePersistentArrays, generate_real_part, generate_virt_part, &
            update_virt_part

contains

   !====================================================================================================================
   function return_ntotal() result(ntotal)

      implicit none
      integer:: ntotal

      ntotal = mp*np*op

   end function

   !====================================================================================================================
   function return_nvirt() result(nvirt)

      implicit none
      integer:: nvirt

      nvirt = (rp + 2*nlayer)*(2*nlayer + pp)*(2*nlayer + qp) - pp*qp*rp

   end function

   !====================================================================================================================
   subroutine allocatePersistentArrays(ntotal, nvirt, parts, pairs, nexti, maxnloc, maxinter)

      implicit none
      integer, intent(in):: ntotal, nvirt
      integer, intent(out), allocatable:: nexti(:)
      type(particles), intent(inout), allocatable, codimension[:]:: parts(:)
      type(interactions), intent(out), allocatable:: pairs(:)
      integer, intent(out):: maxnloc, maxinter

      maxnloc = ntotal + nvirt
      maxinter = 262*maxnloc

      allocate (parts(maxnloc) [*], pairs(maxinter), nexti(maxnloc + 1))

   end subroutine allocatePersistentArrays

   !====================================================================================================================
   subroutine generate_real_part(thisImage, numImages, ntotal, ntotal_loc, parts)
      ! Generates initial physical particle configuration.
      ! 2 cases: return only number of particles retrieved, or generating the particles

      implicit none
      integer, intent(in):: thisImage, numImages, ntotal
      integer, intent(out):: ntotal_loc
      type(particles), intent(out):: parts(:)
      integer:: i, j, k, n, n_loc, n_loc_i, n_start, n_done

      ! how many particles to generate per process
      n_loc_i = ceiling(dble(ntotal)/numImages)
      if (thisImage .eq. numImages) then
         n_loc = ntotal - (numImages - 1)*n_loc_i
      else
         n_loc = n_loc_i
      end if
      n_start = (thisImage - 1)*n_loc_i + 1
      n_done = n_start + n_loc_i - 1

      ! stopping program if array bounds are exceeded
      !if ((procid .eq. 0) .and. (n_loc .gt. maxnloc)) call error_msg(1, 1)

      ! intitial setup
      n = 0
      ntotal_loc = 0
      do i = 1, mp
         do j = 1, np
            do k = 1, op
               n = n + 1 ! tracking total number of particles generated
               ! Only generating particles assigned to process
               if ((n .ge. n_start) .and. (n .le. n_done)) then
                  ntotal_loc = ntotal_loc + 1
                  parts(ntotal_loc)%indglob = n
                  parts(ntotal_loc)%indloc = ntotal_loc
                  parts(ntotal_loc)%x(1) = rxmin + (i - 0.5_f)*dxo
                  parts(ntotal_loc)%x(2) = rymin + (j - 0.5_f)*dxo
                  parts(ntotal_loc)%x(3) = rzmin + (k - 0.5_f)*dxo
                  parts(ntotal_loc)%vx(:) = 0._f
                  parts(ntotal_loc)%itype = 1
                  parts(ntotal_loc)%rho = irho
                  parts(ntotal_loc)%p = 0._f
               end if
            end do
         end do
      end do

   end subroutine generate_real_part

   !====================================================================================================================
   subroutine generate_virt_part(thisImage, numImages, ntotal, ntotal_loc, nvirt, nvirt_loc, parts)

      ! Generates the virtual particle configuration. Can change over time or remain static
      ! 2 cases: return only number of particles retrieved, or generating the particles

      implicit none
      integer, intent(in):: thisImage, numImages, ntotal_loc, ntotal, nvirt
      type(particles), intent(inout):: parts(:)
      integer, intent(out):: nvirt_loc
      integer:: i, j, k, n, n_loc, n_loc_i, n_start, n_done

      ! how many particles to generate per process
      n_loc_i = ceiling(dble(nvirt)/numImages)
      if (thisImage .eq. numImages) then
         n_loc = nvirt - (numImages - 1)*n_loc_i
      else
         n_loc = n_loc_i
      end if
      n_start = (thisImage - 1)*n_loc_i + 1
      n_done = n_start + n_loc_i - 1
      write(*,*) thisImage, n_start, n_done
      nvirt_loc = 0
      n = 0 ! counter used to track particle indices

      !---Virtual particle on the bottom face
      do i = 1 - nlayer, pp + nlayer
         do j = 1 - nlayer, qp + nlayer
            do k = 1 - nlayer, rp + nlayer
               if (i < 1 .or. i > pp .or. j < 1 .or. j > qp .or. k < 1 .or. k > rp) then
                  n = n + 1
                  ! Only generating particles assigned to process
                  if ((n .ge. n_start) .and. (n .le. n_done)) then
                     nvirt_loc = nvirt_loc + 1
                     parts(ntotal_loc + nvirt_loc)%indglob = n + ntotal
                     parts(ntotal_loc + nvirt_loc)%indloc = ntotal_loc + nvirt_loc
                     if (k < 1) parts(ntotal_loc + nvirt_loc)%itype = -1
                     if (k > rp) parts(ntotal_loc + nvirt_loc)%itype = -2
                     if (j < 1 .or. j > qp) parts(ntotal_loc + nvirt_loc)%itype = -4
                     if (i < 1 .or. i > pp) parts(ntotal_loc + nvirt_loc)%itype = -3
                     if ( (i < 1 .and. j < 1) .or. (i < 1 .and. j > qp) .or. (i > pp .and. j < 1) .or. &
                        (i > pp .and. j > qp)) parts(ntotal_loc + nvirt_loc)%itype = -5
                     parts(ntotal_loc + nvirt_loc)%x(1) = vxmin + (i - 0.5_f)*dxo
                     parts(ntotal_loc + nvirt_loc)%x(2) = vymin + (j - 0.5_f)*dxo
                     parts(ntotal_loc + nvirt_loc)%x(3) = vzmin + (k - 0.5_f)*dxo
                     parts(ntotal_loc + nvirt_loc)%vx(:) = 0._f
                     parts(ntotal_loc + nvirt_loc)%rho = irho
                  end if
               end if
            end do
         end do
      end do

   end subroutine generate_virt_part

   !==============================================================================================================================
   subroutine update_virt_part(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, nexti, vw)

      implicit none
      integer, intent(in):: ntotal_loc, nhalo_loc, nvirt_loc, niac, nexti(:)
      type(interactions), intent(in):: pairs(:)
      type(particles), intent(inout):: parts(:)
      real(f), intent(inout):: vw(:)
      integer:: i, j, k, tmptype
      real(f):: tmp

      vw(:) = 0._f

      do i = 1, ntotal_loc + nvirt_loc + nhalo_loc
         if (parts(i)%itype<0) then
            parts(i)%rho = 0._f
            parts(i)%vx(:) = 0._f
         end if
      end do

      do i = 1, ntotal_loc + nhalo_loc + nvirt_loc
         do k = nexti(i), nexti(i + 1) - 1
            j = pairs(k)%j
            if ( parts(i)%itype < 0 .and. parts(j)%itype > 0) then
               if (parts(i)%itype < -halotype) then
                  tmptype = parts(i)%itype + halotype
               else
                  tmptype = parts(i)%itype
               end if
               tmp = mass*pairs(k)%w/parts(j)%rho
               vw(i) = vw(i) + tmp
               parts(i)%rho = parts(i)%rho + mass*pairs(k)%w
               select case (tmptype)
               case default
                  parts(i)%vx(:) = parts(i)%vx(:) - parts(j)%vx(:)*tmp
               case (-2)
                  parts(i)%vx(1) = parts(i)%vx(1) + parts(j)%vx(1)*tmp
                  parts(i)%vx(2) = parts(i)%vx(2) + parts(j)%vx(2)*tmp
                  parts(i)%vx(3) = parts(i)%vx(3) - parts(j)%vx(3)*tmp
               case (-3)
                  parts(i)%vx(1) = parts(i)%vx(1) - parts(j)%vx(1)*tmp
                  parts(i)%vx(2) = parts(i)%vx(2) + parts(j)%vx(2)*tmp
                  parts(i)%vx(3) = parts(i)%vx(3) + parts(j)%vx(3)*tmp
               case (-4)
                  parts(i)%vx(1) = parts(i)%vx(1) + parts(j)%vx(1)*tmp
                  parts(i)%vx(2) = parts(i)%vx(2) - parts(j)%vx(2)*tmp
                  parts(i)%vx(3) = parts(i)%vx(3) + parts(j)%vx(3)*tmp
               case (-5)
                  parts(i)%vx(1) = parts(i)%vx(1) - parts(j)%vx(1)*tmp
                  parts(i)%vx(2) = parts(i)%vx(2) - parts(j)%vx(2)*tmp
                  parts(i)%vx(3) = parts(i)%vx(3) + parts(j)%vx(3)*tmp
               end select
            else if (parts(j)%itype < 0 .and. parts(i)%itype > 0) then
               if (parts(j)%itype < -halotype) then
                  tmptype = parts(j)%itype + halotype
               else
                  tmptype = parts(j)%itype
               end if
               tmp = mass*pairs(k)%w/parts(i)%rho
               vw(j) = vw(j) + tmp
               parts(j)%rho = parts(j)%rho + mass*pairs(k)%w
               select case (tmptype)
               case default
                  parts(j)%vx(:) = parts(j)%vx(:) - parts(i)%vx(:)*tmp
               case (-2)
                  parts(j)%vx(1) = parts(j)%vx(1) + parts(i)%vx(1)*tmp
                  parts(j)%vx(2) = parts(j)%vx(2) + parts(i)%vx(2)*tmp
                  parts(j)%vx(3) = parts(j)%vx(3) - parts(i)%vx(3)*tmp
               case (-3)
                  parts(j)%vx(1) = parts(j)%vx(1) - parts(i)%vx(1)*tmp
                  parts(j)%vx(2) = parts(j)%vx(2) + parts(i)%vx(2)*tmp
                  parts(j)%vx(3) = parts(j)%vx(3) + parts(i)%vx(3)*tmp
               case (-4)
                  parts(j)%vx(1) = parts(j)%vx(1) + parts(i)%vx(1)*tmp
                  parts(j)%vx(2) = parts(j)%vx(2) - parts(i)%vx(2)*tmp
                  parts(j)%vx(3) = parts(j)%vx(3) + parts(i)%vx(3)*tmp
               case (-5)
                  parts(j)%vx(1) = parts(j)%vx(1) - parts(i)%vx(1)*tmp
                  parts(j)%vx(2) = parts(j)%vx(2) - parts(i)%vx(2)*tmp
                  parts(j)%vx(3) = parts(j)%vx(3) + parts(i)%vx(3)*tmp
               end select
            end if
         end do
      end do

      do i = 1, ntotal_loc+nvirt_loc+nhalo_loc
         if (parts(i)%itype<0) then
            if (vw(i) > 0._f) then
               parts(i)%rho = parts(i)%rho/vw(i)
               parts(i)%vx(:) = parts(i)%vx(:)/vw(i)
            else
               parts(i)%rho = irho
               parts(i)%vx(:) = 0._f
            end if
         end if
      end do

   end subroutine update_virt_part

end module input_m
