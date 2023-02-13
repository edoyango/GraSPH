module input_m

   use datatypes, only: particles, interactions
   use param, only: dim, irho, dxo, f, hsml, mp, np, op, pp, qp, rp, nlayer, mass

   private
   real(f), parameter:: vxmin = 0._f, vymin = 0._f, vzmin = 0._f, &
                        vxmax = vxmin + pp*dxo, vymax = vymin + qp*dxo, vzmax = vzmin + rp*dxo
   real(f), parameter:: rxmin = 0._f, rymin = 0._f, rzmin = 0._f, &
                        rxmax = rxmin + mp*dxo, rymax = rymin + np*dxo, rzmax = rzmin + op*dxo

   public:: return_ntotal, return_nvirt, allocatePersistentArrays, generate_real_part

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

      allocate (parts(maxnloc)[*], pairs(maxinter), nexti(maxnloc + 1))

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
      n_start = (thisImage-1)*n_loc_i + 1
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

end module input_m
