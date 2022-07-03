module input_m

   use datatypes, only: particles, interactions
   use mpi_f08
   use param, only: dim, irho, dxo, f, hsml, mp, np, op, pp, qp, rp, nlayer

   private
   real(f), parameter:: vxmin = 0._f, vymin = 0._f, vzmin = 0._f, &
                        vxmax = vxmin + pp*dxo, vymax = vymin + qp*dxo, vzmax = vzmin + rp*dxo
   real(f), parameter:: rxmin = 0._f, rymin = 0._f, rzmin = 0._f, &
                        rxmax = rxmin + mp*dxo, rymax = rymin + np*dxo, rzmax = rzmin + op*dxo

   public:: return_ntotal, return_nvirt, allocatePersistentArrays, generate_real_part, generate_virt_part, generate_ghost_part, &
            update_ghost_part, virt_mirror

contains

   !==============================================================================================================================
   function return_ntotal() result(ntotal)

      implicit none
      integer:: ntotal

      ntotal = mp*np*op

   end function

   !==============================================================================================================================
   function return_nvirt() result(nvirt)

      implicit none
      integer:: nvirt

      nvirt = nlayer*(2*nlayer + pp)*(2*nlayer + qp)

   end function

   !==============================================================================================================================
   subroutine allocatePersistentArrays(ntotal, nvirt, parts, pairs, nexti, gind, maxnloc, maxinter)

      implicit none
      integer, intent(in):: ntotal, nvirt
      integer, intent(out), allocatable:: gind(:), nexti(:)
      type(particles), intent(out), allocatable:: parts(:)
      type(interactions), intent(out), allocatable:: pairs(:)
      integer, intent(out):: maxnloc, maxinter

      maxnloc = 2*ntotal + nvirt
      maxinter = 262*maxnloc

      allocate (parts(maxnloc), pairs(maxinter), nexti(maxnloc + 1), gind(ntotal))

   end subroutine allocatePersistentArrays

   !==============================================================================================================================
   subroutine generate_real_part(procid, numprocs, ntotal, ntotal_loc, parts)
      ! Generates initial physical particle configuration.
      ! 2 cases: return only number of particles retrieved, or generating the particles

      implicit none
      integer, intent(in):: procid, numprocs, ntotal
      integer, intent(out):: ntotal_loc
      type(particles), intent(out):: parts(:)
      integer:: i, j, k, n, n_loc, n_loc_i, n_start, n_done

      ! how many particles to generate per process
      n_loc_i = ceiling(dble(ntotal)/numprocs)
      if (procid .eq. numprocs - 1) then
         n_loc = ntotal - (numprocs - 1)*n_loc_i
      else
         n_loc = n_loc_i
      end if
      n_start = procid*n_loc_i + 1
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

   !==============================================================================================================================
   subroutine generate_virt_part(procid, bounds_loc, scale_k, ntotal, ntotal_loc, nhalo_loc, nvirt_loc, parts)
      ! Generates the virtual particle configuration. Can change over time or remain static
      ! 2 cases: return only number of particles retrieved, or generating the particles

      implicit none
      integer, intent(in):: procid, ntotal_loc, nhalo_loc, ntotal
      real(f), intent(in):: bounds_loc(2*dim), scale_k
      type(particles), intent(inout):: parts(:)
      integer, intent(out):: nvirt_loc
      integer:: i, j, k, n
      real(f):: xi(dim), xmin_loc(dim), xmax_loc(dim)

      xmin_loc(:) = bounds_loc(1:dim) - scale_k*hsml
      xmax_loc(:) = bounds_loc(dim + 1:2*dim) + scale_k*hsml

      nvirt_loc = 0
      n = ntotal ! counter used to track particle indices

      !---Virtual particle on the bottom face
      do i = 1 - nlayer, pp + nlayer
         do j = 1 - nlayer, qp + nlayer
            do k = 1, nlayer
               n = n + 1
               xi(1) = vxmin + (i - 0.5_f)*dxo
               xi(2) = vymin + (j - 0.5_f)*dxo
               xi(3) = vzmin - (k - 0.5_f)*dxo
               if (all(xi(:) .ge. xmin_loc(:) .and. xi(:) .le. xmax_loc(:))) then
                  nvirt_loc = nvirt_loc + 1
                  parts(ntotal_loc + nhalo_loc + nvirt_loc)%indglob = n
                  parts(ntotal_loc + nhalo_loc + nvirt_loc)%indloc = ntotal_loc + nhalo_loc + nvirt_loc
                  parts(ntotal_loc + nhalo_loc + nvirt_loc)%itype = -1
                  parts(ntotal_loc + nhalo_loc + nvirt_loc)%x(:) = xi(:)
                  parts(ntotal_loc + nhalo_loc + nvirt_loc)%vx(:) = 0._f
                  parts(ntotal_loc + nhalo_loc + nvirt_loc)%rho = irho
               end if
            end do
         end do
      end do

   end subroutine generate_virt_part

   !==============================================================================================================================
   pure subroutine generate_ghost_part(scale_k, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, parts, gind)

      implicit none
      integer, intent(in):: ntotal_loc, nhalo_loc, nvirt_loc
      real(f), intent(in):: scale_k
      type(particles), intent(inout):: parts(:)
      integer, intent(out):: nghos_loc, gind(:)
      integer:: i, ig

      nghos_loc = 0

      do i = 1, ntotal_loc + nhalo_loc
         if (abs(parts(i)%x(1) - vxmin) < scale_k*hsml .and. parts(i)%x(1) > vxmin) then
            nghos_loc = nghos_loc + 1
            ig = ntotal_loc + nhalo_loc + nvirt_loc + nghos_loc
            gind(nghos_loc) = i
            parts(ig) = parts(i)
            parts(ig)%indloc = ig
            parts(ig)%itype = 99
            parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmin
            parts(ig)%vx(1) = -parts(ig)%vx(1)
         end if
         if (abs(parts(i)%x(1) - vxmax) < scale_k*hsml .and. parts(i)%x(1) < vxmax) then
            nghos_loc = nghos_loc + 1
            ig = ntotal_loc + nhalo_loc + nvirt_loc + nghos_loc
            gind(nghos_loc) = i
            parts(ig) = parts(i)
            parts(ig)%indloc = ig
            parts(ig)%itype = 99
            parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmax
            parts(ig)%vx(1) = -parts(ig)%vx(1)
         end if
         if (abs(parts(i)%x(2) - vymin) < scale_k*hsml .and. parts(i)%x(2) > vymin) then
            nghos_loc = nghos_loc + 1
            ig = ntotal_loc + nhalo_loc + nvirt_loc + nghos_loc
            gind(nghos_loc) = i
            parts(ig) = parts(i)
            parts(ig)%indloc = ig
            parts(ig)%itype = 98
            parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymin
            parts(ig)%vx(2) = -parts(ig)%vx(2)
         end if
         if (abs(parts(i)%x(2) - vymax) < scale_k*hsml .and. parts(i)%x(2) < vymax) then
            nghos_loc = nghos_loc + 1
            ig = ntotal_loc + nhalo_loc + nvirt_loc + nghos_loc
            gind(nghos_loc) = i
            parts(ig) = parts(i)
            parts(ig)%indloc = ig
            parts(ig)%itype = 98
            parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymax
            parts(ig)%vx(2) = -parts(ig)%vx(2)
         end if
         if ((parts(i)%x(1) - vxmin)**2 + (parts(i)%x(2) - vymin)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) > vxmin) then
            nghos_loc = nghos_loc + 1
            ig = ntotal_loc + nhalo_loc + nvirt_loc + nghos_loc
            gind(nghos_loc) = i
            parts(ig) = parts(i)
            parts(ig)%indloc = ig
            parts(ig)%itype = 97
            parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmin
            parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymin
            parts(ig)%vx(1) = -parts(ig)%vx(1)
            parts(ig)%vx(2) = -parts(ig)%vx(2)
         end if
         if ((parts(i)%x(1) - vxmin)**2 + (parts(i)%x(2) - vymax)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) > vxmin) then
            nghos_loc = nghos_loc + 1
            ig = ntotal_loc + nhalo_loc + nvirt_loc + nghos_loc
            gind(nghos_loc) = i
            parts(ig) = parts(i)
            parts(ig)%indloc = ig
            parts(ig)%itype = 97
            parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmin
            parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymax
            parts(ig)%vx(1) = -parts(ig)%vx(1)
            parts(ig)%vx(2) = -parts(ig)%vx(2)
         end if
         if ((parts(i)%x(1) - vxmax)**2 + (parts(i)%x(2) - vymax)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) < vxmax) then
            nghos_loc = nghos_loc + 1
            ig = ntotal_loc + nhalo_loc + nvirt_loc + nghos_loc
            gind(nghos_loc) = i
            parts(ig) = parts(i)
            parts(ig)%indloc = ig
            parts(ig)%itype = 97
            parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmax
            parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymax
            parts(ig)%vx(1) = -parts(ig)%vx(1)
            parts(ig)%vx(2) = -parts(ig)%vx(2)
         end if
         if ((parts(i)%x(1) - vxmax)**2 + (parts(i)%x(2) - vymin)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) < vxmax) then
            nghos_loc = nghos_loc + 1
            ig = ntotal_loc + nhalo_loc + nvirt_loc + nghos_loc
            gind(nghos_loc) = i
            parts(ig) = parts(i)
            parts(ig)%indloc = ig
            parts(ig)%itype = 97
            parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmax
            parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymin
            parts(ig)%vx(1) = -parts(ig)%vx(1)
            parts(ig)%vx(2) = -parts(ig)%vx(2)
         end if
      end do

   end subroutine generate_ghost_part

   !==============================================================================================================================
   pure subroutine update_ghost_part(gind, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, parts)

      implicit none
      integer, intent(in):: ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, gind(:)
      type(particles), intent(inout):: parts(:)
      integer:: i, ig, ir

      do i = 1, nghos_loc
         ig = ntotal_loc + nhalo_loc + nvirt_loc + i
         ir = gind(i)
         select case (parts(ig)%itype)
         case (99)
            parts(ig)%rho = parts(ir)%rho
            parts(ig)%p = parts(ir)%p
            parts(ig)%vx(1) = -parts(ir)%vx(1)
            parts(ig)%vx(2) = parts(ir)%vx(2)
            parts(ig)%vx(3) = parts(ir)%vx(3)
         case (98)
            parts(ig)%rho = parts(ir)%rho
            parts(ig)%p = parts(ir)%p
            parts(ig)%vx(1) = parts(ir)%vx(1)
            parts(ig)%vx(2) = -parts(ir)%vx(2)
            parts(ig)%vx(3) = parts(ir)%vx(3)
         case (97)
            parts(ig)%rho = parts(ir)%rho
            parts(ig)%p = parts(ir)%p
            parts(ig)%vx(1) = -parts(ir)%vx(1)
            parts(ig)%vx(2) = -parts(ir)%vx(2)
            parts(ig)%vx(3) = parts(ir)%vx(3)
         end select
      end do

   end subroutine update_ghost_part

   !==============================================================================================================================
   pure subroutine virt_mirror(pr, pv)

      implicit none
      type(particles), intent(in):: pr
      type(particles), intent(inout):: pv
      real(f):: da, db, beta
      real(f), parameter:: beta_max = 5._f

      da = ABS(pr%x(3) - vzmin)
      db = ABS(pv%x(3) - vzmin)

      beta = MIN(1._f + db/da, beta_max)
      if (ISNAN(beta)) beta = beta_max

      pv%rho = pr%rho
      pv%p = pr%p
      pv%vx(:) = (1._f - beta)*pr%vx(:)

   end subroutine virt_mirror

end module input_m
