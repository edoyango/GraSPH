module input_m

   use datatypes, only: particles, interactions
   use param, only: dim, f, dxo, mp, np, op, pp, qp, rp, nlayer, irho, hsml, mass, rh0, gamma, c

   real(f), parameter:: vxmin = 0._f, vymin = 0._f, vzmin = 0._f, &
                        vxmax = vxmin + pp*dxo, vymax = vymin + qp*dxo, vzmax = vzmin + rp*dxo
   real(f), parameter:: rxmin = 0._f, rymin = 0._f, rzmin = 0._f, &
                        rxmax = rxmin + mp*dxo, rymax = rymin + np*dxo, rzmax = rzmin + op*dxo

   public:: generate_real_part, generate_virt_part, return_ntotal, return_nvirt
   private:: vxmin, vymin, vzmin, vxmax, vymax, vzmax, rxmin, rymin, rzmin, rxmax, rymax, rzmax

contains

   !==============================================================================================================================
   pure function return_ntotal() result(ntotal)

      implicit none
      integer:: ntotal

      ntotal = mp*np*op

   end function return_ntotal

   !==============================================================================================================================
   pure function return_nvirt() result(nvirt)

      implicit none
      integer:: nvirt, i, j, k

      nvirt = 0
      do i = 1 - nlayer, pp + nlayer
         do j = 1 - nlayer, qp + nlayer
            do k = 1 - nlayer, rp + nlayer
               if (i < 1 .or. i > pp .or. j < 1 .or. j > qp .or. k < 1 .or. k > rp) then
                  nvirt = nvirt + 1
               end if
            end do
         end do
      end do

   end function return_nvirt

   !==============================================================================================================================
   subroutine allocatePersistentArrays(ntotal, nvirt, parts, pairs, nexti, gind)

      implicit none
      integer, intent(in):: ntotal, nvirt
      type(particles), allocatable, intent(out):: parts(:)
      type(interactions), allocatable, intent(out):: pairs(:)
      integer, allocatable, intent(out):: nexti(:), gind(:)
      integer:: maxn, maxinter

      maxn = 2*ntotal + nvirt
      maxinter = 262*maxn

      allocate (parts(maxn), pairs(maxinter), nexti(maxn + 1), gind(ntotal))

   end subroutine allocatePersistentArrays

   !==============================================================================================================================
   pure subroutine generate_real_part(ntotal, parts)

      implicit none
      integer, intent(in):: ntotal
      type(particles), intent(out):: parts(:)
      integer:: i, j, k, n

      n = 0
      do i = 1, mp
         do j = 1, np
            do k = 1, op
               n = n + 1
               parts(n)%indloc = n
               parts(n)%indglob = n
               parts(n)%x(1) = (i - 0.5_f)*dxo
               parts(n)%x(2) = (j - 0.5_f)*dxo
               parts(n)%x(3) = (k - 0.5_f)*dxo
               parts(n)%vx(:) = 0_f
               parts(n)%itype = 1
               parts(n)%rho = irho
               parts(n)%p = 0_f
            end do
         end do
      end do

   end subroutine generate_real_part

   !==============================================================================================================================
   pure subroutine generate_virt_part(ntotal, nvirt, parts)

      implicit none
      integer, intent(in):: ntotal, nvirt
      type(particles), intent(out):: parts(:)
      integer:: i, j, k, n

      n = ntotal
      do i = 1 - nlayer, pp + nlayer
         do j = 1 - nlayer, qp + nlayer
            do k = 1 - nlayer, rp + nlayer
               if (i < 1 .or. i > pp .or. j < 1 .or. j > qp .or. k < 1 .or. k > rp) then
                  n = n + 1
                  parts(n)%indloc = n
                  parts(n)%indglob = n
                  parts(n)%x(1) = vxmin + (i - 0.5_f)*dxo
                  parts(n)%x(2) = vymin + (j - 0.5_f)*dxo
                  parts(n)%x(3) = vzmin + (k - 0.5_f)*dxo
                  parts(n)%vx(:) = 0._f
                  if (k < 1) parts(n)%itype = -1
                  if (k > rp) parts(n)%itype = -2
                  if (i < 1 .or. i > pp) parts(n)%itype = -3
                  if (j < 1 .or. j > qp) parts(n)%itype = -4
               end if
            end do
         end do
      end do

      parts(ntotal + 1:ntotal + nvirt)%rho = irho
      parts(ntotal + 1:ntotal + nvirt)%p = 0_f

   end subroutine generate_virt_part

   !==============================================================================================================================
    subroutine update_virt_part(ki, ntotal, nvirt, nghos, parts, niac, pairs, nexti, vw)

       implicit none
       integer, intent(in):: ki, ntotal, nvirt, nghos, niac, nexti(:)
       type(interactions), intent(in):: pairs(:)
       type(particles), intent(inout):: parts(:)
       real(f), intent(inout):: vw(:)
       integer:: i, j, k
       real(f):: tmp

       vw(:) = 0._f

       do i = 1, nvirt
          parts(ntotal + i)%rho = 0._f
          parts(ntotal + i)%vx(:) = 0._f
       end do

       do i = 1, ntotal + nvirt + nghos
         do k = nexti(i), nexti(i + 1) - 1
            j = pairs(k)%j
          if (parts(i)%itype < 0 .and. parts(j)%itype > 0) then
             tmp = mass*pairs(k)%w/parts(j)%rho
             vw(i - ntotal) = vw(i - ntotal) + tmp
             parts(i)%rho = parts(i)%rho + mass*pairs(k)%w
             select case (parts(i)%itype)
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
             end select
          else if (parts(j)%itype < 0 .and. parts(i)%itype > 0) then
             tmp = mass*pairs(k)%w/parts(i)%rho
             vw(j - ntotal) = vw(j - ntotal) + tmp
             parts(j)%rho = parts(j)%rho + mass*pairs(k)%w
             select case (parts(j)%itype)
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
             end select
          end if
         end do
       end do

       do i = 1, nvirt
          if (vw(i) > 0._f) then
             parts(ntotal + i)%rho = parts(ntotal + i)%rho/vw(i)
             parts(ntotal + i)%vx(:) = parts(ntotal + i)%vx(:)/vw(i)
          else
             parts(ntotal + i)%rho = irho
             parts(ntotal + i)%vx(:) = 0._f
          end if
          parts(ntotal + i)%p = rh0*c**2*((parts(ntotal + i)%rho/rh0)**gamma - 1_f)/gamma
       end do

    end subroutine update_virt_part

end module input_m
