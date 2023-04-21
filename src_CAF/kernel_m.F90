module kernel_m

   use param, only: f

   private
   public:: kernel_k, kernel_w, kernel_dwdx

contains

   !====================================================================================================================
   pure function kernel_w(r, thsml) result(tw)
      ! Contains the kernels
      ! NB: DIM is a gfortran function where DIM(x,y) = MAX(0,x-y). DIM is slightly faster.

      use param, only: skf, dims => dim, pi

      implicit none
      real(f), intent(in):: r, thsml
      real(f):: q, factor, tw

      q = r/thsml

      ! SELECT CASE (SKF)
      ! CASE (1) ! cubic
#if defined(KERNEL_CUBIC)
#ifdef D2
         factor = 10._f/(7._f*pi*thsml*thsml)
#else
         factor = 1._f/(pi*thsml*thsml)
#endif
         tw = factor*(0.25_f*DIM(2._f, q)**3 - DIM(1._f, q)**3)
#elif defined(KERNEL_QUARTIC)
#ifdef D2
         factor = 96._f/(1199._f*pi*thsml*thsml)
#else
         factor = 1._f/(20._f*pi*thsml*thsml*thsml)
#endif
         tw = factor*(DIM(2.5_f, q)**4 - 5._f*DIM(1.5_f, q)**4 + 10._f*DIM(0.5_f, q)**4)
#elif defined(KERNEL_QUINTIC)
#ifdef D2
         factor = 7._f/(478._f*pi*thsml*thsml)
#else
         factor = 1._f/(120._f*pi*thsml*thsml*thsml)
#endif
         tw = factor*(DIM(3._f, q)**5 - 6._f*DIM(2._f, q)**5 + 15._f*DIM(1._f, q)**5)
#elif defined(KERNEL_WQUINTIC_C4)
#ifdef D2
         factor = 3._f/(1024._f*pi*thsml*thsml)
#else
         factor = 165._f/(65536._f*pi*thsml*thsml*thsml)
#endif
         tw = factor*DIM(2._f, q)**6*(35._f*q**2 + 36._f*q + 12._f)
      ! CASE (6) ! Wenland Quintic C6
#elif defined(KERNEL_WQUINTIC_C6)
#ifdef D2
         factor = 39._f/(14336._f*pi*thsml*thsml)
#else
         factor = 1365._f/(524288._f*pi*thsml*thsml*thsml)
#endif
         tw = factor*DIM(2._f, q)**8*(16._f*q**3 + 25_f*q**2 + 16._f*q + 4._f)
#elif defined(KERNEL_GUASSIAN)
         factor = 1._f/(thsml**dims*pi**(0.5_f*dims))
         if (q .ge. 0._f .and. q .le. 3._f) then
            tw = factor*exp(-q*q)
         end if
#else
#ifdef D2
         factor = 7._f/(64._f*pi*thsml*thsml)
#else
         factor = 21._f/(256._f*pi*thsml*thsml*thsml)
#endif
         tw = factor*DIM(2._f, q)**4*(2._f*q + 1._f)
#endif

   end function kernel_w

   !====================================================================================================================
   pure function kernel_dwdx(dx, thsml) result(tdwdx)
      ! Contains the kernels
      ! NB: DIM is a gfortran function where DIM(x,y) = MAX(0,x-y). DIM is slightly faster.

      use param, only: skf, dims => dim, pi

      implicit none
      real(f), intent(in):: dx(dims), thsml
      real(f):: q, factor, tdwdx(dims), r

      r = sqrt(sum(dx*dx))
      q = r/thsml

      ! SELECT CASE (SKF)
      ! CASE (1) ! cubic
#if defined(KERNEL_CUBIC)
#ifdef D2
         factor = 10._f/(7._f*pi*thsml*thsml)
#else
         factor = 1._f/(pi*thsml*thsml)
#endif
         tdwdx = -factor*3._f*(0.25_f*DIM(2._f, q)**2 - DIM(1._f, q)**2)*dx(:)/(r*thsml)
      ! CASE (2) ! quartic
#elif defined(KERNEL_QUARTIC)
#ifdef D2
         factor = 96._f/(1199._f*pi*thsml*thsml)
#else
         factor = 1._f/(20._f*pi*thsml*thsml*thsml)
#endif
         tdwdx(:) = -factor*4._f*(DIM(2.5_f, q)**3 - 5._f*DIM(1.5_f, q)**3 + 10._f*DIM(0.5_f, q)**3)*dx(:)/(r*thsml)
      ! CASE (3) ! quintic
#elif defined(KERNEL_QUINTIC)
#ifdef D2
         factor = 7._f/(478._f*pi*thsml*thsml)
#else
         factor = 1._f/(120._f*pi*thsml*thsml*thsml)
#endif
         tdwdx(:) = -factor*5._f*(DIM(3._f, q)**4 - 6._f*DIM(2._f, q)**4 + 15._f*DIM(1._f, q)**4)*dx(:)/(r*thsml)
      ! CASE (5) ! Wenland Quintic C4
#elif defined(KERNEL_WQUINTIC_C4)
#ifdef D2
         factor = 3._f/(1024._f*pi*thsml*thsml)
#else
         factor = 165._f/(65536._f*pi*thsml*thsml*thsml)
#endif
         tdwdx(:) = -factor*56._f*q*DIM(2._f, q)**5*(5._f*q + 2._f)*dx(:)/(r*thsml)
      ! CASE (6) ! Wenland Quintic C6
#elif defined(KERNEL_WQUINTIC_C6)
#ifdef D2
         factor = 39._f/(14336._f*pi*thsml*thsml)
#else
         factor = 1365._f/(524288._f*pi*thsml*thsml*thsml)
#endif
         tdwdx(:) = -factor*22._f*q*(8._f*q**2 + 7._f*q + 2._f)*DIM(2._f, q)**7*dx(:)/(r*thsml)
      ! CASE (7) ! gaussian
#elif defined(KERNEL_GUASSIAN)
         factor = 1._f/(thsml**dims*pi**(0.5_f*dims))
         if (q .ge. 0._f .and. q .le. 3._f) then
            tdwdx(:) = factor*exp(-q*q)*2._f*dx(:)/(thsml*thsml)
         end if
      ! CASE (4) ! Wenland Quintic C2
#else
#ifdef D2
         factor = 7._f/(64._f*pi*thsml*thsml)
#else
         factor = 21._f/(256._f*pi*thsml*thsml*thsml)
#endif
         tdwdx(:) = -factor*10._f*q*DIM(2._f, q)**3*dx(:)/(r*thsml)
#endif
      ! END SELECT

   end function kernel_dwdx

   !====================================================================================================================
   pure function kernel_k(skf) result(scale_k)
      ! setting k parameter for kernel radius (kh)

      implicit none
      integer, intent(in):: skf
      real(f):: scale_k

      select case (skf)
      case default
         scale_k = 2._f
      case (2)
         scale_k = 2.5_f
      case (3, 7)
         scale_k = 3._f
      end select

   end function kernel_k

end module kernel_m
