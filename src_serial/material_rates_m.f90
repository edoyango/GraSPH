module material_rates_m

   use datatypes, only: particles, interactions
   use param, only: mass, dim, f

   public:: con_density, art_visc_coeff, int_force_coeff

contains

   !=================================================================================
   pure subroutine ext_force(ki, p_i, p_j, dwdx, exdvxdti, exdvxdtj)
      ! This exists more as reference in case it ever needs to be implemted again

      use param, only: p1, p2, rr0, dd

      implicit none
      integer, intent(in):: ki
      type(particles), intent(in):: p_i, p_j
      real(f), intent(in):: dwdx(dim)
      real(f), intent(inout):: exdvxdti(dim), exdvxdtj(dim)
      real(f):: dx(dim), rr, f

      dx(:) = p_i%x(:) - p_j%x(:)
      rr = SQRT(SUM(dx(:)*dx(:)))

      if (rr .lt. rr0) then
         f = ((rr0/rr)**p1 - (rr0/rr)**p2)/rr**2
         exdvxdti(:) = exdvxdti(:) + dd*dx(:)*f
         exdvxdtj(:) = exdvxdtj(:) - dd*dx(:)*f
      end if

   end subroutine ext_force

   !=================================================================================
   pure subroutine con_density(ki, p_i, p_j, dwdx, codrhodti, codrhodtj)

      implicit none
      integer, intent(in):: ki
      type(particles), intent(in):: p_i, p_j
      real(f), intent(in):: dwdx(dim)
      real(f), intent(inout):: codrhodti, codrhodtj
      real(f):: dvx(dim), vcc

      dvx(:) = p_i%vx(:) - p_j%vx(:)

      vcc = DOT_PRODUCT(dvx(:), dwdx(:))

      codrhodti = codrhodti + mass*vcc
      codrhodtj = codrhodtj + mass*vcc

   end subroutine con_density

   !=================================================================================
   pure function art_visc_coeff(ki, p_i, p_j) result(coeff)

      use param, only: alpha, beta, etq, hsml, c

      implicit none
      integer, intent(in):: ki
      type(particles), intent(in):: p_i, p_j
      real(f):: dx(dim), piv(dim), muv, vr, rr, mrho, coeff

      dx(:) = p_i%x(:) - p_j%x(:)
      vr = DOT_PRODUCT(p_i%vx(:) - p_j%vx(:), dx(:))
      if (vr > 0_f) vr = 0_f

      !Artificial viscous force only if v_ij * r_ij < 0
      rr = DOT_PRODUCT(dx(:), dx(:))
      muv = hsml*vr/(rr + hsml*hsml*etq*etq)
      mrho = 0.5_f*(p_i%rho + p_j%rho)
      coeff = -mass*(beta*muv - alpha*c)*muv/mrho

   end function art_visc_coeff

   !=================================================================================
   pure function int_force_coeff(ki, prhoi, prhoj) result(coeff)

      implicit none
      integer, intent(in):: ki
      real(f), intent(in):: prhoi, prhoj
      real(f):: coeff

      coeff = -mass*(prhoi + prhoj)

   end function int_force_coeff

end module material_rates_m
