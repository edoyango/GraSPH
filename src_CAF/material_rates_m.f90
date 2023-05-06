module material_rates_m

   use datatypes, only: particles
   use param, only: mass, dim, f

   private
   public:: con_density, art_visc_coeff, int_force_coeff

contains

   !=================================================================================
   pure subroutine ext_force(x_i, x_j, dwdx, exdvxdti, exdvxdtj)

      use param, only: p1, p2, rr0, dd

      implicit none
      real(f), intent(in):: x_i(dim), x_j(dim), dwdx(dim)
      real(f), intent(inout):: exdvxdti(dim), exdvxdtj(dim)
      real(f):: dx(dim), rr, f

      dx(:) = x_i(:) - x_j(:)
      rr = SQRT(SUM(dx(:)*dx(:)))

      if (rr .lt. rr0) then
         f = ((rr0/rr)**p1 - (rr0/rr)**p2)/rr**2
         exdvxdti(:) = exdvxdti(:) + dd*dx(:)*f
         exdvxdtj(:) = exdvxdtj(:) - dd*dx(:)*f
      end if

   end subroutine ext_force

   !=================================================================================
   pure subroutine con_density(vx_i, vx_j, dwdx, codrhodti, codrhodtj)

      implicit none
      real(f), intent(in):: vx_i(dim), vx_j(dim), dwdx(dim)
      real(f), intent(inout):: codrhodti, codrhodtj
      real(f):: dvx(dim), vcc

      dvx(:) = vx_i(:) - vx_j(:)

      vcc = DOT_PRODUCT(dvx(:), dwdx(:))

      codrhodti = codrhodti + mass*vcc
      codrhodtj = codrhodtj + mass*vcc

   end subroutine con_density

   !=================================================================================
   pure function art_visc_coeff(rho_i, vx_i, rho_j, vx_j, dx) result(coeff)

      use param, only: alpha, beta, etq, hsml, c

      implicit none
      real(f), intent(in):: dx(dim), rho_i, vx_i(dim), rho_j, vx_j(dim)
      real(f):: piv(dim), muv, vr, rr, mrho, coeff

      ! dx(:) = p_i%x(:) - p_j%x(:)
      vr = DOT_PRODUCT(vx_i(:) - vx_j(:), dx(:))
      if (vr > 0._f) vr = 0._f

      !Artificial viscous force only if v_ij * r_ij < 0
      rr = DOT_PRODUCT(dx(:), dx(:))
      muv = hsml*vr/(rr + hsml*hsml*etq*etq)
      mrho = 0.5_f*(rho_i + rho_j)
      coeff = -mass*(beta*muv - alpha*c)*muv/mrho

   end function art_visc_coeff

   !=================================================================================
   pure function int_force_coeff(prhoi, prhoj) result(coeff)

      implicit none
      real(f), intent(in):: prhoi, prhoj
      real(f):: coeff

      coeff = -mass*(prhoi + prhoj)

   end function int_force_coeff

end module material_rates_m
