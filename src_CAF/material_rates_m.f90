module material_rates_m

    use datatypes, only: particles
    use param, only: mass, dim, f
 
    private
    public:: art_visc, ext_force, int_force, con_density, art_visc_coeff, int_force_coeff
 
 contains
 
    !=================================================================================
    pure subroutine art_visc(p_i, p_j, dwdx, ardvxdti, ardvxdtj)
 
       use param, only: alpha, beta, etq, hsml, c
 
       implicit none
       type(particles), intent(in):: p_i, p_j
       real(f), intent(in):: dwdx(dim)
       real(f), intent(inout):: ardvxdti(dim), ardvxdtj(dim)
       real(f):: dx(dim), piv(dim), muv, vr, rr, mrho
 
       dx(:) = p_i%x(:) - p_j%x(:)
       vr = DOT_PRODUCT(p_i%vx(:) - p_j%vx(:), dx(:))
       if (vr > 0_f) vr = 0._f
 
       !Artificial viscous force only if v_ij * r_ij < 0
       rr = DOT_PRODUCT(dx(:), dx(:))
       muv = hsml*vr/(rr + hsml*hsml*etq*etq)
       mrho = 0.5_f*(p_i%rho + p_j%rho)
       piv = (beta*muv - alpha*c)*muv/mrho*dwdx(:)
 
       ardvxdti(:) = ardvxdti(:) - mass*piv(:)
       ardvxdtj(:) = ardvxdtj(:) + mass*piv(:)
 
    end subroutine art_visc
 
    !=================================================================================
    pure subroutine ext_force(p_i, p_j, dwdx, exdvxdti, exdvxdtj)
 
       use param, only: p1, p2, rr0, dd
 
       implicit none
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
    pure subroutine int_force(p_i, p_j, dwdx, indvxdti, indvxdtj)
 
       implicit none
       type(particles), intent(in):: p_i, p_j
       real(f), intent(in):: dwdx(dim)
       real(f), intent(inout):: indvxdti(dim), indvxdtj(dim)
       real(f):: h(dim)
 
       h = -(p_i%p/p_i%rho**2 + p_j%p/p_j%rho**2)*dwdx(:)
       indvxdti(:) = indvxdti(:) + mass*h(:)
       indvxdtj(:) = indvxdtj(:) - mass*h(:)
 
    end subroutine int_force
 
    !=================================================================================
    pure subroutine con_density(p_i, p_j, dwdx, codrhodti, codrhodtj)
 
       implicit none
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
    pure function art_visc_coeff(p_i, p_j) result(coeff)
 
       use param, only: alpha, beta, etq, hsml, c
 
       implicit none
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
    pure function int_force_coeff(prhoi, prhoj) result(coeff)
 
       implicit none
       real(f), intent(in):: prhoi, prhoj
       real(f):: coeff
 
       coeff = -mass*(prhoi + prhoj)
 
    end function int_force_coeff
 
 end module material_rates_m
 