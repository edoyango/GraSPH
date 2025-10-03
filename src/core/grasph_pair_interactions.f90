!> @file grasph_pair_interactions.f90
!> @brief Module containing subroutines describing common particle interactions.
!> @author Edward Yang
!> @date 2025-06-09
module grasph_pair_interactions_m

    use grasph_constants_m, only: fp, ndims

    implicit none

contains

    !> @brief Artificial viscosity as described in http://dx.doi.org/10.1006/jcph.1994.1034
    !> @param xi The LHS particle's position.
    !> @param xj The RHS "                 ".
    !> @param vi The LHS particle's velocity.
    !> @param vj The RHS "                 ".
    !> @param rhoi The LHS particle's density.
    !> @param rhoj The RHS "                ".
    !> @param hsmli The LHS particle's smoothing length.
    !> @param hsmlj The RHS "                         ".
    !> @param ci The LHS particle's local speed of sound.
    !> @param cj The RHS particle's "                  ".
    !> @param massi The LHS particle's mass.
    !> @param massj The RHS "             ".
    !> @param dvxdti The LHS particle's acceleration.
    !> @param dvxdtj The RHS "                     ".
    !> @param dwdx The kernel gradient, relative to the LHS particle.
    !> @param alpha The viscous damping coefficient.
    !> @param beta Another damping coefficient.
    pure subroutine artificial_viscosity_monaghan1994(xi, xj, vi, vj, rhoi, rhoj, hsmli, hsmlj, ci, &
                                                      cj, massi, massj, dvxdti, dvxdtj, dwdx, alpha, beta)

        real(fp), intent(in):: xi(ndims), xj(ndims), vi(ndims), vj(ndims), rhoi, rhoj, hsmli, hsmlj, &
                               massi, massj, ci, cj, dwdx(ndims), alpha, beta
        real(fp), intent(inout):: dvxdti(ndims), dvxdtj(ndims)
        real(fp):: dx(ndims), dv(ndims), vr, muv, mhsml, mc, mrho, piv(ndims)

        dx(:) = xi(:) - xj(:)
        dv(:) = vi(:) - vj(:)
        vr = sum(dx(:)*dv(:))
        ! apply viscous damping force only for divergent particles
        ! set vr to 0._fp as it's faster than using an if statement
        vr = min(vr, 0._fp)
        mhsml = 0.5_fp*hsmli*hsmlj
        muv = mhsml*vr/(sum(dx(:)*dx(:)) + mhsml*mhsml*0.01_fp)
        mc = 0.5_fp*ci*cj
        mrho = 0.5_fp*rhoi*rhoj
        piv(:) = (beta*muv - alpha*mc)*muv/mrho*dwdx(:)

        dvxdti(:) = dvxdti(:) - massj*piv(:)
        dvxdtj(:) = dvxdtj(:) + massi*piv(:)

    end subroutine artificial_viscosity_monaghan1994

    !> @brief Continuity density as described in http://dx.doi.org/10.1006/jcph.1994.1034
    !> @param vi The LHS particle's velocity.
    !> @param vj The RHS "                 ".
    !> @param massi The LHS particle's mass.
    !> @param massj The RHS "             ".
    !> @param drhodti The LHS particle's density rate of change.
    !> @param drhodtj The RHS "                               ".
    !> @param dwdx The kernel gradient, relative to the LHS particle.
    pure subroutine continuity_density(vi, vj, massi, massj, drhodti, drhodtj, dwdx)

        real(fp), intent(in):: vi(ndims), vj(ndims), massi, massj, dwdx(ndims)
        real(fp), intent(inout):: drhodti, drhodtj
        real(fp):: vcc

        vcc = dot_product(vi(:) - vj(:), dwdx(:))

        drhodti = drhodti + massj*vcc
        drhodtj = drhodtj + massi*vcc

    end subroutine continuity_density

    !> @brief Isotropic pressure force as described in http://dx.doi.org/10.1006/jcph.1994.1034
    !> @param pi The LHS particle's pressure.
    !> @param pj The RHS "                 ".
    !> @param rhoi The LHS particle's density.
    !> @param rhoj The RHS "                ".
    !> @param massi The LHS particle's mass.
    !> @param massj The RHS "             ".
    !> @param dvxdti The LHS particle's acceleration.
    !> @param dvxdtj The RHS "                     ".
    !> @param dwdx The kernel gradient, relative to the LHS particle.
    pure subroutine isotropic_pressure_force(pi, pj, rhoi, rhoj, massi, massj, dvxdti, dvxdtj, dwdx)

        real(fp), intent(in):: pi, pj, rhoi, rhoj, massi, massj, dwdx(ndims)
        real(fp), intent(inout):: dvxdti(ndims), dvxdtj(ndims)
        real(fp):: h(ndims)

        h = -(pi/(rhoi*rhoi) + pj/(rhoj*rhoj))*dwdx(:)
        dvxdti(:) = dvxdti(:) + massj*h(:)
        dvxdtj(:) = dvxdtj(:) - massi*h(:)

    end subroutine isotropic_pressure_force

    !> @brief repulsive used for boundary particles as described in http://dx.doi.org/10.1006/jcph.1994.1034.
    !>        Unlike other interactions, this assumes the "real" particles are always on the LHS and the RHS
    !>        particles aren't evolved using the governing equations.
    !> @param cutoff The cut-off distance to apply the repulsive force.
    !> @param ci The LHS particle's local speed of sound.
    !> @param xi The LHS particle's position.
    !> @param xj The RHS "                 ".
    !> @param dvxdti The LHS particle's acceleration.
    pure subroutine repulsive_force(cutoff, ci, xi, xj, dvxdti)
        real(fp), intent(in):: cutoff, ci, xi(ndims), xj(ndims)
        real(fp), intent(inout):: dvxdti(ndims)
        real(fp):: dx(ndims), f, r
        integer, parameter:: p1 = 4, p2 = 2

        dx(:) = xi(:) - xj(:)
        r = sqrt(sum(dx(:)*dx(:)))
        if (r < cutoff) then
            f = ((cutoff/r)**p1 - (cutoff/r)**p2)/(r*r)
            dvxdti(:) = dvxdti(:) + 0.01_fp*ci*ci*f*dx(:)
        end if
    end subroutine repulsive_force

end module grasph_pair_interactions_m
