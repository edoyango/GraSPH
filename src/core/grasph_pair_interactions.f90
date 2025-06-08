module grasph_pair_interactions

    use grasph_constants, only: fp

    implicit none

contains

    pure subroutine artificial_viscosity_monaghan1994(ndims, xi, xj, vi, vj, rhoi, rhoj, hsmli, hsmlj, ci, &
        cj, massi, massj, dvxdti, dvxdtj, dwdx, alpha, beta)

        integer, intent(in):: ndims
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

    pure subroutine continuity_density(ndims, vi, vj, massi, massj, drhodti, drhodtj, dwdx)

        integer, intent(in):: ndims
        real(fp), intent(in):: vi(ndims), vj(ndims), massi, massj, dwdx(ndims)
        real(fp), intent(inout):: drhodti, drhodtj
        real(fp):: vcc

        vcc = dot_product(vi(:)-vj(:), dwdx(:))

        drhodti = drhodti + massj*vcc
        drhodtj = drhodtj + massi*vcc

    end subroutine continuity_density

    pure subroutine isotropic_pressure_force(ndims, pi, pj, rhoi, rhoj, massi, massj, dvxdti, dvxdtj, dwdx)

        integer, intent(in):: ndims
        real(fp), intent(in):: pi, pj, rhoi, rhoj, massi, massj, dwdx(ndims)
        real(fp), intent(inout):: dvxdti(ndims), dvxdtj(ndims)
        real(fp):: h(ndims)

        h = -(pi/(rhoi*rhoi) + pj/(rhoj*rhoj))*dwdx(:)
        dvxdti(:) = dvxdti(:) + massj*h(:)
        dvxdtj(:) = dvxdtj(:) - massi*h(:)

    end subroutine isotropic_pressure_force

    pure subroutine repulsive_force(ndims, cutoff, ci, xi, xj, dvxdti)
        integer, intent(in):: ndims
        real(fp), intent(in):: cutoff, ci, xi(ndims), xj(ndims)
        real(fp), intent(inout):: dvxdti(ndims)
        real(fp):: dx(ndims), f, r
        integer, parameter:: p1 = 4, p2 = 2

        dx(:) = xi(:) - xj(:)
        r = sqrt(sum(dx(:)*dx(:)))
        if (r < cutoff) then
            f = ((cutoff/r)**p1 - (cutoff/r)**p2)/(r*r)
            dvxdti(:) = dvxdti(:) + 0.01_fp*ci*ci*f*dx(:)
        endif
    end subroutine repulsive_force


end module grasph_pair_interactions