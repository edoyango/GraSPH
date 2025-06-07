module grasph_pair_interactions

    use grasph_constants, only: fp, ndims

    implicit none

contains

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

    end subroutine artificial_viscosity

    pure subroutine continuity_density(vi, vj, massi, massj, drhodti, drhodtj, dwdx)

        real(fp), intent(in):: vi(ndims), vj(ndims), massi, massj, dwdx(ndims)
        real(fp), intent(inout):: drhodti, drhodtj
        real(fp):: vcc

        vcc = dot_product(vi(:)-vj(:), dwdx(:))

        drhodti = drhodti + massj*vcc
        drhodtj = drhodtj + massi*vcc

    end subroutine continuity_density

    pure subroutine isotropic_pressure_force(pi, pj, rhoi, rhoj, massi, massj, dvxdti, dvxdtj, dwdx)

        real(fp), intent(in):: pi, pj, rhoi, rhoj, massi, massj, dwdx(ndims)
        real(fp), intent(inout):: dvxdti(ndims), dvxdtj(ndims)
        real(fp):: h(ndims)

        h = -(pi/(rhoi*rhoi) + pj/(rhoj*rhoj))*dwdx(:)
        dvxdti(:) = dvxdti(:) + mass(j)*h(:)
        dvxdtj(:) = dvxdtj(:) - mass(i)*h(:)

    end subroutine isotropic_pressure_force

end module grasph_pair_interactions