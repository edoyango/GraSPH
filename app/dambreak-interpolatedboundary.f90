module grasph_monaghan1994_2

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles
    use weakly_compressible_particles, only: wcp => tait_eos_particles
    use weakly_compressible_interactions, only: fluid_sweeper
    use grasph_pair_sets, only: particle_interactions, base_sweeper
    use grasph_pair_interactions, only: artificial_viscosity_monaghan1994, continuity_density, isotropic_pressure_force, &
                                        repulsive_force
    use grasph_pairs, only: particle_pairs
    use grasph_particle_shifting, only: xsph_shifter

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.5_fp, g = -9.81_fp
    integer, parameter:: nfx = 25._fp/dx, nfy = 25._fp/dx, nbx = 75._fp/dx, nby = 40._fp/dx

    ! define how fluid particles interact with boundary
    type, extends(base_sweeper):: boundary_update_sweeper
    contains
        procedure:: sweep => boundary_update_sweep
    end type boundary_update_sweeper

contains

    subroutine boundary_update_sweep(self, pairs, ps_lhs, ps_rhs)
        class(boundary_update_sweeper), intent(in):: self
        type(particle_pairs), intent(in):: pairs
        class(base_particles), intent(inout):: ps_lhs
        class(base_particles), optional, intent(inout):: ps_rhs
        integer:: i, j, k
        real(fp):: mw, vw
        real(fp), allocatable:: wsum(:)

        if (present(ps_rhs)) then
            allocate (wsum(ps_rhs%size), source=0._fp)
        else
            error stop "Expected ps_rhs to be passed in."
        end if

        do i = 1, ps_rhs%size
            ps_rhs%ps(i)%rho = 0._fp
            ps_rhs%ps(i)%v(:) = 0._fp
        end do

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            mw = ps_lhs%ps(i)%mass*pairs%w(k)
            vw = mw/ps_lhs%ps(i)%rho
            wsum(j) = wsum(j) + vw
            ps_rhs%ps(j)%rho = ps_rhs%ps(j)%rho + mw
            ps_rhs%ps(j)%v(:) = ps_rhs%ps(j)%v(:) + ps_lhs%ps(i)%v(:)*vw
        end do

        do j = 1, ps_rhs%size
            if (wsum(j) > 0._fp) then
                ps_rhs%ps(j)%v(:) = -ps_rhs%ps(j)%v(:)/wsum(j)
                ps_rhs%ps(j)%rho = ps_rhs%ps(j)%rho/wsum(j)
            end if
        end do

    end subroutine boundary_update_sweep

end module grasph_monaghan1994_2

program main

    use grasph_monaghan1994_2

    use grasph_particles, only: particles_container
    use weakly_compressible_particles, only: wcp => tait_eos_particles
    use grasph_pair_sets, only: particle_interactions
    use grasph_time_integration, only: leap_frog_time_integration
    use grasph_kernels, only: grasph_cubic_bspline_kernel

    implicit none
    type(particles_container):: ps(2)
    type(particle_interactions):: pic(2)
    type(grasph_cubic_bspline_kernel):: kernel
    type(fluid_sweeper):: sweeper
    type(boundary_update_sweeper):: boundary_sweeper
    type(xsph_shifter):: shifter
    integer:: i, j, k, nlayer, nvirt

    ! declare particles - fluid and boundary
    allocate (wcp::ps(1)%p)
    allocate (wcp::ps(2)%p)

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    ! init fluid particles
    select type (ps => ps(1)%p) ! specialise for weakly compressible particles
    class is (wcp)
        call ps%init(n=2500, name="fluid", rho_ref=1000._fp)
        call ps%register_x%register(ps%ps(1), "x", ps%ps(1)%x, ps%ps(1)%v)
        call ps%register_v%register(ps%ps(1), "v", ps%ps(1)%v, ps%ps(1)%dvxdt)
        call ps%register_v%register(ps%ps(1), "rho", ps%ps(1)%rho, ps%ps(1)%drhodt)
    end select
    do i = 0, nfx - 1
        do j = 0, nfy - 1
            k = j*nfx + i + 1
            ps(1)%p%ps(k)%id = k
            ps(1)%p%ps(k)%type = 1 ! not sure if type is needed anymore
            ps(1)%p%ps(k)%x(1) = (i + 0.5_fp)*dx
            ps(1)%p%ps(k)%x(2) = (j + 0.5_fp)*dx
            ps(1)%p%ps(k)%rho = 1000._fp
            ps(1)%p%ps(k)%mass = 1000._fp*dx*dx
            ps(1)%p%ps(k)%c = 10._fp*sqrt(490.5_fp) ! 10*sqrt(2gH)
            ps(1)%p%ps(k)%v(:) = 0._fp
        end do
    end do

    ! init boundary particles
    ! use base_init since we're using the base type
    ! only need to initialize metadata and position as only position is used to calculate repulsive force
    nlayer = ceiling(kernel%cutoff/dx)
    nvirt = nlayer*2*(nbx + nby) + 4*nlayer*nlayer
    select type (ps => ps(2)%p)
    type is (wcp)
        call ps%init(n=nvirt, name="boundary", rho_ref=1000._fp)
    end select
    k = 0
    ! bottom layer and corners
    do i = -nlayer, nbx + nlayer - 1
        do j = 0, nlayer - 1
            k = k + 1
            ps(2)%p%ps(k)%x(1) = (i + 0.5_fp)*dx
            ps(2)%p%ps(k)%x(2) = -(j + 0.5_fp)*dx
        end do
    end do
    ! top layer and corners
    do i = -nlayer, nbx + nlayer - 1
        do j = 0, nlayer - 1
            k = k + 1
            ps(2)%p%ps(k)%x(1) = (i + 0.5_fp)*dx
            ps(2)%p%ps(k)%x(2) = 40._fp + (j + 0.5_fp)*dx
        end do
    end do
    ! left wall
    do j = 0, nby - 1
        do i = 0, nlayer - 1
            k = k + 1
            ps(2)%p%ps(k)%x(1) = -(i + 0.5_fp)*dx
            ps(2)%p%ps(k)%x(2) = (j + 0.5_fp)*dx
        end do
    end do
    ! right wall
    do j = 0, nby - 1
        do i = 0, nlayer - 1
            k = k + 1
            ps(2)%p%ps(k)%x(1) = 75._fp + (i + 0.5_fp)*dx
            ps(2)%p%ps(k)%x(2) = (j + 0.5_fp)*dx
        end do
    end do
    do i = 1, k
        ps(2)%p%ps(i)%id = i
        ps(2)%p%ps(i)%type = -1
        ps(2)%p%ps(i)%rho = 1000._fp
        ps(2)%p%ps(i)%mass = 1000._fp*dx*dx
        ps(2)%p%ps(i)%c = 10._fp*sqrt(490.5_fp)
    end do

    sweeper%artvisc_alpha = 0.01_fp
    sweeper%artvisc_beta = 0._fp
    sweeper%g = g
    sweeper%h = 1.2_fp*dx
    sweeper%update_rhs = .false.

    shifter%epsilon = 0.5_fp
    shifter%update_rhs = .false.

    ! init interactions
    call pic(1)%init(30, ps(1)%p, sweeper=sweeper, shifter=shifter)
    sweeper%initialize = .false. ! second sweeper doesn't need to zero acceleation arrays
    call pic(2)%init(30, ps(1)%p, ps(2)%p, prologue_sweeper=boundary_sweeper, sweeper=sweeper, shifter=shifter)

    call leap_frog_time_integration(100000, 1000, 1000, ps, pic, 0.05_fp, kernel, "/home/edwardy/test", output_comp_level=4)

end program main
