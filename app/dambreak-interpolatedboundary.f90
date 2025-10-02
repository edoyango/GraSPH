module grasph_monaghan1994_2

    use grasph_constants, only: fp
    use grasph_particles, only: particle_system_t
    use weakly_compressible_interactions, only: fluid_sweeper
    use grasph_pair_sets, only: particle_interactions, base_sweeper
    use grasph_pair_interactions, only: artificial_viscosity_monaghan1994, continuity_density, isotropic_pressure_force, &
                                        repulsive_force
    use grasph_pairs, only: particle_pairs
    use grasph_particle_shifting, only: xsph_shifter

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.5_fp, g = -9.81_fp, rho0 = 1000._fp
    integer, parameter:: nfx = 25._fp/dx, nfy = 25._fp/dx, nbx = 75._fp/dx, nby = 40._fp/dx

    ! define how fluid particles interact with boundary
    type, extends(base_sweeper):: boundary_update_sweeper
    contains
        procedure:: sweep => boundary_update_sweep
    end type boundary_update_sweeper

contains

    subroutine boundary_update_sweep(self, pairs, psys_lhs, psys_rhs)
        class(boundary_update_sweeper), intent(in):: self
        type(particle_pairs), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        integer:: i, j, k
        real(fp):: mw, vw
        real(fp), allocatable:: wsum(:)

        if (present(psys_rhs)) then
            allocate (wsum(psys_rhs%size), source=0._fp)
        else
            error stop "Expected psys_rhs to be passed in."
        end if

        do i = 1, psys_rhs%size
            psys_rhs%particles(i)%rho = 0._fp
            psys_rhs%particles(i)%v(:) = 0._fp
        end do

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            mw = psys_lhs%particles(i)%mass*pairs%w(k)
            vw = mw/psys_lhs%particles(i)%rho
            wsum(j) = wsum(j) + vw
            psys_rhs%particles(j)%rho = psys_rhs%particles(j)%rho + mw
            psys_rhs%particles(j)%v(:) = psys_rhs%particles(j)%v(:) + psys_lhs%particles(i)%v(:)*vw
        end do

        do j = 1, psys_rhs%size
            if (wsum(j) > 0._fp) then
                psys_rhs%particles(j)%v(:) = -psys_rhs%particles(j)%v(:)/wsum(j)
                psys_rhs%particles(j)%rho = psys_rhs%particles(j)%rho/wsum(j)
            end if
        end do

    end subroutine boundary_update_sweep

end module grasph_monaghan1994_2

program main

    use grasph_monaghan1994_2

    use grasph_particles, only: particle_system_t
    use weakly_compressible_particles, only: tait_eos_state_updater, eos_particle
    use grasph_pair_sets, only: particle_interactions
    use grasph_time_integration, only: leap_frog_time_integration
    use grasph_kernels, only: grasph_cubic_bspline_kernel

    implicit none
    type(particle_system_t):: psys(2)
    type(particle_interactions):: pic(2)
    type(grasph_cubic_bspline_kernel):: kernel
    type(fluid_sweeper):: sweeper
    type(boundary_update_sweeper):: boundary_sweeper
    type(xsph_shifter):: shifter
    type(tait_eos_state_updater):: state_updater
    type(eos_particle):: ps_template
    integer:: i, j, k, nlayer, nvirt

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    ! init fluid particles
    state_updater%rho_ref = rho0
    call psys(1)%base_init(n=2500, name="fluid", state_updater=state_updater, particle_template=ps_template)

    ! register variables for time-update
    call psys(1)%register_x%register(psys(1)%particles(1), "x", psys(1)%particles(1)%x, psys(1)%particles(1)%v)
    call psys(1)%register_v%register(psys(1)%particles(1), "v", psys(1)%particles(1)%v, psys(1)%particles(1)%dvxdt)
    call psys(1)%register_v%register(psys(1)%particles(1), "rho", psys(1)%particles(1)%rho, psys(1)%particles(1)%drhodt)

    ! register variables for io
    select type (p => psys(1)%particles)
    class is (eos_particle)
        call psys(1)%register_io%register_variable(p(1), "x", p(1)%x)
        call psys(1)%register_io%register_variable(p(1), "v", p(1)%v)
        call psys(1)%register_io%register_variable(p(1), "rho", p(1)%rho)
        call psys(1)%register_io%register_variable(p(1), "mass", p(1)%mass)
        call psys(1)%register_io%register_variable(p(1), "c", p(1)%c)
        call psys(1)%register_io%register_variable(p(1), "dvxdt", p(1)%dvxdt)
        call psys(1)%register_io%register_variable(p(1), "drhodt", p(1)%drhodt)
        call psys(1)%register_io%register_variable(p(1), "p", p(1)%p)
    class default
        error stop "Expected eos_particle for psys(1)%p."
    end select

    do i = 0, nfx - 1
        do j = 0, nfy - 1
            k = j*nfx + i + 1
            psys(1)%particles(k)%id = k
            psys(1)%particles(k)%type = 1 ! not sure if type is needed anymore
            psys(1)%particles(k)%x(1) = (i + 0.5_fp)*dx
            psys(1)%particles(k)%x(2) = (j + 0.5_fp)*dx
            psys(1)%particles(k)%rho = rho0
            psys(1)%particles(k)%mass = rho0*dx*dx
            psys(1)%particles(k)%c = 10._fp*sqrt(490.5_fp) ! 10*sqrt(2gH)
            psys(1)%particles(k)%v(:) = 0._fp
        end do
    end do

    ! init boundary particles
    ! use base_init since we're using the base type
    ! only need to initialize metadata and position as only position is used to calculate repulsive force
    nlayer = ceiling(kernel%cutoff/dx)
    nvirt = nlayer*2*(nbx + nby) + 4*nlayer*nlayer
    call psys(2)%base_init(n=nvirt, name="boundary", state_updater=state_updater, particle_template=ps_template)

    ! register variables for io
    select type (p => psys(2)%particles)
    class is (eos_particle)
        call psys(2)%register_io%register_variable(p(1), "x", p(1)%x)
        call psys(2)%register_io%register_variable(p(1), "v", p(1)%v)
        call psys(2)%register_io%register_variable(p(1), "rho", p(1)%rho)
        call psys(2)%register_io%register_variable(p(1), "mass", p(1)%mass)
        call psys(2)%register_io%register_variable(p(1), "c", p(1)%c)
        call psys(2)%register_io%register_variable(p(1), "dvxdt", p(1)%dvxdt)
        call psys(2)%register_io%register_variable(p(1), "drhodt", p(1)%drhodt)
        call psys(2)%register_io%register_variable(p(1), "p", p(1)%p)
    class default
        error stop "Expected eos_particle for psys(1)%p."
    end select
    k = 0
    ! bottom layer and corners
    do i = -nlayer, nbx + nlayer - 1
        do j = 0, nlayer - 1
            k = k + 1
            psys(2)%particles(k)%x(1) = (i + 0.5_fp)*dx
            psys(2)%particles(k)%x(2) = -(j + 0.5_fp)*dx
        end do
    end do
    ! top layer and corners
    do i = -nlayer, nbx + nlayer - 1
        do j = 0, nlayer - 1
            k = k + 1
            psys(2)%particles(k)%x(1) = (i + 0.5_fp)*dx
            psys(2)%particles(k)%x(2) = 40._fp + (j + 0.5_fp)*dx
        end do
    end do
    ! left wall
    do j = 0, nby - 1
        do i = 0, nlayer - 1
            k = k + 1
            psys(2)%particles(k)%x(1) = -(i + 0.5_fp)*dx
            psys(2)%particles(k)%x(2) = (j + 0.5_fp)*dx
        end do
    end do
    ! right wall
    do j = 0, nby - 1
        do i = 0, nlayer - 1
            k = k + 1
            psys(2)%particles(k)%x(1) = 75._fp + (i + 0.5_fp)*dx
            psys(2)%particles(k)%x(2) = (j + 0.5_fp)*dx
        end do
    end do
    do i = 1, k
        psys(2)%particles(i)%id = i
        psys(2)%particles(i)%type = -1
        psys(2)%particles(i)%rho = rho0
        psys(2)%particles(i)%mass = rho0*dx*dx
        psys(2)%particles(i)%c = 10._fp*sqrt(490.5_fp)
    end do

    sweeper%artvisc_alpha = 0.01_fp
    sweeper%artvisc_beta = 0._fp
    sweeper%g = g
    sweeper%h = 1.2_fp*dx
    sweeper%update_rhs = .false.

    shifter%epsilon = 0.5_fp
    shifter%update_rhs = .false.

    ! init interactions
    call pic(1)%init(30, psys(1), sweeper=sweeper, shifter=shifter)
    sweeper%initialize = .false. ! second sweeper doesn't need to zero acceleation arrays
    call pic(2)%init(30, psys(1), psys(2), prologue_sweeper=boundary_sweeper, sweeper=sweeper, shifter=shifter)

    call leap_frog_time_integration(100000, 1000, 1000, psys, pic, 0.05_fp, kernel, "/home/edwardy/test", output_comp_level=4)

end program main
