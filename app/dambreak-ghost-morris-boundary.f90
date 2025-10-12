module grasph_dambreak_ghost_morris_boundary_m

    use grasph_constants_m, only: fp, ndims
    use grasph_particle_m, only: base_particle_t
    use grasph_particle_system_m, only: particle_system_t, base_state_updater_t
    use weakly_compressible_particles_m, only: eos_particle_t
    use weakly_compressible_interactions_m, only: fluid_sweeper_t
    use grasph_system_interactions_m, only: base_sweeper_t
    use grasph_pair_interactions_m, only: artificial_viscosity_monaghan1994, continuity_density, isotropic_pressure_force, &
                                          repulsive_force
    use grasph_pairs_m, only: particle_pairs_t
    use grasph_particle_shifting_m, only: xsph_shifter_t

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.5_fp, g = -9.81_fp, rho0 = 1000._fp

    ! define how fluid particles interact with boundary
    type, extends(base_sweeper_t):: boundary_update_sweeper_t
    contains
        procedure:: sweep => boundary_update_sweep
    end type boundary_update_sweeper_t

    type, extends(eos_particle_t):: eos_ghost_particle_t
        class(eos_particle_t), pointer:: original
    end type eos_ghost_particle_t

    type, extends(base_sweeper_t):: ghost_timestep_setuper_t
        real(fp):: cutoff
        real(fp):: surface_normal(ndims), point(ndims)
    contains
        procedure:: sweep => ghost_timestep_setup_sweep
    end type

    type, extends(base_state_updater_t):: ghost_state_updater_t
        real(fp):: surface_normal(ndims)
    contains
        procedure:: update_state => ghost_state_update
    end type ghost_state_updater_t

contains

    subroutine boundary_update_sweep(self, pairs, psys_lhs, psys_rhs)
        class(boundary_update_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
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

    subroutine ghost_timestep_setup_sweep(self, pairs, psys_lhs, psys_rhs)
        class(ghost_timestep_setuper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        integer:: i
        class(eos_particle_t), pointer:: ps_real(:)
        class(eos_ghost_particle_t), pointer:: ps_ghost(:)
        real(fp):: dx(ndims), dr

        if (.not. present(psys_rhs)) error stop "Expected psys_rhs to be passed in."

        select type (ps => psys_lhs%particles)
        class is (eos_particle_t)
            ps_real => ps
        class default
            error stop "Expected psys_lhs to be eos_particle_t."
        end select

        select type (ps => psys_rhs%particles)
        class is (eos_ghost_particle_t)
            ps_ghost => ps
        class default
            error stop "Expected psys_rhs to be eos_ghost_particle_t."
        end select

        psys_rhs%size = 0

        do i = 1, psys_lhs%size
            dr = dot_product(ps_real(i)%x(:) - self%point(:), self%surface_normal(:))
            ! ghost particle if within cutoff and inside the modelled region.
            if (abs(dr) <= self%cutoff .and. dr > 0._fp) then
                psys_rhs%size = psys_rhs%size + 1
                ps_ghost(psys_rhs%size)%id = psys_rhs%size
                ps_ghost(psys_rhs%size)%original => ps_real(i)
                ps_ghost(psys_rhs%size)%x(:) = ps_real(i)%x(:) - 2._fp*dr*self%surface_normal(:)
            end if

        end do

    end subroutine ghost_timestep_setup_sweep

    subroutine ghost_state_update(self, ps, n, dt)
        class(ghost_state_updater_t), intent(in):: self
        integer, intent(in):: n
        class(base_particle_t), intent(inout):: ps(n)
        real(fp), optional, intent(in):: dt
        integer:: i
        real(fp):: projection(ndims)

        select type (ps_ghost => ps)
        class is (eos_ghost_particle_t)
            do i = 1, n
                projection(:) = dot_product(ps_ghost(i)%original%v(:), self%surface_normal(:))*self%surface_normal(:)
                ps_ghost(i)%v(:) = ps_ghost(i)%original%v(:) - 2._fp*projection(:)
                ps_ghost(i)%rho = ps_ghost(i)%original%rho
                ps_ghost(i)%mass = ps_ghost(i)%original%mass
                ps_ghost(i)%p = ps_ghost(i)%original%p
                ps_ghost(i)%c = ps_ghost(i)%original%c
            end do
        class default
            error stop "Expected self%particles to be eos_ghost_particle_t."
        end select

    end subroutine ghost_state_update

end module grasph_dambreak_ghost_morris_boundary_m

program main

    use grasph_dambreak_ghost_morris_boundary_m

    use grasph_particle_system_m, only: particle_system_t
    use weakly_compressible_particles_m, only: tait_eos_state_updater_t, eos_particle_t
    use grasph_system_interactions_m, only: system_interaction_t
    use grasph_time_integration_m, only: leap_frog_time_integration
    use grasph_kernels_m, only: cubic_bspline_kernel_t

    implicit none
    type(particle_system_t):: psys(4)
    type(system_interaction_t):: psys_interactions(4)
    type(cubic_bspline_kernel_t):: kernel
    type(fluid_sweeper_t):: sweeper
    type(boundary_update_sweeper_t):: boundary_sweeper
    type(xsph_shifter_t):: shifter
    type(tait_eos_state_updater_t):: state_updater
    type(eos_particle_t):: ps_template
    type(ghost_state_updater_t):: ghost_state_updater
    type(ghost_timestep_setuper_t):: ghost_timestep_setuper
    type(eos_ghost_particle_t):: ghost_ps_template
    integer:: i, j, k, nlayer, nfx, nfy

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    ! init fluid particles
    state_updater%rho_ref = rho0
    call psys(1)%init(n=2500, name="fluid", state_updater=state_updater, particle_template=ps_template)

    ! register variables for time-update
    call psys(1)%register_x%register(psys(1)%particles(1), "x", psys(1)%particles(1)%x, psys(1)%particles(1)%v)
    call psys(1)%register_v%register(psys(1)%particles(1), "v", psys(1)%particles(1)%v, psys(1)%particles(1)%dvxdt)
    call psys(1)%register_v%register(psys(1)%particles(1), "rho", psys(1)%particles(1)%rho, psys(1)%particles(1)%drhodt)

    ! register variables for io
    select type (p => psys(1)%particles)
    class is (eos_particle_t)
        call psys(1)%register_io%register_variable(p(1), "x", p(1)%x)
        call psys(1)%register_io%register_variable(p(1), "v", p(1)%v)
        call psys(1)%register_io%register_variable(p(1), "rho", p(1)%rho)
        call psys(1)%register_io%register_variable(p(1), "mass", p(1)%mass)
        call psys(1)%register_io%register_variable(p(1), "c", p(1)%c)
        call psys(1)%register_io%register_variable(p(1), "dvxdt", p(1)%dvxdt)
        call psys(1)%register_io%register_variable(p(1), "drhodt", p(1)%drhodt)
        call psys(1)%register_io%register_variable(p(1), "p", p(1)%p)
    class default
        error stop "Expected eos_particle_t for psys(1)%p."
    end select

    nfx = 25._fp/dx
    nfy = 25._fp/dx

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
    ! only need to initialize metadata and position as only position is used to calculate repulsive force
    nlayer = ceiling(kernel%cutoff/dx)

    call generate_boundary(psys(2), 25._fp, 40._fp)

    ghost_state_updater%surface_normal(:) = [1._fp, 0._fp]
    call psys(3)%init(n=2500, name="ghost_boundary_left", particle_template=ghost_ps_template, state_updater=ghost_state_updater)
    call psys(4)%init(n=2500, name="ghost_boundary_right", particle_template=ghost_ps_template, state_updater=ghost_state_updater)
    ! register variables for io
    select type (p => psys(3)%particles)
    class is (eos_ghost_particle_t)
        call psys(3)%register_io%register_variable(p(1), "x", p(1)%x)
        call psys(3)%register_io%register_variable(p(1), "v", p(1)%v)
        call psys(3)%register_io%register_variable(p(1), "rho", p(1)%rho)
        call psys(3)%register_io%register_variable(p(1), "mass", p(1)%mass)
        call psys(3)%register_io%register_variable(p(1), "c", p(1)%c)
        call psys(3)%register_io%register_variable(p(1), "dvxdt", p(1)%dvxdt)
        call psys(3)%register_io%register_variable(p(1), "drhodt", p(1)%drhodt)
        call psys(3)%register_io%register_variable(p(1), "p", p(1)%p)
    class default
        error stop "Expected eos_ghost_particle_t for psys(4)%p."
    end select

    select type (p => psys(4)%particles)
    class is (eos_ghost_particle_t)
        call psys(4)%register_io%register_variable(p(1), "x", p(1)%x)
        call psys(4)%register_io%register_variable(p(1), "v", p(1)%v)
        call psys(4)%register_io%register_variable(p(1), "rho", p(1)%rho)
        call psys(4)%register_io%register_variable(p(1), "mass", p(1)%mass)
        call psys(4)%register_io%register_variable(p(1), "c", p(1)%c)
        call psys(4)%register_io%register_variable(p(1), "dvxdt", p(1)%dvxdt)
        call psys(4)%register_io%register_variable(p(1), "drhodt", p(1)%drhodt)
        call psys(4)%register_io%register_variable(p(1), "p", p(1)%p)
    class default
        error stop "Expected eos_ghost_particle_t for psys(4)%p."
    end select

    sweeper%artvisc_alpha = 0.01_fp
    sweeper%artvisc_beta = 0._fp
    sweeper%g = g
    sweeper%h = 1.2_fp*dx
    sweeper%update_rhs = .false.

    shifter%epsilon = 0.5_fp
    shifter%update_rhs = .false.

    ! init interactions
    call psys_interactions(1)%init(30, psys(1), sweeper=sweeper, shifter=shifter)
    sweeper%initialize = .false. ! second sweeper doesn't need to zero acceleration arrays
    call psys_interactions(2)%init(30, psys(1), psys(2), prologue_sweeper=boundary_sweeper, sweeper=sweeper, shifter=shifter)
    ghost_timestep_setuper%cutoff = kernel%cutoff
    ghost_timestep_setuper%surface_normal(:) = [1._fp, 0._fp]
    ghost_timestep_setuper%point(:) = [0._fp, 0._fp]
    call psys_interactions(3)%init(30, psys(1), psys(3), timestep_setuper=ghost_timestep_setuper, sweeper=sweeper)
    ghost_timestep_setuper%surface_normal(:) = [-1._fp, 0._fp]
    ghost_timestep_setuper%point(:) = [25._fp, 0._fp]
    call psys_interactions(4)%init(30, psys(1), psys(4), timestep_setuper=ghost_timestep_setuper, sweeper=sweeper)

    ! start time-evolution with damping for setting up of initial conditions for fluid.
    call leap_frog_time_integration( &
        maxtimestep=50000, &
        print_step=1000, &
        save_step=1000, &
        psystems=psys, &
        interactions=psys_interactions, &
        CFL=0.05_fp, &
        kernel=kernel, &
        output_path="/home/edwardy/test", &
        output_prefix="damping", &
        output_comp_level=4, &
        damping_coef=390._fp &
        )

    ! re-initialize boundary conditions so that geometry matches the dambreak setup (without ramp).
    call generate_boundary(psys(2), 75._fp, 40._fp)

    ghost_timestep_setuper%cutoff = kernel%cutoff
    ghost_timestep_setuper%point(:) = [75._fp, 0._fp]
    call psys_interactions(4)%init(30, psys(1), psys(4), timestep_setuper=ghost_timestep_setuper, sweeper=sweeper)

    ! start time-evolution what dambreak setup and without damping.
    call leap_frog_time_integration( &
        maxtimestep=100000, &
        print_step=1000, &
        save_step=1000, &
        psystems=psys, &
        interactions=psys_interactions, &
        CFL=0.05_fp, &
        kernel=kernel, &
        output_path="/home/edwardy/test", &
        output_comp_level=4 &
        )

contains

    subroutine generate_boundary(psys_boundary, extx, exty)

        type(particle_system_t), intent(out):: psys_boundary
        real(fp), intent(in):: extx, exty
        integer:: nbx, nby, nvirt

        nbx = extx/dx
        nby = exty/dx

        nvirt = 2*nlayer*nbx + 4*nlayer*nlayer

        call psys_boundary%init(nvirt, name="boundary", state_updater=state_updater, particle_template=ps_template)

        select type (p => psys_boundary%particles)
        class is (eos_particle_t)

            call psys_boundary%register_io%register_variable(p(1), "x", p(1)%x)
            call psys_boundary%register_io%register_variable(p(1), "v", p(1)%v)
            call psys_boundary%register_io%register_variable(p(1), "rho", p(1)%rho)
            call psys_boundary%register_io%register_variable(p(1), "mass", p(1)%mass)
            call psys_boundary%register_io%register_variable(p(1), "c", p(1)%c)
            call psys_boundary%register_io%register_variable(p(1), "dvxdt", p(1)%dvxdt)
            call psys_boundary%register_io%register_variable(p(1), "drhodt", p(1)%drhodt)
            call psys_boundary%register_io%register_variable(p(1), "p", p(1)%p)
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
                psys(2)%particles(k)%x(2) = exty + (j + 0.5_fp)*dx
            end do
        end do
        do i = 1, k
            psys(2)%particles(i)%id = i
            psys(2)%particles(i)%type = -1
            psys(2)%particles(i)%rho = rho0
            psys(2)%particles(i)%mass = rho0*dx*dx
            psys(2)%particles(i)%c = 10._fp*sqrt(490.5_fp)
        end do

    end subroutine generate_boundary

end program main
