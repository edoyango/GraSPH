module grasph_dambreak_ghost_morris_boundary_m

    use grasph_constants_m, only: fp, ndims
    use grasph_particle_m, only: base_particle_t
    use grasph_particle_system_m, only: particle_system_t, base_state_updater_t
    use weakly_compressible_particles_m, only: eos_particle_t, eos_ghost_particle_t, ghost_state_updater_t
    use weakly_compressible_interactions_m, only: fluid_sweeper_t, ghost_timestep_setuper_t, morris_boundary_sweeper_t
    use grasph_system_interactions_m, only: base_sweeper_t
    use grasph_pair_interactions_m, only: artificial_viscosity_monaghan1994, continuity_density, isotropic_pressure_force, &
                                          repulsive_force
    use grasph_pairs_m, only: particle_pairs_t
    use grasph_particle_shifting_m, only: xsph_shifter_t

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.5_fp, g = -9.81_fp, rho0 = 1000._fp

contains

end module grasph_dambreak_ghost_morris_boundary_m

program main

    use grasph_dambreak_ghost_morris_boundary_m

    use grasph_particle_system_m, only: particle_system_t
    use weakly_compressible_particles_m, only: tait_eos_state_updater_t, eos_particle_t
    use grasph_system_interactions_m, only: system_interaction_t
    use grasph_time_integration_m, only: leap_frog_time_integration
    use grasph_kernels_m, only: cubic_bspline_kernel_t

    implicit none
    type(particle_system_t):: psys(5)
    type(system_interaction_t):: psys_interactions(5)
    type(cubic_bspline_kernel_t):: kernel
    type(fluid_sweeper_t):: sweeper
    type(morris_boundary_sweeper_t):: boundary_sweeper
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
    call psys(1)%init(n=2500, name="fluid", state_updater_1=state_updater, particle_template=ps_template)

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

    call generate_boundary(psys(2), 25._fp, .true.)
    call generate_boundary(psys(3), 25._fp, .false.)

    ghost_state_updater%surface_normal(:) = [1._fp, 0._fp]
    call psys(4)%init(n=2500, name="ghost_boundary_left", particle_template=ghost_ps_template, state_updater_1=ghost_state_updater)
    call psys(5)%init(n=2500, name="ghost_boundary_right", particle_template=ghost_ps_template, state_updater_1=ghost_state_updater)
    ! register variables for io
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

    select type (p => psys(4)%particles)
    class is (eos_ghost_particle_t)
        call psys(5)%register_io%register_variable(p(1), "x", p(1)%x)
        call psys(5)%register_io%register_variable(p(1), "v", p(1)%v)
        call psys(5)%register_io%register_variable(p(1), "rho", p(1)%rho)
        call psys(5)%register_io%register_variable(p(1), "mass", p(1)%mass)
        call psys(5)%register_io%register_variable(p(1), "c", p(1)%c)
        call psys(5)%register_io%register_variable(p(1), "dvxdt", p(1)%dvxdt)
        call psys(5)%register_io%register_variable(p(1), "drhodt", p(1)%drhodt)
        call psys(5)%register_io%register_variable(p(1), "p", p(1)%p)
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

    boundary_sweeper%artvisc_alpha = 0.01_fp
    boundary_sweeper%artvisc_beta = 0._fp
    boundary_sweeper%h = 1.2_fp*dx

    ! init interactions
    call psys_interactions(1)%init(30, psys(1), sweeper=sweeper, shifter=shifter)
    sweeper%initialize = .false. ! second sweeper doesn't need to zero acceleration arrays
    boundary_sweeper%point(:) = [0._fp, 0._fp]
    boundary_sweeper%normal(:) = [0._fp, 1._fp]
    call psys_interactions(2)%init(30, psys(1), psys(2), sweeper=boundary_sweeper, shifter=shifter)
    boundary_sweeper%point(:) = [0._fp, 40._fp]
    boundary_sweeper%normal(:) = [0._fp, -1._fp]
    call psys_interactions(3)%init(30, psys(1), psys(3), sweeper=boundary_sweeper, shifter=shifter)
    ghost_timestep_setuper%cutoff = kernel%cutoff
    ghost_timestep_setuper%surface_normal(:) = [1._fp, 0._fp]
    ghost_timestep_setuper%point(:) = [0._fp, 0._fp]
    call psys_interactions(4)%init(30, psys(1), psys(4), timestep_setuper=ghost_timestep_setuper, sweeper=sweeper)
    ghost_timestep_setuper%surface_normal(:) = [-1._fp, 0._fp]
    ghost_timestep_setuper%point(:) = [25._fp, 0._fp]
    call psys_interactions(5)%init(30, psys(1), psys(5), timestep_setuper=ghost_timestep_setuper, sweeper=sweeper)

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
    call generate_boundary(psys(2), 75._fp, .true.)
    call generate_boundary(psys(3), 75._fp, .false.)

    ghost_timestep_setuper%cutoff = kernel%cutoff
    ghost_timestep_setuper%point(:) = [75._fp, 0._fp]
    call psys_interactions(5)%init(30, psys(1), psys(5), timestep_setuper=ghost_timestep_setuper, sweeper=sweeper)

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

    subroutine generate_boundary(psys_boundary, extx, bottom)

        type(particle_system_t), intent(out):: psys_boundary
        real(fp), intent(in):: extx
        logical, intent(in):: bottom
        integer:: nbx, nvirt
        type(base_particle_t):: p

        nbx = extx/dx

        nvirt = nlayer*nbx + 2*nlayer*nlayer

        k = 0
        if (bottom) then
            ! bottom layer and corners
            call psys_boundary%init(nvirt, name="bottom_boundary")
            do i = -nlayer, nbx + nlayer - 1
                do j = 0, nlayer - 1
                    k = k + 1
                    psys_boundary%particles(k)%x(1) = (i + 0.5_fp)*dx
                    psys_boundary%particles(k)%x(2) = -(j + 0.5_fp)*dx
                end do
            end do
        else
            ! top layer and corners
            call psys_boundary%init(nvirt, name="top_boundary")
            do i = -nlayer, nbx + nlayer - 1
                do j = 0, nlayer - 1
                    k = k + 1
                    psys_boundary%particles(k)%x(1) = (i + 0.5_fp)*dx
                    psys_boundary%particles(k)%x(2) = 40._fp + (j + 0.5_fp)*dx
                end do
            end do
        end if
        do i = 1, k
            psys_boundary%particles(i)%id = i
            psys_boundary%particles(i)%type = -1
        end do

        psys_boundary%to_print_summary = .false.
        call psys_boundary%register_io%register_variable(psys_boundary%particles(1), "x", psys_boundary%particles(1)%x)

    end subroutine generate_boundary

end program main
