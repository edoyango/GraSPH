program main

    use grasph_constants_m, only: fp, pi
    use grasph_particle_system_m, only: particle_system_t, state_updater_container_t
    use weakly_compressible_particles_m, only: eos_viscous_stress_particles_t, eos_viscous_stress_ghost_particles_t, &
                                               eos_viscous_stress_ghost_state_updater_t, dp_visco_elastic_state_updater_t
    use weakly_compressible_interactions_m, only: viscous_stress_fluid_sweeper_t, eos_viscous_stress_ghost_timestep_setuper_t, &
                                                  eos_viscous_stress_morris_boundary_sweeper_t, strain_rate_sweeper_t, &
                                                  strain_rate_morris_boundary_sweeper_t

    use grasph_system_interactions_m, only: system_interaction_t, sweeper_container_t
    use grasph_time_integration_m, only: leap_frog_time_integration
    use grasph_kernels_m, only: cubic_bspline_kernel_t

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.002_fp, g = -9.81_fp, rho0 = 1850._fp
    real(fp), parameter:: soil_maxx = 0.2_fp, soil_maxy = 0.1_fp, boundary_maxx = 0.8_fp
    type(particle_system_t):: psys(3)
    type(system_interaction_t):: psys_interactions(3)
    type(cubic_bspline_kernel_t):: kernel
    type(viscous_stress_fluid_sweeper_t):: sweeper
    type(eos_viscous_stress_morris_boundary_sweeper_t):: boundary_sweeper
    type(dp_visco_elastic_state_updater_t):: state_updater
    type(eos_viscous_stress_particles_t):: ps_template
    type(eos_viscous_stress_ghost_state_updater_t):: ghost_state_updater
    type(eos_viscous_stress_ghost_timestep_setuper_t):: ghost_timestep_setuper
    type(eos_viscous_stress_ghost_particles_t):: ghost_ps_template
    type(strain_rate_sweeper_t):: strain_rate_sweeper
    type(strain_rate_morris_boundary_sweeper_t):: strain_rate_boundary_sweeper
    type(sweeper_container_t):: soil_soil_sweepers(2), soil_ghost_sweepers(2), soil_virtual_sweepers(2)
    type(state_updater_container_t):: soil_state_updaters(2), ghost_state_updaters(2), virtual_state_updaters(2)
    integer:: i, j, k, nlayer, nvirt, nfx, nfy, nbx

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    ! init fluid particles
    state_updater%rho_ref = rho0
    state_updater%friction_angle = 19.8_fp*pi/180._fp
    allocate (soil_state_updaters(1)%updater)
    allocate (soil_state_updaters(2)%updater, source=state_updater)
    call psys(1)%init(n=5000, name="soil", state_updaters=soil_state_updaters, particle_template=ps_template)

    ! register variables for time-update
    call psys(1)%register_x%register("x", psys(1)%particles%x, psys(1)%particles%v)
    call psys(1)%register_v%register("v", psys(1)%particles%v, psys(1)%particles%dvxdt)
    call psys(1)%register_v%register("rho", psys(1)%particles%rho, psys(1)%particles%drhodt)

    ! register variables for io
    select type (p => psys(1)%particles)
    class is (eos_viscous_stress_particles_t)
        call psys(1)%register_io%register_variable("p", p%p)
        call psys(1)%register_io%register_variable("stress", p%stress)
    class default
        error stop "Expected eos_viscous_stress_particles_t for psys(1)%p."
    end select

    nfx = nint(soil_maxx/dx)
    nfy = nint(soil_maxy/dx)

    do i = 0, nfx - 1
        do j = 0, nfy - 1
            k = j*nfx + i + 1
            psys(1)%particles%id(k) = k
            psys(1)%particles%type(k) = 1 ! not sure if type is needed anymore
            psys(1)%particles%x(1, k) = (i + 0.5_fp)*dx
            psys(1)%particles%x(2, k) = (j + 0.5_fp)*dx
            psys(1)%particles%rho(k) = rho0
            psys(1)%particles%mass(k) = rho0*dx*dx
            psys(1)%particles%c(k) = 20._fp
            psys(1)%particles%v(:, k) = 0._fp
        end do
    end do

    ! init boundary particles
    ! only need to initialize metadata and position as only position is used to calculate repulsive force
    nlayer = ceiling(kernel%cutoff/dx)

    nbx = nint(boundary_maxx/dx)

    nvirt = nlayer*nbx + nlayer*nlayer
    ! bottom layer and corners
    call psys(2)%init(nvirt, name="bottom_boundary", particle_template=ps_template)
    ! deregister most variables from IO for boundary
    call psys(2)%register_io%deregister("v")
    call psys(2)%register_io%deregister("rho")
    call psys(2)%register_io%deregister("mass")
    call psys(2)%register_io%deregister("c")
    call psys(2)%register_io%deregister("dvxdt")
    call psys(2)%register_io%deregister("drhodt")

    k = 0
    do i = -nlayer, nbx - 1
        do j = 0, nlayer - 1
            k = k + 1
            psys(2)%particles%x(1, k) = (i + 0.5_fp)*dx
            psys(2)%particles%x(2, k) = -(j + 0.5_fp)*dx
            psys(2)%particles%id(k) = i
            psys(2)%particles%type(k) = -1
        end do
    end do

    psys(2)%to_print_summary = .false.

    ghost_state_updater%surface_normal(:) = [1._fp, 0._fp]
    allocate (ghost_state_updaters(1)%updater, source=ghost_state_updater)
    ! second ghost updater because real particles' stress is updated in second stage.
    allocate (ghost_state_updaters(2)%updater, source=ghost_state_updater)
    call psys(3)%init( &
        n=5000, &
        name="ghost_boundary_left", &
        particle_template=ghost_ps_template, &
        state_updaters=ghost_state_updaters &
        )

    ! register variables for io
    select type (p => psys(3)%particles)
    class is (eos_viscous_stress_ghost_particles_t)
        call psys(3)%register_io%register_variable("p", p%p)
        call psys(3)%register_io%register_variable("stress", p%stress)
    class default
        error stop "Expected eos_ghost_particle_t for psys(4)%p."
    end select
    ! deregister rate-of-change arrays from IO
    call psys(3)%register_io%deregister("dvxdt")
    call psys(3)%register_io%deregister("drhodt")

    strain_rate_sweeper%update_rhs = .false.
    strain_rate_boundary_sweeper%point(:) = 0._fp
    strain_rate_boundary_sweeper%normal(:) = [0._fp, 1._fp]

    sweeper%artvisc_alpha = 0.1_fp
    sweeper%artvisc_beta = 0.1_fp
    sweeper%g = g
    sweeper%h = 1.2_fp*dx
    sweeper%update_rhs = .false.

    boundary_sweeper%artvisc_alpha = 0.1_fp
    boundary_sweeper%artvisc_beta = 0.1_fp
    boundary_sweeper%h = 1.2_fp*dx

    ! init interactions
    allocate (soil_soil_sweepers(1)%sweeper, source=strain_rate_sweeper)
    allocate (soil_soil_sweepers(2)%sweeper, source=sweeper)
    call psys_interactions(1)%init( &
        30, &
        psys(1), &
        sweepers=soil_soil_sweepers &
        )
    strain_rate_sweeper%initialise = .false.
    sweeper%initialise = .false. ! second sweeper doesn't need to zero acceleration arrays
    boundary_sweeper%initialise = .false.
    boundary_sweeper%point(:) = [0._fp, 0._fp]
    boundary_sweeper%normal(:) = [0._fp, 1._fp]
    allocate (soil_virtual_sweepers(1)%sweeper, source=strain_rate_boundary_sweeper)
    allocate (soil_virtual_sweepers(2)%sweeper, source=boundary_sweeper)
    call psys_interactions(2)%init( &
        30, &
        psys(1), &
        psys(2), &
        sweepers=soil_virtual_sweepers &
        )
    ghost_timestep_setuper%cutoff = kernel%cutoff
    ghost_timestep_setuper%surface_normal(:) = [1._fp, 0._fp]
    ghost_timestep_setuper%point(:) = [0._fp, 0._fp]
    allocate (soil_ghost_sweepers(1)%sweeper, source=strain_rate_sweeper)
    allocate (soil_ghost_sweepers(2)%sweeper, source=sweeper)
    call psys_interactions(3)%init( &
        30, &
        psys(1), &
        psys(3), &
        timestep_setuper=ghost_timestep_setuper, &
        sweepers=soil_ghost_sweepers &
        )

    ! start time-evolution what dambreak setup and without damping.
    call leap_frog_time_integration( &
        maxtimestep=50000, &
        print_step=1000, &
        save_step=1000, &
        psystems=psys, &
        interactions=psys_interactions, &
        CFL=0.1_fp, &
        kernel=kernel, &
        output_path="/home/edwardy/test", &
        output_comp_level=4 &
        )

end program main
