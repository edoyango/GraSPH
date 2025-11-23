!> @file dambreak-ghost-morris-boundary.f90
!> @brief This program sets up and runs the classic dambreak experiment which consists of a 25m by 25m 2D block of water,
!>        which is released to move round a box of size 75m by 40m. This version of the experiment uses ghost boundary particles to
!>        enforce the free-slip condition at the vertical walls, and the morris boundary description to enforce the fully-fixed
!>        condition at the horizontal walls.
!> @author Edward Yang
!> @date 2025-10-11
program main

    use grasph_constants_m, only: fp
    use grasph_particle_system_m, only: particle_system_t, base_state_updater_t, state_updater_container_t
    use weakly_compressible_particles_m, only: tait_eos_state_updater_t, ghost_state_updater_t
    use grasph_system_interactions_m, only: system_interaction_t, default_sweeper_t, sweeper_container_t
    use grasph_time_integration_m, only: leap_frog_time_integration
    use grasph_kernels_m, only: cubic_bspline_kernel_t
    use grasph_particle_shifting_m, only: xsph_shifter_t
    use weakly_compressible_interactions_m, only: fluid_sweeper_t, ghost_timestep_setuper_t, morris_boundary_sweeper_t

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.5_fp, g = -9.81_fp, rho0 = 1000._fp
    type(particle_system_t):: psys(5)
    type(system_interaction_t):: psys_interactions(5)
    type(cubic_bspline_kernel_t):: kernel
    type(fluid_sweeper_t):: sweeper
    type(morris_boundary_sweeper_t):: boundary_sweeper
    type(xsph_shifter_t):: shifter
    type(tait_eos_state_updater_t):: state_updater
    type(ghost_state_updater_t):: ghost_state_updater
    type(ghost_timestep_setuper_t):: ghost_timestep_setuper
    type(default_sweeper_t):: donothing_sweeper
    type(sweeper_container_t):: fluid_fluid_sweepers(1), fluid_left_wall_sweepers(1), fluid_right_wall_sweepers(1), &
                                fluid_bot_wall_sweepers(1), fluid_top_wall_sweepers(1)
    type(state_updater_container_t):: fluid_state_updaters(1), left_wall_state_updaters(1), right_wall_state_updaters(1)
    integer:: i, j, k, nlayer, nfx, nfy

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    ! init fluid particles
    state_updater%rho_ref = rho0 ! EOS only needs to know the reference density.
    allocate (fluid_state_updaters(1)%updater, source=state_updater)
    call psys(1)%init( &
        n=2500, &
        name="fluid", &
        state_updaters=fluid_state_updaters &
        )

    ! register variables for time-update
    call psys(1)%register_x%register(psys(1)%particles(1), "x", psys(1)%particles(1)%x, psys(1)%particles(1)%v)
    call psys(1)%register_v%register(psys(1)%particles(1), "v", psys(1)%particles(1)%v, psys(1)%particles(1)%dvxdt)
    call psys(1)%register_v%register(psys(1)%particles(1), "rho", psys(1)%particles(1)%rho, psys(1)%particles(1)%drhodt)

    nfx = nint(25._fp/dx)
    nfy = nint(25._fp/dx)

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
    ! only need to initialize metadata and position as only position is used to calculate repulsive force.
    nlayer = ceiling(kernel%cutoff/dx)

    call generate_boundary(psys(2), 25._fp, .true.)
    call generate_boundary(psys(3), 25._fp, .false.)

    ! initialize ghost boundaries
    ghost_state_updater%surface_normal(:) = [1._fp, 0._fp]
    allocate (left_wall_state_updaters(1)%updater, source=ghost_state_updater)
    call psys(4)%init( &
        n=2500, & ! allocate 2500 particles of space. Realistically, only 240 is neeeded.
        name="ghost_boundary_left", &
        state_updaters=left_wall_state_updaters &
        )
    ghost_state_updater%surface_normal(:) = [-1._fp, 0._fp] ! point normal leftward.
    allocate (right_wall_state_updaters(1)%updater, source=ghost_state_updater)
    call psys(5)%init( &
        n=2500, &
        name="ghost_boundary_right", &
        state_updaters=right_wall_state_updaters &
        )

    ! declare params used for fluid-fluid system interactions.
    sweeper%artvisc_alpha = 0.01_fp
    sweeper%artvisc_beta = 0._fp
    sweeper%g = g
    sweeper%h = 1.2_fp*dx
    sweeper%update_rhs = .false. ! rhs particles are boundary.

    ! use same params for boundary sweeper.
    boundary_sweeper%artvisc_alpha = 0.01_fp
    boundary_sweeper%artvisc_beta = 0._fp
    boundary_sweeper%h = 1.2_fp*dx
    ! update_rhs isn't used in this sweeper.

    ! declare XSPH shifter params.
    shifter%epsilon = 0.5_fp
    shifter%update_rhs = .false. ! rhs particles are boundary.

    ! init interactions.
    ! fluid self interaction.
    allocate (fluid_fluid_sweepers(1)%sweeper, source=sweeper)
    call psys_interactions(1)%init( &
        npairs_per_particle=30, & ! number of predicted interactions.
        psys_lhs=psys(1), & ! fluid particle system interacting with itself.
        sweepers=fluid_fluid_sweepers, & ! sweeper to describe interaction.
        shifter=shifter & ! shifter to perform XSPH shifting.
        )

    ! interaction between fluid and bottom wall (morris boundary).
    sweeper%initialize = .false. ! second sweeper doesn't need to zero acceleration arrays.
    boundary_sweeper%point(:) = [0._fp, 0._fp]
    boundary_sweeper%normal(:) = [0._fp, 1._fp] ! normal pointing upward.
    allocate (fluid_bot_wall_sweepers(1)%sweeper, source=boundary_sweeper)
    call psys_interactions(2)%init( &
        npairs_per_particle=30, &
        psys_lhs=psys(1), & ! fluid particle system.
        psys_rhs=psys(2), & ! boundary that fluid is interacting with.
        sweepers=fluid_bot_wall_sweepers, &
        shifter=shifter &
        )

    ! interaction between fluid and top wall (morris boundary).
    boundary_sweeper%point(:) = [0._fp, 40._fp]
    boundary_sweeper%normal(:) = [0._fp, -1._fp] ! normal pointing downward.
    allocate (fluid_top_wall_sweepers(1)%sweeper, source=boundary_sweeper)
    call psys_interactions(3)%init( &
        npairs_per_particle=30, &
        psys_lhs=psys(1), &
        psys_rhs=psys(3), &
        sweepers=fluid_top_wall_sweepers, &
        shifter=shifter &
        )

    ! interaction between fluid and left wall (ghost boundary).
    ghost_timestep_setuper%cutoff = kernel%cutoff
    ghost_timestep_setuper%surface_normal(:) = [1._fp, 0._fp] ! normal pointing rightward.
    ghost_timestep_setuper%point(:) = [0._fp, 0._fp]
    allocate (fluid_left_wall_sweepers(1)%sweeper, source=sweeper)
    call psys_interactions(4)%init( &
        npairs_per_particle=30, &
        psys_lhs=psys(1), &
        psys_rhs=psys(4), &
        timestep_setuper=ghost_timestep_setuper, &
        sweepers=fluid_left_wall_sweepers &
        )

    ! interaction between fluid and right wall (ghost boundary).
    ! setuper cutoff already set.
    ghost_timestep_setuper%surface_normal(:) = [-1._fp, 0._fp] ! normal pointing leftward.
    ghost_timestep_setuper%point(:) = [25._fp, 0._fp] ! right wall passes through (25, 0) for density initialization step.
    allocate (fluid_right_wall_sweepers(1)%sweeper, source=sweeper)
    call psys_interactions(5)%init( &
        npairs_per_particle=30, &
        psys_lhs=psys(1), &
        psys_rhs=psys(5), &
        timestep_setuper=ghost_timestep_setuper, &
        sweepers=fluid_right_wall_sweepers &
        )

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

    ! re-initialize interaction describing fluid and right wall since the setuper requires changing.
    ghost_timestep_setuper%cutoff = kernel%cutoff
    ghost_timestep_setuper%point(:) = [75._fp, 0._fp] ! update wall to pass through (75, 0) instead of (25, 0)
    call psys_interactions(5)%init( &
        npairs_per_particle=30, &
        psys_lhs=psys(1), &
        psys_rhs=psys(5), &
        timestep_setuper=ghost_timestep_setuper, &
        sweepers=fluid_right_wall_sweepers &
        )

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

    !> @brief Helper subroutine to generate upper and lower boundaries used in this simulation. The boundary is generated between
    !>        x = 0 and x = xext. Also registers the relevant variables for IO.
    !> @param psys_boundary The particle system object to initialize with boundary particles.
    !> @param extx The x-extent of the bottom boundary.
    !> @param bottom Whether the boundary being generated should be the bottom (y = 0) boundary, or the top (y = 40) boundary.
    subroutine generate_boundary(psys_boundary, extx, bottom)

        type(particle_system_t), intent(out):: psys_boundary
        real(fp), intent(in):: extx
        logical, intent(in):: bottom
        integer:: nbx, nvirt

        nbx = nint(extx/dx)

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

        ! morris boundary particles don't have persistent properties, so besides position, other data isn't needed.
        psys_boundary%to_print_summary = .false.

        ! deregister most variables from IO for boundary
        call psys_boundary%register_io%deregister("v")
        call psys_boundary%register_io%deregister("rho")
        call psys_boundary%register_io%deregister("mass")
        call psys_boundary%register_io%deregister("c")
        call psys_boundary%register_io%deregister("dvxdt")
        call psys_boundary%register_io%deregister("drhodt")

    end subroutine generate_boundary

end program main
