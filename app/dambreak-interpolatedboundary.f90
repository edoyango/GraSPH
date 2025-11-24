!> @file dambreak-interpolatedboundary.f90
!> @brief This program sets up and runs the classic dambreak experiment which consists of a 25m by 25m 2D block of water,
!>        which is released to move round a box of size 75m by 40m. This version of the experiment uses virtual particles to
!>        enforce fully-fixed boundaries for all the walls. These virtual particles use kernel interpolation to enforice this
!>        condition.
!> @author Edward Yang
!> @date 2025-09-23
program main

    use grasph_constants_m, only: fp
    use grasph_particle_system_m, only: particle_system_t, state_updater_container_t
    use weakly_compressible_interactions_m, only: fluid_sweeper_t, boundary_update_sweeper_t
    use weakly_compressible_particles_m, only: tait_eos_state_updater_t, eos_particles_t
    use grasph_system_interactions_m, only: system_interaction_t, default_sweeper_t, sweeper_container_t
    use grasph_time_integration_m, only: leap_frog_time_integration
    use grasph_kernels_m, only: cubic_bspline_kernel_t
    use grasph_particle_shifting_m, only: xsph_shifter_t

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.5_fp, g = -9.81_fp, rho0 = 1000._fp
    type(particle_system_t):: psys(2)
    type(system_interaction_t):: psys_interactions(2)
    type(cubic_bspline_kernel_t):: kernel
    type(fluid_sweeper_t):: sweeper
    type(boundary_update_sweeper_t):: boundary_sweeper
    type(xsph_shifter_t):: shifter
    type(tait_eos_state_updater_t):: state_updater
    type(eos_particles_t):: ps_template
    type(default_sweeper_t):: donothing_sweeper
    type(sweeper_container_t):: fluid_fluid_sweepers(2), fluid_boundary_sweepers(2)
    type(state_updater_container_t):: fluid_state_updaters(2)
    integer:: i, j, k, nlayer, nfx, nfy

    ! init kernel
    call kernel%init(1.2_fp*dx)

    ! init fluid particles
    state_updater%rho_ref = rho0 ! EOS only needs to know the reference density.
    allocate (fluid_state_updaters(1)%updater)
    allocate (fluid_state_updaters(2)%updater, source=state_updater)
    call psys(1)%init( &
        n=2500, &
        name="fluid", &
        state_updaters=fluid_state_updaters, &
        particle_template=ps_template &
        )

    ! register variables for time-update
    call psys(1)%register_x%register("x", psys(1)%particles%x, psys(1)%particles%v)
    call psys(1)%register_v%register("v", psys(1)%particles%v, psys(1)%particles%dvxdt)
    call psys(1)%register_v%register("rho", psys(1)%particles%rho, psys(1)%particles%drhodt)

    ! register variables for io
    select type (p => psys(1)%particles)
    class is (eos_particles_t)
        call psys(1)%register_io%register_variable("p", p%p)
    class default
        error stop "Expected eos_particles_t for psys(1)%p."
    end select

    nfx = nint(25._fp/dx)
    nfy = nint(25._fp/dx)

    do i = 0, nfx - 1
        do j = 0, nfy - 1
            k = j*nfx + i + 1
            psys(1)%particles%id(k) = k
            psys(1)%particles%type(k) = 1 ! not sure if type is needed anymore
            psys(1)%particles%x(1, k) = (i + 0.5_fp)*dx
            psys(1)%particles%x(2, k) = (j + 0.5_fp)*dx
            psys(1)%particles%rho(k) = rho0
            psys(1)%particles%mass(k) = rho0*dx*dx
            psys(1)%particles%c(k) = 10._fp*sqrt(490.5_fp) ! 10*sqrt(2gH)
            psys(1)%particles%v(:, k) = 0._fp
        end do
    end do

    ! init boundary particles.
    ! only need to initialize metadata and position as only position is used to calculate repulsive force.
    nlayer = ceiling(kernel%cutoff/dx)

    call generate_boundary(psys(2), 25._fp, 40._fp)

    ! declare params used for system interactions. This applies to interactions between both fluid and fluid, and fluid and
    ! boundary.
    sweeper%artvisc_alpha = 0.01_fp
    sweeper%artvisc_beta = 0._fp
    sweeper%g = g
    sweeper%h = 1.2_fp*dx
    sweeper%update_rhs = .false. ! rhs particles are boundary.

    ! declare XSPH shifter params.
    shifter%epsilon = 0.5_fp
    shifter%update_rhs = .false.

    ! init interactions
    allocate (fluid_fluid_sweepers(1)%sweeper, source=donothing_sweeper)
    allocate (fluid_fluid_sweepers(2)%sweeper, source=sweeper)
    call psys_interactions(1)%init( &
        npairs_per_particle=30, &
        psys_lhs=psys(1), &
        sweepers=fluid_fluid_sweepers, &
        shifter=shifter &
        )
    sweeper%initialise = .false. ! second sweeper doesn't need to zero acceleration arrays
    allocate (fluid_boundary_sweepers(1)%sweeper, source=boundary_sweeper)
    allocate (fluid_boundary_sweepers(2)%sweeper, source=sweeper)
    call psys_interactions(2)%init( &
        npairs_per_particle=30, &
        psys_lhs=psys(1), &
        psys_rhs=psys(2), &
        sweepers=fluid_boundary_sweepers, &
        shifter=shifter &
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
    call generate_boundary(psys(2), 75._fp, 40._fp)

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

    !> @brief Helper function to generate the walls used in the dambreak experiment. The corner of the walls will be at
    !>        (0, 0), (extx, exty). Also registers the relevant variables for IO.
    !> @param psys_boundary The particle system object to initialize with boundary particles.
    !> @param extx The x-extent of the boundary.
    !> @param exty The y-extent of the boundary.
    subroutine generate_boundary(psys_boundary, extx, exty)

        type(particle_system_t), intent(out):: psys_boundary
        real(fp), intent(in):: extx, exty
        integer:: nbx, nby, nvirt
        type(state_updater_container_t):: boundary_state_updaters(2)

        nbx = nint(extx/dx)
        nby = nint(exty/dx)

        nvirt = 2*nlayer*(nbx + nby) + 4*nlayer*nlayer

        allocate (boundary_state_updaters(1)%updater)
        allocate (boundary_state_updaters(2)%updater, source=state_updater)
        call psys_boundary%init(nvirt, name="boundary", state_updaters=boundary_state_updaters, particle_template=ps_template)

        select type (p => psys_boundary%particles)
        class is (eos_particles_t)
            call psys_boundary%register_io%register_variable("p", p%p)
        end select

        k = 0
        ! bottom layer and corners
        do i = -nlayer, nbx + nlayer - 1
            do j = 0, nlayer - 1
                k = k + 1
                psys(2)%particles%x(1, k) = (i + 0.5_fp)*dx
                psys(2)%particles%x(2, k) = -(j + 0.5_fp)*dx
            end do
        end do
        ! top layer and corners
        do i = -nlayer, nbx + nlayer - 1
            do j = 0, nlayer - 1
                k = k + 1
                psys(2)%particles%x(1, k) = (i + 0.5_fp)*dx
                psys(2)%particles%x(2, k) = exty + (j + 0.5_fp)*dx
            end do
        end do
        ! left wall
        do j = 0, nby - 1
            do i = 0, nlayer - 1
                k = k + 1
                psys(2)%particles%x(1, k) = -(i + 0.5_fp)*dx
                psys(2)%particles%x(2, k) = (j + 0.5_fp)*dx
            end do
        end do
        ! right wall
        do j = 0, nby - 1
            do i = 0, nlayer - 1
                k = k + 1
                psys(2)%particles%x(1, k) = extx + (i + 0.5_fp)*dx
                psys(2)%particles%x(2, k) = (j + 0.5_fp)*dx
            end do
        end do
        do i = 1, k
            psys(2)%particles%id(i) = i
            psys(2)%particles%type(i) = -1
            psys(2)%particles%rho(i) = rho0
            psys(2)%particles%mass(i) = rho0*dx*dx
            psys(2)%particles%c(i) = 10._fp*sqrt(490.5_fp)
        end do

    end subroutine generate_boundary

end program main
