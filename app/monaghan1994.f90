!> @file monaghan1994.f90
!> @brief Module and program containing setup code to replicate the classic dambreak experiment from
!>        Monaghan (1994) (https://doi.org/10.1006/jcph.1994.1034). The only known difference is that the Leap-Frog time-integration
!>        is used here instead of the predictor-corrector scheme used in the paper.
!> @author Edward Yang
!> @date 2025-9-22
program main

    use grasph_constants_m, only: fp
    use grasph_particle_system_m, only: particle_system_t, state_updater_container_t
    use weakly_compressible_particles_m, only: eos_particle_t, tait_eos_state_updater_t
    use weakly_compressible_interactions_m, only: fluid_sweeper_t, fluid_boundary_sweeper_monaghan1994_t
    use grasph_system_interactions_m, only: system_interaction_t, sweeper_container_t, default_sweeper_t
    use grasph_time_integration_m, only: leap_frog_time_integration
    use grasph_kernels_m, only: cubic_bspline_kernel_t
    use grasph_particle_shifting_m, only: xsph_shifter_t

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.5_fp, g = -9.81_fp, rho0 = 1000._fp
    type(particle_system_t):: psys(2)
    type(system_interaction_t):: psys_interactions(2)
    type(cubic_bspline_kernel_t):: kernel
    integer:: i, j, k
    type(fluid_sweeper_t):: self_sweeper
    type(fluid_boundary_sweeper_monaghan1994_t):: boundary_sweeper
    type(xsph_shifter_t):: shifter
    type(eos_particle_t):: ps_template
    type(tait_eos_state_updater_t):: state_updater
    type(sweeper_container_t):: fluid_fluid_sweepers(2), fluid_boundary_sweepers(2)
    type(default_sweeper_t):: donothing_sweeper
    type(state_updater_container_t):: fluid_state_updaters(2)
    integer:: nfx, nfy

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

    ! number of fluid particles in the x/y direction
    nfx = nint(25._fp/dx)
    nfy = nint(25._fp/dx)

    do i = 0, nfx - 1
        do j = 0, nfy - 1
            k = j*nfx + i + 1
            psys(1)%particles(k)%id = k
            psys(1)%particles(k)%type = 1 ! not sure if type is needed anymore
            psys(1)%particles(k)%x(1) = (i + 0.5_fp)*dx
            psys(1)%particles(k)%x(2) = (j + 0.5_fp)*dx
            psys(1)%particles(k)%c = 10._fp*sqrt(2._fp*abs(g)*25._fp) ! 10*sqrt(2gH) eqn 3.3
            psys(1)%particles(k)%rho = rho0 ! in the paper eqn 5.1 is used to initialize density, but doesn't seem to improve results.
            psys(1)%particles(k)%mass = rho0*dx*dx
            psys(1)%particles(k)%v(:) = 0._fp
        end do
    end do

    ! first generate boundary for setting up initial conditions where
    ! fluid is confined in a box and damping is applied.
    call generate_boundary(psys(2), 25._fp, 40._fp)! init interactions

    ! setup interactions between the fluid-fluid and fluid-boundary systems.
    self_sweeper%artvisc_alpha = 0.01_fp
    self_sweeper%artvisc_beta = 0._fp
    self_sweeper%h = 1.2_fp*dx
    self_sweeper%g = g
    boundary_sweeper%artvisc_alpha = 0.01_fp
    boundary_sweeper%artvisc_beta = 0._fp
    boundary_sweeper%h = 1.2_fp*dx
    boundary_sweeper%g = g
    boundary_sweeper%cutoff = dx ! disance at which lennard-jones repulsive force applies
    ! declare XSPH shifter params.
    shifter%epsilon = 0.5_fp
    shifter%update_rhs = .true. ! ensure that rhs particles of fluid-fluid interaction are updated.
    allocate (fluid_fluid_sweepers(1)%sweeper, source=donothing_sweeper)
    allocate (fluid_fluid_sweepers(2)%sweeper, source=self_sweeper)
    call psys_interactions(1)%init( &
        npairs_per_particle=30, &
        psys_lhs=psys(1), &
        sweepers=fluid_fluid_sweepers, &
        shifter=shifter &
        )
    shifter%update_rhs = .false. ! ensure that rhs particles of fluid-boundary interaction aren't updated.
    allocate (fluid_boundary_sweepers(1)%sweeper, source=donothing_sweeper)
    allocate (fluid_boundary_sweepers(2)%sweeper, source=boundary_sweeper)
    call psys_interactions(2)%init( &
        npairs_per_particle=30, &
        psys_lhs=psys(1), &
        psys_rhs=psys(2), &
        sweepers=fluid_boundary_sweepers, &
        shifter=shifter &
        )

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    ! start time-evolution with damping for setting up of initial conditions for fluid.
    call leap_frog_time_integration( &
        maxtimestep=6000, &
        print_step=1000, &
        save_step=1000, &
        psystems=psys, &
        interactions=psys_interactions, &
        CFL=0.05_fp, &
        kernel=kernel, &
        output_path="/home/edwardy/test", &
        output_prefix="damping", &
        output_comp_level=4, &
        damping_coef=390._fp & ! this should be enough to ensure damping_coef*dt ~= 0.05 (see comment at end of section 5 of paper).
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
        integer:: nbx, nby
        type(state_updater_container_t):: boundary_state_updaters(2)

        nbx = nint(extx/dx)
        nby = nint(exty/dx)

        allocate (boundary_state_updaters(1)%updater)
        allocate (boundary_state_updaters(2)%updater)
        call psys_boundary%init(2*(nbx + nby) + 4, name="boundary", state_updaters=boundary_state_updaters)
        psys_boundary%to_print_summary = .false.
        call psys_boundary%register_io%register_variable(psys_boundary%particles(1), "x", psys_boundary%particles(1)%x)
        call psys_boundary%register_io%register_variable(psys_boundary%particles(1), "v", psys_boundary%particles(1)%v)
        call psys_boundary%register_io%register_variable(psys_boundary%particles(1), "rho", psys_boundary%particles(1)%rho)
        call psys_boundary%register_io%register_variable(psys_boundary%particles(1), "mass", psys_boundary%particles(1)%mass)
        call psys_boundary%register_io%register_variable(psys_boundary%particles(1), "c", psys_boundary%particles(1)%c)
        call psys_boundary%register_io%register_variable(psys_boundary%particles(1), "dvxdt", psys_boundary%particles(1)%dvxdt)
        call psys_boundary%register_io%register_variable(psys_boundary%particles(1), "drhodt", psys_boundary%particles(1)%drhodt)

        k = 0
        ! bottom layer and corners
        do i = -1, nbx
            k = k + 1
            psys_boundary%particles(k)%x(1) = (i + 0.5_fp)*dx
            psys_boundary%particles(k)%x(2) = -0.5_fp*dx
        end do
        ! top layer and corners
        do i = -1, nbx
            k = k + 1
            psys_boundary%particles(k)%x(1) = (i + 0.5_fp)*dx
            psys_boundary%particles(k)%x(2) = exty + 0.5_fp*dx
        end do
        ! left wall
        do j = 0, nby - 1
            k = k + 1
            psys_boundary%particles(k)%x(1) = -0.5_fp*dx
            psys_boundary%particles(k)%x(2) = (j + 0.5_fp)*dx
        end do
        ! right wall
        do j = 0, nby - 1
            k = k + 1
            psys_boundary%particles(k)%x(1) = extx + 0.5_fp*dx
            psys_boundary%particles(k)%x(2) = (j + 0.5_fp)*dx
        end do
        do i = 1, k
            psys_boundary%particles(i)%id = i
            psys_boundary%particles(i)%type = -1
            psys_boundary%particles(i)%rho = rho0
            psys_boundary%particles(i)%mass = rho0*dx*dx
            psys_boundary%particles(i)%v(:) = 0._fp
            psys_boundary%particles(i)%c = 10._fp*sqrt(2._fp*abs(g)*25._fp) ! 10*max_speed
        end do

    end subroutine generate_boundary

end program main
