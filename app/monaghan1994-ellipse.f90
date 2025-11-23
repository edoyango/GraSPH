!> @file monaghan1994.f90
!> @brief Module and program containing setup code to replicate the classic dambreak experiment from
!>        Monaghan (1994) (https://doi.org/10.1006/jcph.1994.1034). The only known difference is that the Leap-Frog time-integration
!>        is used here instead of the predictor-corrector scheme used in the paper.
!> @author Edward Yang
!> @date 2025-9-22

program main

    use grasph_constants_m, only: fp, pi
    use grasph_particle_system_m, only: particle_system_t, state_updater_container_t
    use weakly_compressible_particles_m, only: tait_eos_state_updater_t, eos_particles_t
    use weakly_compressible_interactions_m, only: fluid_sweeper_t
    use grasph_system_interactions_m, only: system_interaction_t, sweeper_container_t, default_sweeper_t
    use grasph_time_integration_m, only: leap_frog_time_integration
    use grasph_kernels_m, only: cubic_bspline_kernel_t
    use grasph_particle_shifting_m, only: xsph_shifter_t

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.04_fp, g = 0._fp, rho0 = 1000._fp
    ! no. of particles in x, y direction in initial geometry of fluid
    integer, parameter:: nfx = 2._fp/dx, nfy = 2._fp/dx
    type(particle_system_t):: psys(1)
    type(system_interaction_t):: psys_interactions(1)
    type(cubic_bspline_kernel_t):: kernel
    integer:: i, j, k
    type(fluid_sweeper_t):: sweeper
    type(xsph_shifter_t):: shifter
    type(tait_eos_state_updater_t):: state_updater
    type(eos_particles_t):: ps_template
    type(sweeper_container_t):: sweepers(1)
    type(default_sweeper_t):: donothing_sweeper
    type(state_updater_container_t):: state_updaters(1)
    real(fp):: x, y

    ! init fluid particles
    state_updater%rho_ref = rho0

    ! register variables for time-update
    allocate (state_updaters(1)%updater, source=state_updater)
    call psys(1)%init(n=1976, name="fluid", state_updaters=state_updaters, particle_template=ps_template)
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

    ! initialize geometry
    ! generate particles in a grid and save only the ones within the 1 radius circle
    k = 0
    do i = 0, nfx - 1
        do j = 0, nfy - 1
            x = -1._fp + (i + 0.5_fp)*dx
            y = -1._fp + (j + 0.5_fp)*dx
            if (x*x + y*y < 1._fp) then
                k = k + 1
                psys(1)%particles%id(k) = k
                psys(1)%particles%type(k) = 1
                psys(1)%particles%x(1, k) = x
                psys(1)%particles%x(2, k) = y
                psys(1)%particles%c(k) = 1400._fp
                psys(1)%particles%rho(k) = rho0
                psys(1)%particles%mass(k) = pi*rho0/1976
                psys(1)%particles%v(1, k) = -100._fp*x
                psys(1)%particles%v(2, k) = 100._fp*y
            end if
        end do
    end do

    ! init interactions
    sweeper%artvisc_alpha = 0.01_fp
    sweeper%artvisc_beta = 0._fp
    sweeper%h = 1.2_fp*dx
    sweeper%g = g
    shifter%epsilon = 0.5_fp
    shifter%update_rhs = .true.
    allocate (sweepers(1)%sweeper, source=sweeper)
    call psys_interactions(1)%init(30, psys(1), sweepers=sweepers, shifter=shifter)

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    call leap_frog_time_integration( &
        maxtimestep=5000, &
        print_step=10, &
        save_step=10, &
        psystems=psys, &
        interactions=psys_interactions, &
        CFL=0.05_fp, &
        kernel=kernel, &
        output_path="/home/edwardy/test", &
        output_comp_level=4 &
        )

end program main
