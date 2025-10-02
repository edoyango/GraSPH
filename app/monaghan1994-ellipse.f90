!> @file monaghan1994.f90
!> @brief Module and program containing setup code to replicate the classic dambreak experiment from
!>        Monaghan (1994) (https://doi.org/10.1006/jcph.1994.1034). The only known difference is that the Leap-Frog time-integration
!>        is used here instead of the predictor-corrector scheme used in the paper.
!> @author Edward Yang
!> @date 2025/9/22

program main

    use grasph_constants, only: fp, pi
    use grasph_particle_system_m, only: particle_system_t
    use weakly_compressible_particles, only: tait_eos_state_updater, eos_particle
    use weakly_compressible_interactions, only: fluid_sweeper
    use grasph_system_interactions_m, only: system_interaction_t
    use grasph_time_integration, only: leap_frog_time_integration
    use grasph_kernels_m, only: cubic_bspline_kernel_t
    use grasph_particle_shifting, only: xsph_shifter

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.04_fp, g = 0._fp, rho0 = 1000._fp
    ! no. of particles in x, y direction in initial geometry of fluid
    integer, parameter:: nfx = 2._fp/dx, nfy = 2._fp/dx
    type(particle_system_t):: psys(1)
    type(system_interaction_t):: psys_interactions(1)
    type(cubic_bspline_kernel_t):: kernel
    integer:: i, j, k
    type(fluid_sweeper):: sweeper
    type(xsph_shifter):: shifter
    type(tait_eos_state_updater):: state_updater
    type(eos_particle):: ps_template
    real(fp):: x, y

    ! init fluid particles
    state_updater%rho_ref = rho0

    ! register variables for time-update
    call psys(1)%base_init(n=1976, name="fluid", state_updater=state_updater, particle_template=ps_template)
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

    ! initialize geometry
    ! generate particles in a grid and save only the ones within the 1 radius circle
    k = 0
    do i = 0, nfx - 1
        do j = 0, nfy - 1
            x = -1._fp + (i + 0.5_fp)*dx
            y = -1._fp + (j + 0.5_fp)*dx
            if (x*x + y*y < 1._fp) then
                k = k + 1
                psys(1)%particles(k)%id = k
                psys(1)%particles(k)%type = 1
                psys(1)%particles(k)%x(1) = x
                psys(1)%particles(k)%x(2) = y
                psys(1)%particles(k)%c = 1400._fp
                psys(1)%particles(k)%rho = rho0
                psys(1)%particles(k)%mass = pi*rho0/1976
                psys(1)%particles(k)%v(1) = -100._fp*x
                psys(1)%particles(k)%v(2) = 100._fp*y
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
    call psys_interactions(1)%init(30, psys(1), sweeper=sweeper, shifter=shifter)

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
