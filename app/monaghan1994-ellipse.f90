!> @file monaghan1994.f90
!> @brief Module and program containing setup code to replicate the classic dambreak experiment from
!>        Monaghan (1994) (https://doi.org/10.1006/jcph.1994.1034). The only known difference is that the Leap-Frog time-integration
!>        is used here instead of the predictor-corrector scheme used in the paper.
!> @author Edward Yang
!> @date 2025/9/22

program main

    use grasph_constants, only: fp, pi
    use grasph_particles, only: base_particles
    use weakly_compressible_particles, only: tait_eos_state_updater, linear_eos_particle
    use weakly_compressible_interactions, only: fluid_sweeper
    use grasph_pair_sets, only: particle_interactions
    use grasph_time_integration, only: leap_frog_time_integration
    use grasph_kernels, only: grasph_cubic_bspline_kernel
    use grasph_particle_shifting, only: xsph_shifter

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.04_fp, g = 0._fp, rho0 = 1000._fp
    ! no. of particles in x, y direction in initial geometry of fluid
    integer, parameter:: nfx = 2._fp/dx, nfy = 2._fp/dx
    type(base_particles):: ps(1)
    type(particle_interactions):: pic(1)
    type(grasph_cubic_bspline_kernel):: kernel
    integer:: i, j, k
    type(fluid_sweeper):: sweeper
    type(xsph_shifter):: shifter
    type(tait_eos_state_updater):: state_updater
    type(linear_eos_particle):: ps_template
    real(fp):: x, y

    ! init fluid particles
    state_updater%rho_ref = rho0

    ! register variables for time-update
    call ps(1)%base_init(n=1976, name="fluid", state_updater=state_updater, ps_template=ps_template)
    call ps(1)%register_x%register(ps(1)%ps(1), "x", ps(1)%ps(1)%x, ps(1)%ps(1)%v)
    call ps(1)%register_v%register(ps(1)%ps(1), "v", ps(1)%ps(1)%v, ps(1)%ps(1)%dvxdt)
    call ps(1)%register_v%register(ps(1)%ps(1), "rho", ps(1)%ps(1)%rho, ps(1)%ps(1)%drhodt)

    ! register variables for io
    select type (p => ps(1)%ps)
    class is (linear_eos_particle)
        call ps(1)%register_io%register_variable(p(1), "x", p(1)%x)
        call ps(1)%register_io%register_variable(p(1), "v", p(1)%v)
        call ps(1)%register_io%register_variable(p(1), "rho", p(1)%rho)
        call ps(1)%register_io%register_variable(p(1), "mass", p(1)%mass)
        call ps(1)%register_io%register_variable(p(1), "c", p(1)%c)
        call ps(1)%register_io%register_variable(p(1), "dvxdt", p(1)%dvxdt)
        call ps(1)%register_io%register_variable(p(1), "drhodt", p(1)%drhodt)
        call ps(1)%register_io%register_variable(p(1), "p", p(1)%p)
    class default
        error stop "Expected linear_eos_particle for ps(1)%p."
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
                ps(1)%ps(k)%id = k
                ps(1)%ps(k)%type = 1
                ps(1)%ps(k)%x(1) = x
                ps(1)%ps(k)%x(2) = y
                ps(1)%ps(k)%c = 1400._fp
                ps(1)%ps(k)%rho = rho0
                ps(1)%ps(k)%mass = pi*rho0/1976
                ps(1)%ps(k)%v(1) = -100._fp*x
                ps(1)%ps(k)%v(2) = 100._fp*y
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
    call pic(1)%init(30, ps(1), sweeper=sweeper, shifter=shifter)

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    call leap_frog_time_integration(5000, 10, 10, ps, pic, 0.05_fp, kernel, "/home/edwardy/test", output_comp_level=4)

end program main
