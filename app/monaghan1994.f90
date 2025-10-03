!> @file monaghan1994.f90
!> @brief Module and program containing setup code to replicate the classic dambreak experiment from
!>        Monaghan (1994) (https://doi.org/10.1006/jcph.1994.1034). The only known difference is that the Leap-Frog time-integration
!>        is used here instead of the predictor-corrector scheme used in the paper.
!> @author Edward Yang
!> @date 2025/9/22

module grasph_monaghan1994

    use grasph_constants_m, only: fp
    use grasph_particle_system_m, only: particle_system_t
    use weakly_compressible_particles_m, only: eos_particle_t, tait_eos_state_updater_t
    use grasph_pairs_m, only: particle_pairs_t
    use grasph_system_interactions_m, only: base_sweeper_t
    use weakly_compressible_interactions_m, only: fluid_sweeper_t
    use grasph_pair_interactions, only: artificial_viscosity_monaghan1994, continuity_density, repulsive_force
    use grasph_particle_shifting_m, only: xsph_shifter_t

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.5_fp, g = -9.81_fp, rho0 = 1000._fp
    ! no. of particles in x, y direction in initial geometry of fluid
    integer, parameter:: nfx = 25._fp/dx, nfy = 25._fp/dx
    ! no. of particles in x, y direction for boundary
    integer, parameter:: nbx = 75._fp/dx, nby = 40._fp/dx

    type, extends(base_sweeper_t):: fluid_boundary_sweeper_t
        !> @brief Acceleration due to gravity (m/s)
        real(fp):: g = -9.81_fp
        !> @brief Alpha coefficient for artificial viscosity.
        real(fp):: artvisc_alpha = 0.1_fp
        !> @brief Beta coefficient for artificial viscosity.
        real(fp):: artvisc_beta = 0.1_fp
        !> @brief Smoothing length to use fo artificial viscosity.
        real(fp):: h = 0._fp
    contains
        procedure:: sweep => fluid_boundary_sweep_new
    end type fluid_boundary_sweeper_t

contains

    subroutine fluid_boundary_sweep_new(self, pairs, psys_lhs, psys_rhs)
        class(fluid_boundary_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        integer:: i, j, k
        real(fp):: dummy_dvxdt(2)

        if (.not. present(psys_rhs)) error stop "expected psys_rhs to be passed in."

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            ! apply boundary force with eqn 4.1.
            call repulsive_force( &
                dx, psys_lhs%particles(i)%c, psys_lhs%particles(i)%x(:), psys_rhs%particles(j)%x(:), &
                psys_lhs%particles(i)%dvxdt(:) &
                )
            ! boundary particles included in artificial viscosity calculation (start of pg 402), but velocities of boundary
            ! particles aren't updated.
            call artificial_viscosity_monaghan1994( &
                psys_lhs%particles(i)%x(:), psys_rhs%particles(j)%x(:), psys_lhs%particles(i)%v(:), psys_rhs%particles(j)%v(:), &
                psys_lhs%particles(i)%rho, psys_rhs%particles(j)%rho, self%h, self%h, psys_lhs%particles(i)%c, &
                psys_rhs%particles(j)%c, psys_lhs%particles(i)%mass, psys_rhs%particles(j)%mass, psys_lhs%particles(i)%dvxdt(:), &
                dummy_dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                )
        end do

    end subroutine fluid_boundary_sweep_new

end module grasph_monaghan1994

program main

    use grasph_monaghan1994

    use grasph_particle_system_m, only: particle_system_t
    use weakly_compressible_particles_m, only: tait_eos_state_updater_t
    use grasph_system_interactions_m, only: system_interaction_t
    use grasph_time_integration, only: leap_frog_time_integration
    use grasph_kernels_m, only: cubic_bspline_kernel_t

    implicit none
    type(particle_system_t):: psys(2)
    type(system_interaction_t):: psys_interactions(2)
    type(cubic_bspline_kernel_t):: kernel
    integer:: i, j, k
    real(fp):: analytical_pressure
    type(fluid_sweeper_t):: self_sweeper
    type(fluid_boundary_sweeper_t):: boundary_sweeper
    type(xsph_shifter_t):: shifter
    type(eos_particle_t):: ps_template
    type(tait_eos_state_updater_t):: state_updater

    ! init fluid particles
    state_updater%rho_ref = rho0
    call psys(1)%base_init(n=2500, name="fluid", particle_template=ps_template, state_updater=state_updater)

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
    call psys(2)%register_io%register_variable(psys(2)%particles(1), "x", psys(2)%particles(1)%x)
    call psys(2)%register_io%register_variable(psys(2)%particles(1), "v", psys(2)%particles(1)%v)
    call psys(2)%register_io%register_variable(psys(2)%particles(1), "rho", psys(2)%particles(1)%rho)
    call psys(2)%register_io%register_variable(psys(2)%particles(1), "mass", psys(2)%particles(1)%mass)
    call psys(2)%register_io%register_variable(psys(2)%particles(1), "c", psys(2)%particles(1)%c)
    call psys(2)%register_io%register_variable(psys(2)%particles(1), "dvxdt", psys(2)%particles(1)%dvxdt)
    call psys(2)%register_io%register_variable(psys(2)%particles(1), "drhodt", psys(2)%particles(1)%drhodt)

    do i = 0, nfx - 1
        do j = 0, nfy - 1
            k = j*nfx + i + 1
            psys(1)%particles(k)%id = k
            psys(1)%particles(k)%type = 1 ! not sure if type is needed anymore
            psys(1)%particles(k)%x(1) = (i + 0.5_fp)*dx
            psys(1)%particles(k)%x(2) = (j + 0.5_fp)*dx
            psys(1)%particles(k)%c = 10._fp*sqrt(2._fp*abs(g)*25._fp) ! 10*sqrt(2gH) eqn 3.3
            ! initialize density of fluid particles using hydrostatic pressure condition (eqn 5.1)
            analytical_pressure = (25._fp - psys(1)%particles(k)%x(2))*rho0*abs(g)
            psys(1)%particles(k)%rho = rho0*(analytical_pressure*7._fp/(rho0*psys(1)%particles(k)%c**2) + 1._fp)**(1._fp/7._fp)
            psys(1)%particles(k)%mass = rho0*dx*dx
            psys(1)%particles(k)%v(:) = 0._fp
        end do
    end do

    ! init boundary particles
    ! use base_init since we're using the base type
    ! only need to initialize metadata and position as only position is used to calculate repulsive force
    call psys(2)%base_init(n=464, name="boundary")
    psys(2)%to_print_summary = .false.
    k = 0
    ! bottom layer and corners
    do i = -1, nbx
        k = k + 1
        psys(2)%particles(k)%x(1) = (i + 0.5_fp)*dx
        psys(2)%particles(k)%x(2) = -0.5_fp*dx
    end do
    ! top layer and corners
    do i = -1, nbx
        k = k + 1
        psys(2)%particles(k)%x(1) = (i + 0.5_fp)*dx
        psys(2)%particles(k)%x(2) = 40._fp + 0.5_fp*dx
    end do
    ! left wall
    do j = 0, nby - 1
        k = k + 1
        psys(2)%particles(k)%x(1) = -0.5_fp*dx
        psys(2)%particles(k)%x(2) = (j + 0.5_fp)*dx
    end do
    ! right wall
    do j = 0, nby - 1
        k = k + 1
        psys(2)%particles(k)%x(1) = 75._fp + 0.5_fp*dx
        psys(2)%particles(k)%x(2) = (j + 0.5_fp)*dx
    end do
    do i = 1, k
        psys(2)%particles(i)%id = i
        psys(2)%particles(i)%type = -1
        psys(2)%particles(i)%rho = rho0
        psys(2)%particles(i)%mass = rho0*dx*dx
        psys(2)%particles(i)%v(:) = 0._fp
        psys(2)%particles(i)%c = 10._fp*sqrt(2._fp*abs(g)*25._fp) ! 10*max_speed
    end do

    ! init interactions
    self_sweeper%artvisc_alpha = 0.01_fp
    self_sweeper%artvisc_beta = 0._fp
    self_sweeper%h = 1.2_fp*dx
    self_sweeper%g = g
    boundary_sweeper%artvisc_alpha = 0.01_fp
    boundary_sweeper%artvisc_beta = 0._fp
    boundary_sweeper%h = 1.2_fp*dx
    boundary_sweeper%g = g
    shifter%epsilon = 0.5_fp
    shifter%update_rhs = .true.
    call psys_interactions(1)%init(30, psys(1), sweeper=self_sweeper, shifter=shifter)
    shifter%update_rhs = .false.
    call psys_interactions(2)%init(30, psys(1), psys(2), sweeper=boundary_sweeper, shifter=shifter)

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

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

end program main
