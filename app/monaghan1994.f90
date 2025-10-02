!> @file monaghan1994.f90
!> @brief Module and program containing setup code to replicate the classic dambreak experiment from
!>        Monaghan (1994) (https://doi.org/10.1006/jcph.1994.1034). The only known difference is that the Leap-Frog time-integration
!>        is used here instead of the predictor-corrector scheme used in the paper.
!> @author Edward Yang
!> @date 2025/9/22

module grasph_monaghan1994

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles
    use weakly_compressible_particles, only: linear_eos_particle, tait_eos_state_updater
    use grasph_pairs, only: particle_pairs
    use grasph_pair_sets, only: particle_interactions, base_sweeper
    use weakly_compressible_interactions, only: fluid_sweeper
    use grasph_pair_interactions, only: artificial_viscosity_monaghan1994, continuity_density, repulsive_force
    use grasph_particle_shifting, only: xsph_shifter

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.5_fp, g = -9.81_fp, rho0 = 1000._fp
    ! no. of particles in x, y direction in initial geometry of fluid
    integer, parameter:: nfx = 25._fp/dx, nfy = 25._fp/dx
    ! no. of particles in x, y direction for boundary
    integer, parameter:: nbx = 75._fp/dx, nby = 40._fp/dx

    type, extends(base_sweeper):: fluid_boundary_sweeper
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
    end type fluid_boundary_sweeper

contains

    subroutine fluid_boundary_sweep_new(self, pairs, ps_lhs, ps_rhs)
        class(fluid_boundary_sweeper), intent(in):: self
        type(particle_pairs), intent(in):: pairs
        class(base_particles), intent(inout):: ps_lhs
        class(base_particles), optional, intent(inout):: ps_rhs
        integer:: i, j, k
        real(fp):: dummy_dvxdt(2)

        if (.not. present(ps_rhs)) error stop "expected ps_rhs to be passed in."

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            ! apply boundary force with eqn 4.1.
            call repulsive_force(dx, ps_lhs%ps(i)%c, ps_lhs%ps(i)%x(:), ps_rhs%ps(j)%x(:), ps_lhs%ps(i)%dvxdt(:))
            ! boundary particles included in artificial viscosity calculation (start of pg 402), but velocities of boundary
            ! particles aren't updated.
            call artificial_viscosity_monaghan1994(ps_lhs%ps(i)%x(:), ps_rhs%ps(j)%x(:), ps_lhs%ps(i)%v(:), ps_rhs%ps(j)%v(:), &
                                                   ps_lhs%ps(i)%rho, ps_rhs%ps(j)%rho, self%h, self%h, ps_lhs%ps(i)%c, &
                                                   ps_rhs%ps(j)%c, ps_lhs%ps(i)%mass, ps_rhs%ps(j)%mass, ps_lhs%ps(i)%dvxdt(:), &
                                                   dummy_dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                                                   )
        end do

    end subroutine fluid_boundary_sweep_new

end module grasph_monaghan1994

program main

    use grasph_monaghan1994

    use grasph_particles, only: base_particles
    use weakly_compressible_particles, only: tait_eos_state_updater
    use grasph_pair_sets, only: particle_interactions
    use grasph_time_integration, only: leap_frog_time_integration
    use grasph_kernels, only: grasph_cubic_bspline_kernel

    implicit none
    type(base_particles):: ps(2)
    type(particle_interactions):: pic(2)
    type(grasph_cubic_bspline_kernel):: kernel
    integer:: i, j, k
    real(fp):: analytical_pressure
    type(fluid_sweeper):: self_sweeper
    type(fluid_boundary_sweeper):: boundary_sweeper
    type(xsph_shifter):: shifter
    type(linear_eos_particle):: ps_template
    type(tait_eos_state_updater):: state_updater

    ! init fluid particles
    state_updater%rho_ref = rho0
    call ps(1)%base_init(n=2500, name="fluid", ps_template=ps_template, state_updater=state_updater)

    ! register variables for time-update
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
    call ps(2)%register_io%register_variable(ps(2)%ps(1), "x", ps(2)%ps(1)%x)
    call ps(2)%register_io%register_variable(ps(2)%ps(1), "v", ps(2)%ps(1)%v)
    call ps(2)%register_io%register_variable(ps(2)%ps(1), "rho", ps(2)%ps(1)%rho)
    call ps(2)%register_io%register_variable(ps(2)%ps(1), "mass", ps(2)%ps(1)%mass)
    call ps(2)%register_io%register_variable(ps(2)%ps(1), "c", ps(2)%ps(1)%c)
    call ps(2)%register_io%register_variable(ps(2)%ps(1), "dvxdt", ps(2)%ps(1)%dvxdt)
    call ps(2)%register_io%register_variable(ps(2)%ps(1), "drhodt", ps(2)%ps(1)%drhodt)

    do i = 0, nfx - 1
        do j = 0, nfy - 1
            k = j*nfx + i + 1
            ps(1)%ps(k)%id = k
            ps(1)%ps(k)%type = 1 ! not sure if type is needed anymore
            ps(1)%ps(k)%x(1) = (i + 0.5_fp)*dx
            ps(1)%ps(k)%x(2) = (j + 0.5_fp)*dx
            ps(1)%ps(k)%c = 10._fp*sqrt(2._fp*abs(g)*25._fp) ! 10*sqrt(2gH) eqn 3.3
            ! initialize density of fluid particles using hydrostatic pressure condition (eqn 5.1)
            analytical_pressure = (25._fp - ps(1)%ps(k)%x(2))*rho0*abs(g)
            ps(1)%ps(k)%rho = rho0*(analytical_pressure*7._fp/(rho0*ps(1)%ps(k)%c**2) + 1._fp)**(1._fp/7._fp)
            ps(1)%ps(k)%mass = rho0*dx*dx
            ps(1)%ps(k)%v(:) = 0._fp
        end do
    end do

    ! init boundary particles
    ! use base_init since we're using the base type
    ! only need to initialize metadata and position as only position is used to calculate repulsive force
    call ps(2)%base_init(n=464, name="boundary")
    ps(2)%to_print_summary = .false.
    k = 0
    ! bottom layer and corners
    do i = -1, nbx
        k = k + 1
        ps(2)%ps(k)%x(1) = (i + 0.5_fp)*dx
        ps(2)%ps(k)%x(2) = -0.5_fp*dx
    end do
    ! top layer and corners
    do i = -1, nbx
        k = k + 1
        ps(2)%ps(k)%x(1) = (i + 0.5_fp)*dx
        ps(2)%ps(k)%x(2) = 40._fp + 0.5_fp*dx
    end do
    ! left wall
    do j = 0, nby - 1
        k = k + 1
        ps(2)%ps(k)%x(1) = -0.5_fp*dx
        ps(2)%ps(k)%x(2) = (j + 0.5_fp)*dx
    end do
    ! right wall
    do j = 0, nby - 1
        k = k + 1
        ps(2)%ps(k)%x(1) = 75._fp + 0.5_fp*dx
        ps(2)%ps(k)%x(2) = (j + 0.5_fp)*dx
    end do
    do i = 1, k
        ps(2)%ps(i)%id = i
        ps(2)%ps(i)%type = -1
        ps(2)%ps(i)%rho = rho0
        ps(2)%ps(i)%mass = rho0*dx*dx
        ps(2)%ps(i)%v(:) = 0._fp
        ps(2)%ps(i)%c = 10._fp*sqrt(2._fp*abs(g)*25._fp) ! 10*max_speed
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
    call pic(1)%init(30, ps(1), sweeper=self_sweeper, shifter=shifter)
    shifter%update_rhs = .false.
    call pic(2)%init(30, ps(1), ps(2), sweeper=boundary_sweeper, shifter=shifter)

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    call leap_frog_time_integration(100000, 1000, 1000, ps, pic, 0.05_fp, kernel, "/home/edwardy/test", output_comp_level=4)

end program main
