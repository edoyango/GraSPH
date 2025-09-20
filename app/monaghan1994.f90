module grasph_monaghan1994

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles, wcp => weakly_compressible_particles
    use grasph_pair_sets, only: particle_interactions_base
    use grasph_pair_interactions, only: artificial_viscosity_monaghan1994, continuity_density, isotropic_pressure_force, &
                                        repulsive_force

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.5_fp, g = -9.81_fp
    integer, parameter:: nfx = 25._fp/dx, nfy = 25._fp/dx, nbx = 75._fp/dx, nby = 40._fp/dx

    ! declare how fluid particles interact with themselves
    type, extends(particle_interactions_base):: fluid_self_interaction
    contains
        procedure:: sweep => fluid_sweep
    end type fluid_self_interaction

    ! define how fluid particles interact with boundary
    type, extends(particle_interactions_base):: fluid_boundary_interaction
    contains
        procedure:: sweep => fluid_boundary_sweep
    end type fluid_boundary_interaction

contains

    subroutine fluid_sweep(self)
        class(fluid_self_interaction), intent(inout):: self
        integer:: i, jj, j
        class(wcp), pointer:: ps_fluid

        select type (ps => self%ps_lhs)
        class is (wcp)
            ps_fluid => ps
        class default
            error stop 'Particles are not "weakly_compressible_particles"'
        end select

        do i = 1, ps_fluid%size
            ps_fluid%dvxdt(1, i) = 0._fp
            ps_fluid%dvxdt(2, i) = g
            ps_fluid%drhodt(i) = 0._fp
        end do

        do i = 1, self%pairs%n
            do jj = self%pairs%offsets(i) + 1, self%pairs%offsets(i + 1)
                j = self%pairs%rhs(jj)
                call artificial_viscosity_monaghan1994( &
                    2, ps_fluid%x(:, i), ps_fluid%x(:, j), ps_fluid%v(:, i), ps_fluid%v(:, j), ps_fluid%rho(i), ps_fluid%rho(j), &
                    1.2_fp*dx, 1.2_fp*dx, ps_fluid%c(i), ps_fluid%c(j), ps_fluid%mass(i), &
                    ps_fluid%mass(j), ps_fluid%dvxdt(:, i), ps_fluid%dvxdt(:, j), self%pairs%dwdx(:, jj), &
                    0.1_fp, 0.1_fp)
                call isotropic_pressure_force(2, ps_fluid%p(i), ps_fluid%p(j), ps_fluid%rho(i), ps_fluid%rho(j), &
                                              ps_fluid%mass(i), ps_fluid%mass(j), ps_fluid%dvxdt(:, i), &
                                              ps_fluid%dvxdt(:, j), self%pairs%dwdx(:, jj))
                call continuity_density(2, ps_fluid%v(:, i), ps_fluid%v(:, j), ps_fluid%mass(i), ps_fluid%mass(j), &
                                        ps_fluid%drhodt(i), ps_fluid%drhodt(j), self%pairs%dwdx(:, jj))
            end do
        end do

    end subroutine fluid_sweep

    subroutine fluid_boundary_sweep(self)
        class(fluid_boundary_interaction), intent(inout):: self
        integer:: i, jj, j

        do i = 1, self%pairs%n
            do jj = self%pairs%offsets(i) + 1, self%pairs%offsets(i + 1)
                j = self%pairs%rhs(jj)
                call repulsive_force(2, dx, self%ps_lhs%c(i), self%ps_lhs%x(:, i), self%ps_rhs%x(:, j), self%ps_lhs%dvxdt(:, i))
            end do
        end do

    end subroutine fluid_boundary_sweep

end module grasph_monaghan1994

program main

    use grasph_monaghan1994

    use grasph_particles, only: particles_container, bp => base_particles, wcp => weakly_compressible_particles
    use grasph_pair_sets, only: particle_interactions_container
    use grasph_time_integration, only: leap_frog_time_integration
    use grasph_kernels, only: grasph_cubic_bspline_kernel

    implicit none
    type(particles_container):: ps(2)
    type(particle_interactions_container):: pic(2)
    type(grasph_cubic_bspline_kernel):: kernel, kernel2
    integer:: i, j, k

    ! declare particles - fluid and boundary (repulsive force)
    allocate (wcp::ps(1)%p)
    allocate (bp::ps(2)%p)

    ! describe interacting particles - fluid with themselves, and fluid with the boundary
    allocate (fluid_self_interaction::pic(1)%pi)
    allocate (fluid_boundary_interaction::pic(2)%pi)

    ! init fluid particles
    select type (ps => ps(1)%p) ! specialise for weakly compressible particles
    class is (wcp)
        call ps%init(n=2500, d=2, name="fluid", rho_ref=1000._fp)
    end select
    do i = 0, nfx - 1
        do j = 0, nfy - 1
            k = j*nfx + i + 1
            ps(1)%p%id(k) = k
            ps(1)%p%type = 1 ! not sure if type is needed anymore
            ps(1)%p%x(1, k) = (i + 0.5_fp)*dx
            ps(1)%p%x(2, k) = (j + 0.5_fp)*dx
            ps(1)%p%rho(k) = 1000._fp
            ps(1)%p%mass(k) = 1000._fp*dx*dx
            ps(1)%p%c(k) = 10._fp*2._fp*sqrt(abs(g)*25._fp) ! 10*max_speed
            ps(1)%p%v(:, k) = 0._fp
        end do
    end do

    ! init boundary particles
    ! use base_init since we're using the base type
    ! only need to initialize metadata and position as only position is used to calculate repulsive force
    call ps(2)%p%base_init(n=464, d=2, name="boundary")
    ps(2)%p%to_print_summary = .false.
    k = 0
    ! bottom layer and corners
    do i = -1, nbx
        k = k + 1
        ps(2)%p%id(k) = k
        ps(2)%p%type(k) = -1
        ps(2)%p%x(1, k) = (i + 0.5_fp)*dx
        ps(2)%p%x(2, k) = -0.5_fp*dx
    end do
    ! top layer and corners
    do i = -1, nbx
        k = k + 1
        ps(2)%p%id(k) = k
        ps(2)%p%type(k) = -1
        ps(2)%p%x(1, k) = (i + 0.5_fp)*dx
        ps(2)%p%x(2, k) = 40._fp + 0.5_fp*dx
    end do
    ! left wall
    do j = 0, nby - 1
        k = k + 1
        ps(2)%p%id(k) = k
        ps(2)%p%type(k) = -1
        ps(2)%p%x(1, k) = -0.5_fp*dx
        ps(2)%p%x(2, k) = (j + 0.5_fp)*dx
    end do
    ! right wall
    do j = 0, nby - 1
        k = k + 1
        ps(2)%p%id(k) = k
        ps(2)%p%type(k) = -1
        ps(2)%p%x(1, k) = 75._fp + 0.5_fp*dx
        ps(2)%p%x(2, k) = (j + 0.5_fp)*dx
    end do

    ! init interactions
    call pic(1)%pi%base_init(30, ps(1)%p)
    call pic(2)%pi%base_init(30, ps(1)%p, ps(2)%p)

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    call leap_frog_time_integration(100000, 1000, 1000, ps, pic, 0.1_fp, kernel, "/home/edwardy/test", "", 4)

end program main
