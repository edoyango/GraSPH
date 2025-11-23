!> @file weakly_compressible_interactions.f90
!> @brief Module containing subroutines for describing interactions between weakly compressible particles
!> @author Edward Yang
!> @date 2025-09-21
module weakly_compressible_interactions_m

    use grasph_constants_m, only: fp, ndims, max_name_len
    use grasph_particle_m, only: base_particles_t
    use weakly_compressible_particles_m, only: eos_particles_t, eos_ghost_particles_t, eos_viscous_stress_particles_t, &
                                               eos_viscous_stress_ghost_particles_t
    use grasph_system_interactions_m, only: base_sweeper_t
    use grasph_pairs_m, only: particle_pairs_t
    use grasph_pair_interactions_m, only: artificial_viscosity_monaghan1994, continuity_density, isotropic_pressure_force, &
                                          strain_rate, cauchy_stress_force, diffusion_density, repulsive_force

    implicit none

    private

    !> @brief Sweeper describing interaction between systems of weakly compressible particles.
    type, extends(base_sweeper_t):: fluid_sweeper_t
        !> @brief Acceleration due to gravity (m/s)
        real(fp):: g = -9.81_fp
        !> @brief Alpha coefficient for artificial viscosity.
        real(fp):: artvisc_alpha = 0.1_fp
        !> @brief Beta coefficient for artificial viscosity.
        real(fp):: artvisc_beta = 0.1_fp
        !> @brief Smoothing length to use fo artificial viscosity.
        real(fp):: h = 0._fp
    contains
        procedure:: sweep_1system => fluid_sweep_1system
        procedure:: sweep_2system => fluid_sweep_2system
        procedure:: sweep_2system_norhsupdate => fluid_sweep_2system_norhsupdate
        procedure, nopass:: name => fluid_sweeper_name
    end type fluid_sweeper_t

    !> @brief Sweeper describing interaction between a fluid system and boundary using Lennard-Jones repulsive force only.
    type, extends(fluid_sweeper_t):: fluid_boundary_sweeper_monaghan1994_t
        !> @brief The interaction length of the lennard-jones repulsive force.
        real(fp):: cutoff
    contains
        procedure:: sweep_2system => fluid_boundary_sweep_monaghan1994_2system
        procedure:: sweep_2system_norhsupdate => fluid_boundary_sweep_monaghan1994_2system
        procedure, nopass:: name => fluid_boundary_sweeper_monaghan1994_name
    end type fluid_boundary_sweeper_monaghan1994_t

    !> @brief Sweeper describing how to update weakly compressible virtual boundary particles' velocity and density.
    type, extends(base_sweeper_t):: boundary_update_sweeper_t
    contains
        procedure:: sweep_2system => boundary_update_sweep_2system
        procedure:: sweep_2system_norhsupdate => boundary_update_sweep_2system
        procedure, nopass:: name => boundary_update_sweeper_name
    end type boundary_update_sweeper_t

    !> @brief Sweeper describing how to generate weakly compressible ghost particles at the start of the time-step.
    type, extends(base_sweeper_t):: ghost_timestep_setuper_t
        !> @brief The distance between the boundary face to mirror particles.
        real(fp):: cutoff
        !> @brief The normal to the face that particles should be mirrored on.
        real(fp):: surface_normal(ndims)
        !> @brief The point that the mirroring face passes through.
        real(fp):: point(ndims)
    contains
        procedure:: sweep_2system => ghost_timestep_setup_sweep_2system
        procedure:: sweep_2system_norhsupdate => ghost_timestep_setup_sweep_2system
        procedure, nopass:: name => ghost_timestep_setuper_name
    end type

    !> @brief Sweeper describing the interaction between real (LHS) and virtual (RHS) weakly compressible particles using the morris
    !>        boundary condition.
    type, extends(fluid_sweeper_t):: morris_boundary_sweeper_t
        !> @brief The point that the boundary plane passes through.
        real(fp):: point(ndims)
        !> @brief The unit normal vector pointing towards the "real" domain.
        real(fp):: normal(ndims)
    contains
        procedure:: sweep_2system => morris_boundary_sweep_2system
        procedure:: sweep_2system_norhsupdate => morris_boundary_sweep_2system
        procedure, nopass:: name => morris_boundary_sweeper_name
    end type morris_boundary_sweeper_t

    !> @brief Sweeper that calculates engineering strain rate between two particle systems.
    type, extends(base_sweeper_t):: strain_rate_sweeper_t
    contains
        procedure:: sweep_1system => strain_rate_sweep_1system
        procedure:: sweep_2system => strain_rate_sweep_2system
        procedure:: sweep_2system_norhsupdate => strain_rate_sweep_2system_norhsupdate
        procedure, nopass:: name => strain_rate_sweeper_name
    end type strain_rate_sweeper_t

    !> @brief Sweeper describing how to generate weakly compressible ghost particles with strain rate and cauchy stress tensors.
    type, extends(ghost_timestep_setuper_t):: eos_viscous_stress_ghost_timestep_setuper_t
    contains
        !> @brief Generates ghost particles based on the mirroring boundary and the RHS particles.
        procedure:: sweep_2system => eos_viscous_stress_ghost_timestep_setup_sweep_2system
        procedure:: sweep_2system_norhsupdate => eos_viscous_stress_ghost_timestep_setup_sweep_2system
        procedure, nopass:: name => eos_viscous_stress_ghost_timestep_setuper_name
    end type eos_viscous_stress_ghost_timestep_setuper_t

    !> @brief Sweeper describing interaction between systems of weakly compressible particles with cauchy stress tensors.
    type, extends(fluid_sweeper_t):: viscous_stress_fluid_sweeper_t
    contains
        procedure:: sweep_1system => viscous_stress_fluid_sweep_1system
        procedure:: sweep_2system => viscous_stress_fluid_sweep_2system
        procedure:: sweep_2system_norhsupdate => viscous_stress_fluid_sweep_2system_norhsupdate
        procedure, nopass:: name => viscous_stress_fluid_sweeper_name
    end type viscous_stress_fluid_sweeper_t

    !> @brief Sweeper describing the interaction between real (LHS) and virtual (RHS) weakly compressible articles (with cauchy
    !>        stress tensor) using the morris boundary condition.
    type, extends(viscous_stress_fluid_sweeper_t):: eos_viscous_stress_morris_boundary_sweeper_t
        !> @brief The point that the boundary plane passes through.
        real(fp):: point(ndims)
        !> @brief The unit normal vector pointing towards the "real" domain.
        real(fp):: normal(ndims)
    contains
        procedure:: sweep_2system => viscous_stress_morris_boundary_sweep_2system
        procedure:: sweep_2system_norhsupdate => viscous_stress_morris_boundary_sweep_2system
        procedure, nopass:: name => viscous_stress_morris_boundary_sweeper_name
    end type eos_viscous_stress_morris_boundary_sweeper_t

    !> @brief Sweeper that describes how morris boundary particles contribute to the strain rate calculation of real particles.
    type, extends(strain_rate_sweeper_t):: strain_rate_morris_boundary_sweeper_t
        !> @brief The point that the boundary plane passes through.
        real(fp):: point(ndims)
        !> @brief The unit normal vector pointing towards the "real" domain.
        real(fp):: normal(ndims)
    contains
        procedure:: sweep_2system => strain_rate_morris_boundary_sweep_2system
        procedure:: sweep_2system_norhsupdate => strain_rate_morris_boundary_sweep_2system
        procedure, nopass:: name => strain_rate_morris_boundary_sweeper_name
    end type strain_rate_morris_boundary_sweeper_t

    public:: fluid_sweeper_t, fluid_boundary_sweeper_monaghan1994_t, boundary_update_sweeper_t, ghost_timestep_setuper_t, &
             morris_boundary_sweeper_t, strain_rate_sweeper_t, viscous_stress_fluid_sweeper_t, &
             eos_viscous_stress_ghost_timestep_setuper_t, eos_viscous_stress_morris_boundary_sweeper_t, &
             strain_rate_morris_boundary_sweeper_t

contains

    !> @brief For one weakly-compressible fluid particle systems, calculate acceleration and density rate-of-change
    !>        due to isotropic pressure, artificial viscosity, and mass continuity.
    !> @param self The sweeper class holding artificial viscosity constants and gravity.
    !> @param pairs The class storing particle pair index information.
    !> @param ps the weakly compressible particles involved in the interactions.
    !> @param dt The time-step size.
    subroutine fluid_sweep_1system(self, pairs, ps, dt)
        class(fluid_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps
        real(fp), optional, intent(in):: dt
        class(eos_particles_t), pointer:: fluid
        integer:: i, j, k

        ! point to lhs particlse for access to pressure member
        select type (ps => ps)
        class is (eos_particles_t)
            fluid => ps
        class default
            error stop "Invalid type for ps"
        end select

        ! intialize LHS acceleration and density rate-of-change arrays
        if (self%initialise) then
            do i = 1, fluid%size
                fluid%dvxdt(:, i) = 0._fp
                fluid%dvxdt(ndims, i) = self%g
            end do
            fluid%drhodt(:) = 0._fp
        end if

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call artificial_viscosity_monaghan1994( &
                fluid%x(:, i), fluid%x(:, j), fluid%v(:, i), fluid%v(:, j), fluid%rho(i), fluid%rho(j), self%h, self%h, &
                fluid%c(i), fluid%c(j), fluid%mass(i), fluid%mass(j), fluid%dvxdt(:, i), fluid%dvxdt(:, j), pairs%dwdx(:, k), &
                self%artvisc_alpha, self%artvisc_beta &
                )
            call isotropic_pressure_force( &
                fluid%p(i), fluid%p(j), fluid%rho(i), fluid%rho(j), fluid%mass(i), fluid%mass(j), fluid%dvxdt(:, i), &
                fluid%dvxdt(:, j), pairs%dwdx(:, k) &
                )
            call continuity_density( &
                fluid%v(:, i), fluid%v(:, j), fluid%mass(i), fluid%mass(j), fluid%drhodt(i), fluid%drhodt(j), pairs%dwdx(:, k) &
                )
        end do

    end subroutine fluid_sweep_1system

    !> @brief For two weakly-compressible fluid particle systems, calculate acceleration and density rate-of-change
    !>        due to isotropic pressure, artificial viscosity, and mass continuity.
    !> @param self The sweeper class holding artificial viscosity constants and gravity.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS weakly compressible particles involved in the interactions.
    !> @param ps_rhs Ths RHS "                                                        ".
    !> @param dt The time-step size.
    subroutine fluid_sweep_2system(self, pairs, ps_lhs, ps_rhs, dt)
        class(fluid_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        class(eos_particles_t), pointer:: fluid_lhs, fluid_rhs
        integer:: i, j, k

        ! point to lhs particlse for access to pressure member
        select type (ps => ps_lhs)
        class is (eos_particles_t)
            fluid_lhs => ps
        class default
            error stop "Invalid type for ps_lhs"
        end select

        ! point to rhs particlse for access to pressure member
        select type (ps => ps_rhs)
        class is (eos_particles_t)
            fluid_rhs => ps
        class default
            error stop "Invalid type for ps_rhs"
        end select

        ! intialize L/RHS acceleration and density rate-of-change arrays
        if (self%initialise) then
            do i = 1, fluid_lhs%size
                fluid_lhs%dvxdt(:, i) = 0._fp
                fluid_lhs%dvxdt(ndims, i) = self%g
            end do
            fluid_lhs%drhodt(:) = 0._fp
            do i = 1, fluid_rhs%size
                fluid_rhs%dvxdt(:, i) = 0._fp
                fluid_rhs%dvxdt(ndims, i) = self%g
            end do
            fluid_rhs%drhodt(:) = 0._fp
        end if

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call artificial_viscosity_monaghan1994( &
                fluid_lhs%x(:, i), fluid_rhs%x(:, j), fluid_lhs%v(:, i), fluid_rhs%v(:, j), fluid_lhs%rho(i), &
                fluid_rhs%rho(j), self%h, self%h, fluid_lhs%c(i), fluid_rhs%c(j), fluid_lhs%mass(i), fluid_rhs%mass(j), &
                fluid_lhs%dvxdt(:, i), fluid_rhs%dvxdt(:, j), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                )
            call isotropic_pressure_force( &
                fluid_lhs%p(i), fluid_rhs%p(j), fluid_lhs%rho(i), fluid_rhs%rho(j), fluid_lhs%mass(i), fluid_rhs%mass(j), &
                fluid_lhs%dvxdt(:, i), fluid_rhs%dvxdt(:, j), pairs%dwdx(:, k) &
                )
            call continuity_density( &
                fluid_lhs%v(:, i), fluid_rhs%v(:, j), fluid_lhs%mass(i), fluid_rhs%mass(j), fluid_lhs%drhodt(i), &
                fluid_rhs%drhodt(j), pairs%dwdx(:, k) &
                )
        end do

    end subroutine fluid_sweep_2system

    !> @brief For two weakly-compressible fluid particle systems, calculate acceleration and density rate-of-change
    !>        due to isotropic pressure, artificial viscosity, and mass continuity.
    !> @param self The sweeper class holding artificial viscosity constants and gravity.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS weakly compressible particles involved in the interactions.
    !> @param ps_rhs Ths RHS "                                                        ".
    !> @param dt The time-step size.
    subroutine fluid_sweep_2system_norhsupdate(self, pairs, ps_lhs, ps_rhs, dt)
        class(fluid_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        class(eos_particles_t), pointer:: fluid_lhs, fluid_rhs
        integer:: i, j, k
        real(fp):: dummy_drhodt, dummy_dvxdt(ndims) ! dummy variables for when update_rhs is .false.

        ! point to lhs particlse for access to pressure member
        select type (ps => ps_lhs)
        class is (eos_particles_t)
            fluid_lhs => ps
        class default
            error stop "Invalid type for ps_lhs"
        end select

        ! intialize LHS acceleration and density rate-of-change arrays
        if (self%initialise) then
            do i = 1, fluid_lhs%size
                fluid_lhs%dvxdt(:, i) = 0._fp
                fluid_lhs%dvxdt(ndims, i) = self%g
            end do
            fluid_lhs%drhodt(:) = 0._fp
        end if

        ! point to rhs particlse for access to pressure member
        select type (ps => ps_rhs)
        class is (eos_particles_t)
            fluid_rhs => ps
        class default
            error stop "Invalid type for ps_rhs"
        end select

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call artificial_viscosity_monaghan1994( &
                fluid_lhs%x(:, i), fluid_rhs%x(:, j), fluid_lhs%v(:, i), fluid_rhs%v(:, j), fluid_lhs%rho(i), &
                fluid_rhs%rho(j), self%h, self%h, fluid_lhs%c(i), fluid_rhs%c(j), fluid_lhs%mass(i), fluid_rhs%mass(j), &
                fluid_lhs%dvxdt(:, i), dummy_dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                )
            call isotropic_pressure_force( &
                fluid_lhs%p(i), fluid_rhs%p(j), fluid_lhs%rho(i), fluid_rhs%rho(j), fluid_lhs%mass(i), fluid_rhs%mass(j), &
                fluid_lhs%dvxdt(:, i), dummy_dvxdt(:), pairs%dwdx(:, k) &
                )
            call continuity_density( &
                fluid_lhs%v(:, i), fluid_rhs%v(:, j), fluid_lhs%mass(i), fluid_rhs%mass(j), &
                fluid_lhs%drhodt(i), dummy_drhodt, pairs%dwdx(:, k) &
                )
        end do

    end subroutine fluid_sweep_2system_norhsupdate

    !> @brief For a weakly-compressible fluid particle system and a base particle system, calculate acceleration due to
    !>        lennard-jones repulsive force and artificial viscosity.
    !> @param self The sweeper class holding artificial viscosity constants.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the real particles involved in the interactions.
    !> @param ps_rhs The boundary particles involved in the interactions.
    !> @param dt The time-step size.
    subroutine fluid_boundary_sweep_monaghan1994_2system(self, pairs, ps_lhs, ps_rhs, dt)
        class(fluid_boundary_sweeper_monaghan1994_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        real(fp):: dummy_dvxdt(2)

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            ! apply boundary force with eqn 4.1.
            call repulsive_force(self%cutoff, ps_lhs%c(i), ps_lhs%x(:, i), ps_rhs%x(:, j), ps_lhs%dvxdt(:, i))
            ! boundary particles included in artificial viscosity calculation (start of pg 402), but velocities of boundary
            ! particles aren't updated.
            call artificial_viscosity_monaghan1994( &
                ps_lhs%x(:, i), ps_rhs%x(:, j), ps_lhs%v(:, i), ps_rhs%v(:, j), ps_lhs%rho(i), ps_rhs%rho(j), self%h, self%h, &
                ps_lhs%c(i), ps_rhs%c(j), ps_lhs%mass(i), ps_rhs%mass(j), ps_lhs%dvxdt(:, i), dummy_dvxdt(:), pairs%dwdx(:, k), &
                self%artvisc_alpha, self%artvisc_beta &
                )
        end do

    end subroutine fluid_boundary_sweep_monaghan1994_2system

    !> @brief For two particle systems, calculate the LHS particles' density and velocity using kernel interpolation over the RHS
    !>        particles' density and velocity.
    !> @param self The sweeper class.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS particles whose velocity and density will be updated.
    !> @param ps_rhs Ths RHS particles to calculate velocity and density from.
    !> @param dt The time-step size.
    subroutine boundary_update_sweep_2system(self, pairs, ps_lhs, ps_rhs, dt)
        class(boundary_update_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        real(fp):: mw, vw
        real(fp), allocatable:: wsum(:)

        allocate (wsum(ps_rhs%size), source=0._fp)

        ps_rhs%rho(:) = 0._fp
        ps_rhs%v(:, :) = 0._fp

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            mw = ps_lhs%mass(i)*pairs%w(k)
            vw = mw/ps_lhs%rho(i)
            wsum(j) = wsum(j) + vw
            ps_rhs%rho(j) = ps_rhs%rho(j) + mw
            ps_rhs%v(:, j) = ps_rhs%v(:, j) + ps_lhs%v(:, i)*vw
        end do

        do j = 1, ps_rhs%size
            if (wsum(j) > 0._fp) then
                ps_rhs%v(:, j) = -ps_rhs%v(:, j)/wsum(j)
                ps_rhs%rho(j) = ps_rhs%rho(j)/wsum(j)
            end if
        end do

    end subroutine boundary_update_sweep_2system

    !> @brief For two weakly compressible particle systems, generate ghost particles in the LHS system, using information from the
    !>        RHS system and a defined boundary.
    !> @param self The sweeper class containing the plane boundary using a unit normal vector and point.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS ghost particle system where ghost particles will be generated.
    !> @param ps_rhs Ths RHS particles to generate ghost particles from.
    !> @param dt The time-step size.
    subroutine ghost_timestep_setup_sweep_2system(self, pairs, ps_lhs, ps_rhs, dt)
        class(ghost_timestep_setuper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        integer:: i
        class(eos_particles_t), pointer:: ps_real
        class(eos_ghost_particles_t), pointer:: ps_ghost
        real(fp):: dx(ndims), dr

        select type (ps => ps_lhs)
        class is (eos_particles_t)
            ps_real => ps
        class default
            error stop "Expected ps_lhs to be eos_particles_t."
        end select

        select type (ps => ps_rhs)
        class is (eos_ghost_particles_t)
            ps_ghost => ps
        class default
            error stop "Expected ps_rhs to be eos_ghost_particles_t."
        end select

        ps_ghost%size = 0

        ps_ghost%ps_original => ps_real

        do i = 1, ps_lhs%size
            dr = dot_product(ps_real%x(:, i) - self%point(:), self%surface_normal(:))
            ! ghost particle if within cutoff and inside the modelled region.
            if (abs(dr) <= self%cutoff .and. dr > 0._fp) then
                ps_ghost%size = ps_ghost%size + 1
                ps_ghost%id(ps_ghost%size) = ps_ghost%size
                ps_ghost%x(:, ps_ghost%size) = ps_real%x(:, i) - 2._fp*dr*self%surface_normal(:)
                ps_ghost%idx_original(ps_ghost%size) = i
            end if

        end do

    end subroutine ghost_timestep_setup_sweep_2system

    !> @brief For a weakly compressible particle system and base particle system (aka boundary particles), calculate the
    !>        contribution of the boundary particles to the acceleration and rate of change to the real particles using the morris
    !>        boundary.
    !> @param self The sweeper class holding artificial viscosity constants and boundary surface point and unit normal vector.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS weakly compressible particles involved in the interactions.
    !> @param ps_rhs Ths RHS boundary particles involved in the interactions.
    !> @param dt The time-step size.
    subroutine morris_boundary_sweep_2system(self, pairs, ps_lhs, ps_rhs, dt)
        class(morris_boundary_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        class(eos_particles_t), pointer:: ps_fluid
        real(fp):: dummy_drhodt, dummy_dvxdt(ndims), vb(ndims), da, db

        select type (ps => ps_lhs)
        class is (eos_particles_t)
            ps_fluid => ps
        class default
            error stop "Expected ps_lhs to be eos_particles_t."
        end select

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            da = dot_product(ps_fluid%x(:, i) - self%point(:), self%normal(:))
            db = dot_product(ps_rhs%x(:, j) - self%point(:), self%normal(:))
            vb(:) = -min(3._fp, abs(db/da))*ps_fluid%v(:, i)
            call artificial_viscosity_monaghan1994( &
                ps_fluid%x(:, i), ps_rhs%x(:, j), ps_fluid%v(:, i), vb(:), ps_fluid%rho(i), ps_fluid%rho(i), self%h, self%h, &
                ps_fluid%c(i), ps_fluid%c(i), ps_fluid%mass(i), ps_fluid%mass(i), ps_fluid%dvxdt(:, i), dummy_dvxdt(:), &
                pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                )
            call isotropic_pressure_force( &
                ps_fluid%p(i), ps_fluid%p(i), ps_fluid%rho(i), ps_fluid%rho(i), ps_fluid%mass(i), ps_fluid%mass(i), &
                ps_fluid%dvxdt(:, i), dummy_dvxdt(:), pairs%dwdx(:, k) &
                )
            call continuity_density( &
                ps_fluid%v(:, i), vb(:), ps_fluid%mass(i), ps_fluid%mass(i), &
                ps_fluid%drhodt(i), dummy_drhodt, pairs%dwdx(:, k) &
                )
        end do

    end subroutine morris_boundary_sweep_2system

    !> @brief For either one or two particle systems, calculate strain rate.
    !> @param self The sweeper class.
    !> @param pairs The class storing particle pair index information.
    !> @param ps The particle system involved in the interactions.
    !> @param dt The time-step size.
    subroutine strain_rate_sweep_1system(self, pairs, ps, dt)

        use weakly_compressible_particles_m, only: ntensor_elems_voigt

        class(strain_rate_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps
        real(fp), optional, intent(in):: dt
        class(eos_viscous_stress_particles_t), pointer:: ps_vs
        integer:: i, j, k

        select type (ps_sr => ps)
        class is (eos_viscous_stress_particles_t)
            ps_vs => ps_sr
        class default
            error stop "Expected ps class to be eos_viscous_stress_particles_t."
        end select

        if (self%initialise) then
            ps_vs%strain_rate(:, :) = 0._fp
        end if

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call strain_rate( &
                ps_vs%v(:, i), ps_vs%v(:, j), ps_vs%mass(i), ps_vs%mass(j), ps_vs%rho(i), ps_vs%rho(j), pairs%dwdx(:, k), &
                ps_vs%strain_rate(:, i), ps_vs%strain_rate(:, j) &
                )
        end do

    end subroutine strain_rate_sweep_1system

    !> @brief For either one or two particle systems, calculate strain rate.
    !> @param self The sweeper class.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS particles involved in the interactions.
    !> @param ps_rhs Ths RHS "                                    ".
    !> @param dt The time-step size.
    subroutine strain_rate_sweep_2system(self, pairs, ps_lhs, ps_rhs, dt)

        use weakly_compressible_particles_m, only: ntensor_elems_voigt

        class(strain_rate_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        class(eos_viscous_stress_particles_t), pointer:: ps_vs_lhs, ps_vs_rhs
        real(fp):: dummy_strain_rate(ntensor_elems_voigt)
        integer:: i, j, k

        select type (ps => ps_lhs)
        class is (eos_viscous_stress_particles_t)
            ps_vs_lhs => ps
        class default
            error stop "Expected ps_lhs class to be eos_viscous_stress_particles_t."
        end select

        select type (ps => ps_rhs)
        class is (eos_viscous_stress_particles_t)
            ps_vs_rhs => ps
        class default
            error stop "Expected ps_rhs class to be eos_viscous_stress_particles_t."
        end select

        if (self%initialise) then
            ps_vs_lhs%strain_rate(:, :) = 0._fp
            ps_vs_rhs%strain_rate(:, :) = 0._fp
        end if

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call strain_rate( &
                ps_vs_lhs%v(:, i), ps_vs_rhs%v(:, j), ps_vs_lhs%mass(i), ps_vs_rhs%mass(j), ps_vs_lhs%rho(i), ps_vs_rhs%rho(j), &
                pairs%dwdx(:, k), ps_vs_lhs%strain_rate(:, i), ps_vs_rhs%strain_rate(:, j) &
                )
        end do

    end subroutine strain_rate_sweep_2system

    !> @brief For either one or two particle systems, calculate strain rate.
    !> @param self The sweeper class.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS particles involved in the interactions.
    !> @param ps_rhs Ths RHS "                                    ".
    !> @param dt The time-step size.
    subroutine strain_rate_sweep_2system_norhsupdate(self, pairs, ps_lhs, ps_rhs, dt)

        use weakly_compressible_particles_m, only: ntensor_elems_voigt

        class(strain_rate_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        class(eos_viscous_stress_particles_t), pointer:: ps_vs_lhs, ps_vs_rhs
        real(fp):: dummy_strain_rate(ntensor_elems_voigt)
        integer:: i, j, k

        select type (ps => ps_lhs)
        class is (eos_viscous_stress_particles_t)
            ps_vs_lhs => ps
        class default
            error stop "Expected ps_lhs class to be eos_viscous_stress_particles_t."
        end select

        if (self%initialise) then
            do i = 1, ps_vs_lhs%size
                ps_vs_lhs%strain_rate(:, :) = 0._fp
            end do
        end if

        select type (ps => ps_rhs)
        class is (eos_viscous_stress_particles_t)
            ps_vs_rhs => ps
        class default
            error stop "Expected ps_rhs class to be eos_viscous_stress_particles_t."
        end select

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call strain_rate( &
                ps_vs_lhs%v(:, i), ps_vs_rhs%v(:, j), ps_vs_lhs%mass(i), ps_vs_rhs%mass(j), ps_vs_lhs%rho(i), ps_vs_rhs%rho(j), &
                pairs%dwdx(:, k), ps_vs_lhs%strain_rate(:, i), dummy_strain_rate(:) &
                )
        end do

    end subroutine strain_rate_sweep_2system_norhsupdate

    !> @brief For two weakly compressible particle systems with stress tensor, generate ghost particles in the LHS system, using
    !>        information from the RHS system and a defined boundary.
    !> @param self The sweeper class containing the plane boundary using a unit normal vector and point.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS ghost particle system where ghost particles will be generated.
    !> @param ps_rhs Ths RHS particles to generate ghost particles from.
    !> @param dt The time-step size.
    subroutine eos_viscous_stress_ghost_timestep_setup_sweep_2system(self, pairs, ps_lhs, ps_rhs, dt)
        class(eos_viscous_stress_ghost_timestep_setuper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        integer:: i
        class(eos_viscous_stress_particles_t), pointer:: ps_real
        class(eos_viscous_stress_ghost_particles_t), pointer:: ps_ghost
        real(fp):: dr

        select type (ps => ps_lhs)
        class is (eos_viscous_stress_particles_t)
            ps_real => ps
        class default
            error stop "Expected ps_lhs to be eos_viscous_stress_particles_t."
        end select

        select type (ps => ps_rhs)
        class is (eos_viscous_stress_ghost_particles_t)
            ps_ghost => ps
        class default
            error stop "Expected ps_rhs to be eos_viscous_stress_ghost_particles_t."
        end select

        ps_rhs%size = 0

        ps_ghost%ps_original => ps_real

        do i = 1, ps_lhs%size
            dr = dot_product(ps_real%x(:, i) - self%point(:), self%surface_normal(:))
            ! ghost particle if within cutoff and inside the modelled region.
            if (abs(dr) <= self%cutoff .and. dr > 0._fp) then
                ps_rhs%size = ps_rhs%size + 1
                ps_ghost%id(ps_rhs%size) = ps_rhs%size
                ps_ghost%x(:, ps_rhs%size) = ps_real%x(:, i) - 2._fp*dr*self%surface_normal(:)
                ps_ghost%idx_original(ps_rhs%size) = i
            end if

        end do

    end subroutine eos_viscous_stress_ghost_timestep_setup_sweep_2system

    !> @brief For one weakly-compressible fluid particle system with stress tensor, calculate acceleration and density
    !>        rate-of-change due to cauchy stress, artificial viscosity, mass continuity, and density diffusion.
    !> @param self The sweeper class holding artificial viscosity constants and gravity.
    !> @param pairs The class storing particle pair index information.
    !> @param ps the weakly compressible particles with stress tensor involved in the interactions.
    !> @param dt The time-step size.
    subroutine viscous_stress_fluid_sweep_1system(self, pairs, ps, dt)
        class(viscous_stress_fluid_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps
        class(eos_viscous_stress_particles_t), pointer:: fluid
        real(fp), optional, intent(in):: dt
        integer:: i, j, k

        ! point to lhs particlse for access to pressure member
        select type (ps => ps)
        class is (eos_viscous_stress_particles_t)
            fluid => ps
        class default
            error stop "Invalid type for ps"
        end select

        ! intialize LHS acceleration and density rate-of-change arrays
        if (self%initialise) then
            do i = 1, ps%size
                fluid%dvxdt(:, i) = 0._fp
                fluid%dvxdt(ndims, i) = self%g
            end do
            fluid%drhodt(:) = 0._fp
        end if

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call artificial_viscosity_monaghan1994( &
                fluid%x(:, i), fluid%x(:, j), fluid%v(:, i), fluid%v(:, j), fluid%rho(i), fluid%rho(j), self%h, self%h, &
                fluid%c(i), fluid%c(j), fluid%mass(i), fluid%mass(j), fluid%dvxdt(:, i), fluid%dvxdt(:, j), pairs%dwdx(:, k), &
                self%artvisc_alpha, self%artvisc_beta &
                )
            call cauchy_stress_force( &
                fluid%stress(:, i), fluid%stress(:, j), fluid%rho(i), fluid%rho(j), fluid%mass(i), fluid%mass(j), &
                fluid%dvxdt(:, i), fluid%dvxdt(:, j), pairs%dwdx(:, k) &
                )
            call continuity_density( &
                fluid%v(:, i), fluid%v(:, j), fluid%mass(i), fluid%mass(j), fluid%drhodt(i), fluid%drhodt(j), pairs%dwdx(:, k) &
                )
            call diffusion_density( &
                fluid%rho(i), fluid%rho(j), fluid%x(:, i), fluid%x(:, j), fluid%mass(i), fluid%mass(j), self%h, self%h, &
                fluid%c(i), fluid%c(j), pairs%dwdx(:, k), fluid%drhodt(i), fluid%drhodt(j) &
                )
        end do

    end subroutine viscous_stress_fluid_sweep_1system

    !> @brief For two weakly-compressible fluid particle systems with stress tensor, calculate acceleration and density
    !>        rate-of-change due to cauchy stress, artificial viscosity, mass continuity, and density diffusion.
    !> @param self The sweeper class holding artificial viscosity constants and gravity.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS weakly compressible particles with stress tensor involved in the interactions.
    !> @param ps_rhs Ths RHS "                                                                           ".
    !> @param dt The time-step size.
    subroutine viscous_stress_fluid_sweep_2system(self, pairs, ps_lhs, ps_rhs, dt)
        class(viscous_stress_fluid_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        class(eos_viscous_stress_particles_t), pointer:: fluid_lhs, fluid_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k

        ! point to lhs particlse for access to pressure member
        select type (ps => ps_lhs)
        class is (eos_viscous_stress_particles_t)
            fluid_lhs => ps
        class default
            error stop "Invalid type for ps_lhs"
        end select

        ! point to rhs particlse for access to pressure member
        select type (ps => ps_rhs)
        class is (eos_viscous_stress_particles_t)
            fluid_rhs => ps
        class default
            error stop "Invalid type for ps_rhs"
        end select

        ! intialize L/RHS acceleration and density rate-of-change arrays
        if (self%initialise) then
            do i = 1, ps_lhs%size
                fluid_lhs%dvxdt(:, i) = 0._fp
                fluid_lhs%dvxdt(ndims, i) = self%g
            end do
            fluid_lhs%drhodt(:) = 0._fp
            do i = 1, fluid_rhs%size
                fluid_rhs%dvxdt(:, i) = 0._fp
                fluid_rhs%dvxdt(ndims, i) = self%g
            end do
            fluid_rhs%drhodt(:) = 0._fp
        end if

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call artificial_viscosity_monaghan1994( &
                fluid_lhs%x(:, i), fluid_rhs%x(:, j), fluid_lhs%v(:, i), fluid_rhs%v(:, j), fluid_lhs%rho(i), &
                fluid_rhs%rho(j), self%h, self%h, fluid_lhs%c(i), fluid_rhs%c(j), fluid_lhs%mass(i), fluid_rhs%mass(j), &
                fluid_lhs%dvxdt(:, i), fluid_rhs%dvxdt(:, j), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                )
            call cauchy_stress_force( &
                fluid_lhs%stress(:, i), fluid_rhs%stress(:, j), fluid_lhs%rho(i), fluid_rhs%rho(j), fluid_lhs%mass(i), &
                fluid_rhs%mass(j), fluid_lhs%dvxdt(:, i), fluid_rhs%dvxdt(:, j), pairs%dwdx(:, k) &
                )
            call continuity_density( &
                fluid_lhs%v(:, i), fluid_rhs%v(:, j), fluid_lhs%mass(i), fluid_rhs%mass(j), fluid_lhs%drhodt(i), &
                fluid_rhs%drhodt(j), pairs%dwdx(:, k) &
                )
            call diffusion_density( &
                fluid_lhs%rho(i), fluid_rhs%rho(j), fluid_lhs%x(:, i), fluid_rhs%x(:, j), fluid_lhs%mass(i), &
                fluid_rhs%mass(j), self%h, self%h, fluid_lhs%c(i), fluid_rhs%c(j), pairs%dwdx(:, k), fluid_lhs%drhodt(i), &
                fluid_rhs%drhodt(j) &
                )
        end do

    end subroutine viscous_stress_fluid_sweep_2system

    !> @brief For two weakly-compressible fluid particle systems with stress tensor, calculate acceleration and density
    !>        rate-of-change due to cauchy stress, artificial viscosity, mass continuity, and density diffusion.
    !> @param self The sweeper class holding artificial viscosity constants and gravity.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS weakly compressible particles with stress tensor involved in the interactions.
    !> @param ps_rhs Ths RHS "                                                                           ".
    !> @param dt The time-step size.
    subroutine viscous_stress_fluid_sweep_2system_norhsupdate(self, pairs, ps_lhs, ps_rhs, dt)
        class(viscous_stress_fluid_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        class(eos_viscous_stress_particles_t), pointer:: fluid_lhs, fluid_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        real(fp):: dummy_drhodt, dummy_dvxdt(ndims) ! dummy variables for when update_rhs is .false.

        ! point to lhs particlse for access to pressure member
        select type (ps => ps_lhs)
        class is (eos_viscous_stress_particles_t)
            fluid_lhs => ps
        class default
            error stop "Invalid type for ps_lhs"
        end select

        ! point to rhs particlse for access to pressure member
        select type (ps => ps_rhs)
        class is (eos_viscous_stress_particles_t)
            fluid_rhs => ps
        class default
            error stop "Invalid type for ps_rhs"
        end select

        ! intialize LHS acceleration and density rate-of-change arrays
        if (self%initialise) then
            do i = 1, ps_lhs%size
                fluid_lhs%dvxdt(:, i) = 0._fp
                fluid_lhs%dvxdt(ndims, :) = self%g
            end do
            fluid_lhs%drhodt(:) = 0._fp
        end if

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call artificial_viscosity_monaghan1994( &
                fluid_lhs%x(:, i), fluid_rhs%x(:, j), fluid_lhs%v(:, i), fluid_rhs%v(:, j), fluid_lhs%rho(i), &
                fluid_rhs%rho(j), self%h, self%h, fluid_lhs%c(i), fluid_rhs%c(j), fluid_lhs%mass(i), fluid_rhs%mass(j), &
                fluid_lhs%dvxdt(:, i), dummy_dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                )
            call cauchy_stress_force( &
                fluid_lhs%stress(:, i), fluid_rhs%stress(:, j), fluid_lhs%rho(i), fluid_rhs%rho(j), fluid_lhs%mass(i), &
                fluid_rhs%mass(j), fluid_lhs%dvxdt(:, i), dummy_dvxdt(:), pairs%dwdx(:, k) &
                )
            call continuity_density( &
                fluid_lhs%v(:, i), fluid_rhs%v(:, j), fluid_lhs%mass(i), fluid_rhs%mass(j), &
                fluid_lhs%drhodt(i), dummy_drhodt, pairs%dwdx(:, k) &
                )
            call diffusion_density( &
                fluid_lhs%rho(i), fluid_rhs%rho(j), fluid_lhs%x(:, i), fluid_rhs%x(:, j), fluid_lhs%mass(i), &
                fluid_rhs%mass(j), self%h, self%h, fluid_lhs%c(i), fluid_rhs%c(j), pairs%dwdx(:, k), fluid_lhs%drhodt(i), &
                dummy_drhodt &
                )
        end do

    end subroutine viscous_stress_fluid_sweep_2system_norhsupdate

    !> @brief For a weakly compressible particle system with stress tensor and base particle system (aka boundary particles),
    !>        calculate the contribution of the boundary particles to the acceleration and rate of change to the real particles
    !>        using the morris boundary.
    !> @param self The sweeper class holding artificial viscosity constants and boundary surface point and unit normal vector.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS weakly compressible particles with stress tensor involved in the interactions.
    !> @param ps_rhs Ths RHS boundary particles involved in the interactions.
    !> @param dt The time-step size.
    subroutine viscous_stress_morris_boundary_sweep_2system(self, pairs, ps_lhs, ps_rhs, dt)
        class(eos_viscous_stress_morris_boundary_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        class(eos_viscous_stress_particles_t), pointer:: ps_fluid
        real(fp):: dummy_drhodt, dummy_dvxdt(ndims), vb(ndims), da, db

        select type (ps => ps_lhs)
        class is (eos_viscous_stress_particles_t)
            ps_fluid => ps
        class default
            error stop "Expected ps_lhs to be eos_viscous_stress_particles_t."
        end select

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            da = dot_product(ps_fluid%x(:, i) - self%point(:), self%normal(:))
            db = dot_product(ps_rhs%x(:, j) - self%point(:), self%normal(:))
            vb(:) = -min(3._fp, abs(db/da))*ps_fluid%v(:, i)
            call artificial_viscosity_monaghan1994( &
                ps_fluid%x(:, i), ps_rhs%x(:, j), ps_fluid%v(:, i), vb(:), ps_fluid%rho(i), &
                ps_fluid%rho(i), self%h, self%h, ps_fluid%c(i), ps_fluid%c(i), ps_fluid%mass(i), ps_fluid%mass(i), &
                ps_fluid%dvxdt(:, i), dummy_dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                )
            call cauchy_stress_force( &
                ps_fluid%stress(:, i), ps_fluid%stress(:, i), ps_fluid%rho(i), ps_fluid%rho(i), ps_fluid%mass(i), &
                ps_fluid%mass(i), ps_fluid%dvxdt(:, i), dummy_dvxdt(:), pairs%dwdx(:, k) &
                )
            call continuity_density( &
                ps_fluid%v(:, i), vb(:), ps_fluid%mass(i), ps_fluid%mass(i), &
                ps_fluid%drhodt(i), dummy_drhodt, pairs%dwdx(:, k) &
                )
        end do

    end subroutine viscous_stress_morris_boundary_sweep_2system

    !> @brief For two base particle systems, calculate the contribution of the boundary particles (RHS) to the strain rate of the
    !>        real particles (LHS) using the morris boundary.
    !> @param self The sweeper class holding boundary surface point and unit normal vector.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS particles involved in the interactions.
    !> @param ps_rhs Ths RHS "                                    ".
    !> @param dt The time-step size.
    subroutine strain_rate_morris_boundary_sweep_2system(self, pairs, ps_lhs, ps_rhs, dt)

        use weakly_compressible_particles_m, only: ntensor_elems_voigt

        class(strain_rate_morris_boundary_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        class(eos_viscous_stress_particles_t), pointer:: ps_real
        real(fp):: da, db, vb(ndims), dummy_strain_rate(ntensor_elems_voigt)
        integer:: i, j, k

        select type (ps => ps_lhs)
        class is (eos_viscous_stress_particles_t)
            ps_real => ps
        class default
            error stop "Expected ps_lhs class to be eos_viscous_stress_particles_t."
        end select

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            da = dot_product(ps_real%x(:, i) - self%point(:), self%normal(:))
            db = dot_product(ps_rhs%x(:, j) - self%point(:), self%normal(:))
            vb(:) = -min(3._fp, abs(db/da))*ps_real%v(:, i)
            call strain_rate( &
                ps_real%v(:, i), vb(:), ps_real%mass(i), ps_real%mass(i), ps_real%rho(i), ps_real%rho(i), &
                pairs%dwdx(:, k), ps_real%strain_rate(:, i), dummy_strain_rate(:) &
                )
        end do

    end subroutine strain_rate_morris_boundary_sweep_2system

    pure character(max_name_len) function eos_viscous_stress_ghost_timestep_setuper_name()
        eos_viscous_stress_ghost_timestep_setuper_name = "eos_viscous_stress_ghost_timestep_setuper_t"
    end function eos_viscous_stress_ghost_timestep_setuper_name

    pure character(max_name_len) function viscous_stress_fluid_sweeper_name()
        viscous_stress_fluid_sweeper_name = "viscous_stress_fluid_sweeper_t"
    end function viscous_stress_fluid_sweeper_name

    pure character(max_name_len) function viscous_stress_morris_boundary_sweeper_name()
        viscous_stress_morris_boundary_sweeper_name = "viscous_stress_morris_boundary_sweeper_t"
    end function viscous_stress_morris_boundary_sweeper_name

    pure character(max_name_len) function viscous_stress_boundary_sweeper_name()
        viscous_stress_boundary_sweeper_name = "viscous_stress_boundary_sweeper_t"
    end function viscous_stress_boundary_sweeper_name

    pure character(max_name_len) function fluid_boundary_sweeper_monaghan1994_name()
        fluid_boundary_sweeper_monaghan1994_name = "fluid_boundary_sweeper_monaghan1994_t"
    end function fluid_boundary_sweeper_monaghan1994_name

    pure character(max_name_len) function morris_boundary_sweeper_name()
        morris_boundary_sweeper_name = "morris_boundary_sweeper_t"
    end function morris_boundary_sweeper_name

    pure character(max_name_len) function strain_rate_sweeper_name()
        strain_rate_sweeper_name = "strain_rate_sweeper_t"
    end function strain_rate_sweeper_name

    pure character(max_name_len) function ghost_timestep_setuper_name()
        ghost_timestep_setuper_name = "ghost_timestep_setuper_t"
    end function ghost_timestep_setuper_name

    pure character(max_name_len) function boundary_update_sweeper_name()
        boundary_update_sweeper_name = "boundary_update_sweeper_t"
    end function boundary_update_sweeper_name

    pure character(max_name_len) function fluid_sweeper_name()
        fluid_sweeper_name = "fluid_sweeper_t"
    end function fluid_sweeper_name

    pure character(max_name_len) function strain_rate_morris_boundary_sweeper_name()
        strain_rate_morris_boundary_sweeper_name = "strain_rate_morris_boundary_sweeper_t"
    end function strain_rate_morris_boundary_sweeper_name

end module weakly_compressible_interactions_m
