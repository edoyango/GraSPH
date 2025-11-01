!> @file weakly_compressible_interactions.f90
!> @brief Module containing subroutines for describing interactions between weakly compressible particles
!> @author Edward Yang
!> @date 2025-09-21
module weakly_compressible_interactions_m

    use grasph_constants_m, only: fp, ndims
    use grasph_particle_system_m, only: particle_system_t
    use weakly_compressible_particles_m, only: eos_particle_t, eos_ghost_particle_t, eos_viscous_stress_particle_t, &
                                               eos_viscous_stress_ghost_particle_t
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
        !> @brief For a single set of weakly-compressible fluid particles, calculate acceleration and density rate-of-change due to
        !>        isotropic pressure, artificial viscosity, and mass continuity.
        procedure:: sweep => fluid_sweep
    end type fluid_sweeper_t

    !> @brief Sweeper describing interaction between a fluid system and boundary using Lennard-Jones repulsive force only.
    type, extends(fluid_sweeper_t):: fluid_boundary_sweeper_monaghan1994_t
        !> @brief The interaction length of the lennard-jones repulsive force.
        real(fp):: cutoff
    contains
        !> @brief Calculate acceleration due to lennard-jones repulsive force and artificial viscosity.
        procedure:: sweep => fluid_boundary_sweep_monaghan1994
    end type fluid_boundary_sweeper_monaghan1994_t

    !> @brief Sweeper describing how to update weakly compressible virtual boundary particles' velocity and density.
    type, extends(base_sweeper_t):: boundary_update_sweeper_t
    contains
        !> @brief Calculate weakly compressible virtual boundary particles' velocity and density.
        procedure:: sweep => boundary_update_sweep
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
        !> @brief Generates ghost particles based on the mirroring boundary and the RHS particles.
        procedure:: sweep => ghost_timestep_setup_sweep
    end type

    !> @brief Sweeper describing the interaction between real (LHS) and virtual (RHS) weakly compressible particles using the morris
    !>        boundary condition.
    type, extends(fluid_sweeper_t):: morris_boundary_sweeper_t
        !> @brief The point that the boundary plane passes through.
        real(fp):: point(ndims)
        !> @brief The unit normal vector pointing towards the "real" domain.
        real(fp):: normal(ndims)
    contains
        !> @brief Calculates acceleration and density change of the LHS particles due to interaction with the morris boundary
        !>       particles (RHS).
        procedure:: sweep => morris_boundary_sweep
    end type morris_boundary_sweeper_t

    !> @brief Sweeper that calculates engineering strain rate between two particle systems.
    type, extends(base_sweeper_t):: strain_rate_sweeper_t
    contains
        !> @brief Calculates engineering strain rate.
        procedure:: sweep => strain_rate_sweep
    end type strain_rate_sweeper_t

    !> @brief Sweeper describing how to generate weakly compressible ghost particles with strain rate and cauchy stress tensors.
    type, extends(ghost_timestep_setuper_t):: eos_viscous_stress_ghost_timestep_setuper_t
    contains
        !> @brief Generates ghost particles based on the mirroring boundary and the RHS particles.
        procedure:: sweep => eos_viscous_stress_ghost_timestep_setup_sweep
    end type eos_viscous_stress_ghost_timestep_setuper_t

    !> @brief Sweeper describing interaction between systems of weakly compressible particles with cauchy stress tensors.
    type, extends(fluid_sweeper_t):: viscous_stress_fluid_sweeper_t
    contains
        !> @brief Calculate density change and acceleration between weakly compressible particles with cauchy stress tensor.
        procedure:: sweep => viscous_stress_fluid_sweep
    end type viscous_stress_fluid_sweeper_t

    !> @brief Sweeper describing the interaction between real (LHS) and virtual (RHS) weakly compressible articles (with cauchy
    !>        stress tensor) using the morris boundary condition.
    type, extends(viscous_stress_fluid_sweeper_t):: eos_viscous_stress_morris_boundary_sweeper_t
        !> @brief The point that the boundary plane passes through.
        real(fp):: point(ndims)
        !> @brief The unit normal vector pointing towards the "real" domain.
        real(fp):: normal(ndims)
    contains
        !> @brief Calculates acceleration and density change of the LHS particles due to interaction with the weakly compressible
        !>       morris boundary particles with stress tensor (RHS).
        procedure:: sweep => viscous_stress_morris_boundary_sweep
    end type eos_viscous_stress_morris_boundary_sweeper_t

    !> @brief Sweeper that describes how morris boundary particles contribute to the strain rate calculation of real particles.
    type, extends(strain_rate_sweeper_t):: strain_rate_morris_boundary_sweeper_t
        !> @brief The point that the boundary plane passes through.
        real(fp):: point(ndims)
        !> @brief The unit normal vector pointing towards the "real" domain.
        real(fp):: normal(ndims)
    contains
        !> @brief Calculates Contribution of morris boundary particles (RHS) to real particles (LHS) strain rate.
        procedure:: sweep => strain_rate_morris_boundary_sweep
    end type strain_rate_morris_boundary_sweeper_t

    public:: fluid_sweeper_t, fluid_boundary_sweeper_monaghan1994_t, boundary_update_sweeper_t, ghost_timestep_setuper_t, &
             morris_boundary_sweeper_t, strain_rate_sweeper_t, viscous_stress_fluid_sweeper_t, &
             eos_viscous_stress_ghost_timestep_setuper_t, eos_viscous_stress_morris_boundary_sweeper_t, &
             strain_rate_morris_boundary_sweeper_t

contains

    !> @brief For either one or two weakly-compressible fluid particle systems, calculate acceleration and density rate-of-change
    !>        due to isotropic pressure, artificial viscosity, and mass continuity.
    !> @param self The sweeper class holding artificial viscosity constants and gravity.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS weakly compressible particles involved in the interactions.
    !> @param psys_rhs Ths RHS "                                                        ".
    subroutine fluid_sweep(self, pairs, psys_lhs, psys_rhs, dt)
        class(fluid_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        real(fp), optional, intent(in):: dt
        class(eos_particle_t), pointer:: fluid_lhs(:), fluid_rhs(:)
        integer:: i, j, k
        real(fp):: dummy_drhodt, dummy_dvxdt(ndims) ! dummy variables for when update_rhs is .false.

        ! point to lhs particlse for access to pressure member
        select type (ps => psys_lhs%particles)
        class is (eos_particle_t)
            fluid_lhs => ps
        class default
            error stop "Invalid type for psys_lhs"
        end select

        ! intialize LHS acceleration and density rate-of-change arrays
        if (self%initialize) then
            do i = 1, size(fluid_lhs)
                fluid_lhs(i)%dvxdt(:) = 0._fp
                fluid_lhs(i)%dvxdt(ndims) = self%g
                fluid_lhs(i)%drhodt = 0._fp
            end do
        end if

        ! branch to handle logic for when psys_rhs is present as well as whether to update rhs
        if (present(psys_rhs)) then
            ! point to rhs particlse for access to pressure member
            select type (ps => psys_rhs%particles)
            class is (eos_particle_t)
                fluid_rhs => ps
            class default
                error stop "Invalid type for psys_rhs"
            end select
            if (self%update_rhs) then ! sweep using both lhs and rhs, and updating both
                ! intialize RHS acceleration and density rate-of-change arrays
                if (self%initialize) then
                    do i = 1, size(fluid_rhs)
                        fluid_rhs(i)%dvxdt(:) = 0._fp
                        fluid_rhs(i)%dvxdt(ndims) = self%g
                        fluid_rhs(i)%drhodt = 0._fp
                    end do
                end if

                ! perform sweep
                do k = 1, pairs%npairs_total
                    i = pairs%pair_ij(1, k)
                    j = pairs%pair_ij(2, k)
                    call artificial_viscosity_monaghan1994( &
                        fluid_lhs(i)%x(:), fluid_rhs(j)%x(:), fluid_lhs(i)%v(:), fluid_rhs(j)%v(:), fluid_lhs(i)%rho, &
                        fluid_rhs(j)%rho, self%h, self%h, fluid_lhs(i)%c, fluid_rhs(j)%c, fluid_lhs(i)%mass, fluid_rhs(j)%mass, &
                        fluid_lhs(i)%dvxdt(:), fluid_rhs(j)%dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                        )
                    call isotropic_pressure_force( &
                        fluid_lhs(i)%p, fluid_rhs(j)%p, fluid_lhs(i)%rho, fluid_rhs(j)%rho, fluid_lhs(i)%mass, fluid_rhs(j)%mass, &
                        fluid_lhs(i)%dvxdt(:), fluid_rhs(j)%dvxdt(:), pairs%dwdx(:, k) &
                        )
                    call continuity_density( &
                        fluid_lhs(i)%v(:), fluid_rhs(j)%v(:), fluid_lhs(i)%mass, fluid_rhs(j)%mass, fluid_lhs(i)%drhodt, &
                        fluid_rhs(j)%drhodt, pairs%dwdx(:, k) &
                        )
                end do
            else ! sweep using both lhs and rhs, but updating only lhs

                ! perform sweep
                do k = 1, pairs%npairs_total
                    i = pairs%pair_ij(1, k)
                    j = pairs%pair_ij(2, k)
                    call artificial_viscosity_monaghan1994( &
                        fluid_lhs(i)%x(:), fluid_rhs(j)%x(:), fluid_lhs(i)%v(:), fluid_rhs(j)%v(:), fluid_lhs(i)%rho, &
                        fluid_rhs(j)%rho, self%h, self%h, fluid_lhs(i)%c, fluid_rhs(j)%c, fluid_lhs(i)%mass, fluid_rhs(j)%mass, &
                        fluid_lhs(i)%dvxdt(:), dummy_dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                        )
                    call isotropic_pressure_force( &
                        fluid_lhs(i)%p, fluid_rhs(j)%p, fluid_lhs(i)%rho, fluid_rhs(j)%rho, fluid_lhs(i)%mass, fluid_rhs(j)%mass, &
                        fluid_lhs(i)%dvxdt(:), dummy_dvxdt(:), pairs%dwdx(:, k) &
                        )
                    call continuity_density( &
                        fluid_lhs(i)%v(:), fluid_rhs(j)%v(:), fluid_lhs(i)%mass, fluid_rhs(j)%mass, &
                        fluid_lhs(i)%drhodt, dummy_drhodt, pairs%dwdx(:, k) &
                        )
                end do

            end if

        else ! self-sweep using only psys_lhs

            ! perform sweep
            do k = 1, pairs%npairs_total
                i = pairs%pair_ij(1, k)
                j = pairs%pair_ij(2, k)
                call artificial_viscosity_monaghan1994( &
                    fluid_lhs(i)%x(:), fluid_lhs(j)%x(:), fluid_lhs(i)%v(:), fluid_lhs(j)%v(:), fluid_lhs(i)%rho, &
                    fluid_lhs(j)%rho, self%h, self%h, fluid_lhs(i)%c, fluid_lhs(j)%c, fluid_lhs(i)%mass, fluid_lhs(j)%mass, &
                    fluid_lhs(i)%dvxdt(:), fluid_lhs(j)%dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                    )
                call isotropic_pressure_force( &
                    fluid_lhs(i)%p, fluid_lhs(j)%p, fluid_lhs(i)%rho, fluid_lhs(j)%rho, fluid_lhs(i)%mass, &
                    fluid_lhs(j)%mass, fluid_lhs(i)%dvxdt(:), fluid_lhs(j)%dvxdt(:), pairs%dwdx(:, k) &
                    )
                call continuity_density( &
                    fluid_lhs(i)%v(:), fluid_lhs(j)%v(:), fluid_lhs(i)%mass, fluid_lhs(j)%mass, &
                    fluid_lhs(i)%drhodt, fluid_lhs(j)%drhodt, pairs%dwdx(:, k) &
                    )
            end do

        end if
    end subroutine fluid_sweep

    !> @brief For a weakly-compressible fluid particle system and a base particle system, calculate acceleration due to
    !>        lennard-jones repulsive force and artificial viscosity.
    !> @param self The sweeper class holding artificial viscosity constants.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the real particles involved in the interactions.
    !> @param psys_rhs The boundary particles involved in the interactions.
    subroutine fluid_boundary_sweep_monaghan1994(self, pairs, psys_lhs, psys_rhs, dt)
        class(fluid_boundary_sweeper_monaghan1994_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        real(fp):: dummy_dvxdt(2)

        if (.not. present(psys_rhs)) error stop "expected psys_rhs to be passed in."

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            ! apply boundary force with eqn 4.1.
            call repulsive_force( &
                self%cutoff, psys_lhs%particles(i)%c, psys_lhs%particles(i)%x(:), psys_rhs%particles(j)%x(:), &
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

    end subroutine fluid_boundary_sweep_monaghan1994

    !> @brief For two particle systems, calculate the LHS particles' density and velocity using kernel interpolation over the RHS
    !>        particles' density and velocity.
    !> @param self The sweeper class.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS particles whose velocity and density will be updated.
    !> @param psys_rhs Ths RHS particles to calculate velocity and density from.
    subroutine boundary_update_sweep(self, pairs, psys_lhs, psys_rhs, dt)
        class(boundary_update_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        real(fp):: mw, vw
        real(fp), allocatable:: wsum(:)

        if (present(psys_rhs)) then
            allocate (wsum(psys_rhs%size), source=0._fp)
        else
            error stop "Expected psys_rhs to be passed in."
        end if

        do i = 1, psys_rhs%size
            psys_rhs%particles(i)%rho = 0._fp
            psys_rhs%particles(i)%v(:) = 0._fp
        end do

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            mw = psys_lhs%particles(i)%mass*pairs%w(k)
            vw = mw/psys_lhs%particles(i)%rho
            wsum(j) = wsum(j) + vw
            psys_rhs%particles(j)%rho = psys_rhs%particles(j)%rho + mw
            psys_rhs%particles(j)%v(:) = psys_rhs%particles(j)%v(:) + psys_lhs%particles(i)%v(:)*vw
        end do

        do j = 1, psys_rhs%size
            if (wsum(j) > 0._fp) then
                psys_rhs%particles(j)%v(:) = -psys_rhs%particles(j)%v(:)/wsum(j)
                psys_rhs%particles(j)%rho = psys_rhs%particles(j)%rho/wsum(j)
            end if
        end do

    end subroutine boundary_update_sweep

    !> @brief For two weakly compressible particle systems, generate ghost particles in the LHS system, using information from the
    !>        RHS system and a defined boundary.
    !> @param self The sweeper class containing the plane boundary using a unit normal vector and point.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS ghost particle system where ghost particles will be generated.
    !> @param psys_rhs Ths RHS particles to generate ghost particles from.
    subroutine ghost_timestep_setup_sweep(self, pairs, psys_lhs, psys_rhs, dt)
        class(ghost_timestep_setuper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        real(fp), optional, intent(in):: dt
        integer:: i
        class(eos_particle_t), pointer:: ps_real(:)
        class(eos_ghost_particle_t), pointer:: ps_ghost(:)
        real(fp):: dx(ndims), dr

        if (.not. present(psys_rhs)) error stop "Expected psys_rhs to be passed in."

        select type (ps => psys_lhs%particles)
        class is (eos_particle_t)
            ps_real => ps
        class default
            error stop "Expected psys_lhs to be eos_particle_t."
        end select

        select type (ps => psys_rhs%particles)
        class is (eos_ghost_particle_t)
            ps_ghost => ps
        class default
            error stop "Expected psys_rhs to be eos_ghost_particle_t."
        end select

        psys_rhs%size = 0

        do i = 1, psys_lhs%size
            dr = dot_product(ps_real(i)%x(:) - self%point(:), self%surface_normal(:))
            ! ghost particle if within cutoff and inside the modelled region.
            if (abs(dr) <= self%cutoff .and. dr > 0._fp) then
                psys_rhs%size = psys_rhs%size + 1
                ps_ghost(psys_rhs%size)%id = psys_rhs%size
                ps_ghost(psys_rhs%size)%original => ps_real(i)
                ps_ghost(psys_rhs%size)%x(:) = ps_real(i)%x(:) - 2._fp*dr*self%surface_normal(:)
            end if

        end do

    end subroutine ghost_timestep_setup_sweep

    !> @brief For a weakly compressible particle system and base particle system (aka boundary particles), calculate the
    !>        contribution of the boundary particles to the acceleration and rate of change to the real particles using the morris
    !>        boundary.
    !> @param self The sweeper class holding artificial viscosity constants and boundary surface point and unit normal vector.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS weakly compressible particles involved in the interactions.
    !> @param psys_rhs Ths RHS boundary particles involved in the interactions.
    subroutine morris_boundary_sweep(self, pairs, psys_lhs, psys_rhs, dt)
        class(morris_boundary_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        class(eos_particle_t), pointer:: ps_fluid(:)
        real(fp):: dummy_drhodt, dummy_dvxdt(ndims), vb(ndims), da, db

        if (.not. present(psys_rhs)) then
            error stop "Expected psys_rhs to be passed in."
        end if

        select type (ps => psys_lhs%particles)
        class is (eos_particle_t)
            ps_fluid => ps
        class default
            error stop "Expected psys_lhs%particles to be eos_particle_t."
        end select

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            da = dot_product(ps_fluid(i)%x(:) - self%point(:), self%normal(:))
            db = dot_product(psys_rhs%particles(j)%x(:) - self%point(:), self%normal(:))
            vb(:) = -min(3._fp, abs(db/da))*ps_fluid(i)%v(:)
            call artificial_viscosity_monaghan1994( &
                ps_fluid(i)%x(:), psys_rhs%particles(j)%x(:), ps_fluid(i)%v(:), vb(:), ps_fluid(i)%rho, &
                ps_fluid(i)%rho, self%h, self%h, ps_fluid(i)%c, ps_fluid(i)%c, ps_fluid(i)%mass, ps_fluid(i)%mass, &
                ps_fluid(i)%dvxdt(:), dummy_dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                )
            call isotropic_pressure_force( &
                ps_fluid(i)%p, ps_fluid(i)%p, ps_fluid(i)%rho, ps_fluid(i)%rho, ps_fluid(i)%mass, ps_fluid(i)%mass, &
                ps_fluid(i)%dvxdt(:), dummy_dvxdt(:), pairs%dwdx(:, k) &
                )
            call continuity_density( &
                ps_fluid(i)%v(:), vb(:), ps_fluid(i)%mass, ps_fluid(i)%mass, &
                ps_fluid(i)%drhodt, dummy_drhodt, pairs%dwdx(:, k) &
                )
        end do

    end subroutine morris_boundary_sweep

    !> @brief For either one or two particle systems, calculate strain rate.
    !> @param self The sweeper class.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS particles involved in the interactions.
    !> @param psys_rhs Ths RHS "                                    ".
    subroutine strain_rate_sweep(self, pairs, psys_lhs, psys_rhs, dt)

        use weakly_compressible_particles_m, only: ntensor_elems_voigt

        class(strain_rate_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        real(fp), optional, intent(in):: dt
        class(eos_viscous_stress_particle_t), pointer:: ps_lhs(:), ps_rhs(:)
        real(fp):: dummy_strain_rate(ntensor_elems_voigt)
        integer:: i, j, k

        select type (ps => psys_lhs%particles)
        class is (eos_viscous_stress_particle_t)
            ps_lhs => ps
        class default
            error stop "Expected psys_lhs%particles class to be eos_viscous_stress_particle_t."
        end select

        if (self%initialize) then
            do i = 1, psys_lhs%size
                ps_lhs(i)%strain_rate(:) = 0._fp
            end do
        end if

        if (present(psys_rhs)) then

            select type (ps => psys_rhs%particles)
            class is (eos_viscous_stress_particle_t)
                ps_rhs => ps
            class default
                error stop "Expected psys_rhs%particles class to be eos_viscous_stress_particle_t."
            end select
            if (self%update_rhs) then
                do k = 1, pairs%npairs_total
                    i = pairs%pair_ij(1, k)
                    j = pairs%pair_ij(2, k)
                    call strain_rate( &
                        ps_lhs(i)%v(:), ps_rhs(j)%v(:), ps_lhs(i)%mass, ps_rhs(j)%mass, ps_lhs(i)%rho, ps_rhs(j)%rho, &
                        pairs%dwdx(:, k), ps_lhs(i)%strain_rate(:), ps_rhs(j)%strain_rate(:) &
                        )
                end do
            else
                do k = 1, pairs%npairs_total
                    i = pairs%pair_ij(1, k)
                    j = pairs%pair_ij(2, k)
                    call strain_rate( &
                        ps_lhs(i)%v(:), ps_rhs(j)%v(:), ps_lhs(i)%mass, ps_rhs(j)%mass, ps_lhs(i)%rho, ps_rhs(j)%rho, &
                        pairs%dwdx(:, k), ps_lhs(i)%strain_rate(:), dummy_strain_rate(:) &
                        )
                end do
            end if
        else
            do k = 1, pairs%npairs_total
                i = pairs%pair_ij(1, k)
                j = pairs%pair_ij(2, k)
                call strain_rate( &
                    ps_lhs(i)%v(:), ps_lhs(j)%v(:), ps_lhs(i)%mass, ps_lhs(j)%mass, ps_lhs(i)%rho, ps_lhs(j)%rho, &
                    pairs%dwdx(:, k), ps_lhs(i)%strain_rate(:), ps_lhs(j)%strain_rate(:) &
                    )
            end do
        end if

    end subroutine strain_rate_sweep

    !> @brief For two weakly compressible particle systems with stress tensor, generate ghost particles in the LHS system, using
    !>        information from the RHS system and a defined boundary.
    !> @param self The sweeper class containing the plane boundary using a unit normal vector and point.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS ghost particle system where ghost particles will be generated.
    !> @param psys_rhs Ths RHS particles to generate ghost particles from.
    subroutine eos_viscous_stress_ghost_timestep_setup_sweep(self, pairs, psys_lhs, psys_rhs, dt)
        class(eos_viscous_stress_ghost_timestep_setuper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        real(fp), optional, intent(in):: dt
        integer:: i
        class(eos_viscous_stress_particle_t), pointer:: ps_real(:)
        class(eos_viscous_stress_ghost_particle_t), pointer:: ps_ghost(:)
        real(fp):: dr

        if (.not. present(psys_rhs)) error stop "Expected psys_rhs to be passed in."

        select type (ps => psys_lhs%particles)
        class is (eos_viscous_stress_particle_t)
            ps_real => ps
        class default
            error stop "Expected psys_lhs to be eos_viscous_stress_particle_t."
        end select

        select type (ps => psys_rhs%particles)
        class is (eos_viscous_stress_ghost_particle_t)
            ps_ghost => ps
        class default
            error stop "Expected psys_rhs to be eos_viscous_stress_ghost_particle_t."
        end select

        psys_rhs%size = 0

        do i = 1, psys_lhs%size
            dr = dot_product(ps_real(i)%x(:) - self%point(:), self%surface_normal(:))
            ! ghost particle if within cutoff and inside the modelled region.
            if (abs(dr) <= self%cutoff .and. dr > 0._fp) then
                psys_rhs%size = psys_rhs%size + 1
                ps_ghost(psys_rhs%size)%id = psys_rhs%size
                ps_ghost(psys_rhs%size)%original => ps_real(i)
                ps_ghost(psys_rhs%size)%x(:) = ps_real(i)%x(:) - 2._fp*dr*self%surface_normal(:)
            end if

        end do

    end subroutine eos_viscous_stress_ghost_timestep_setup_sweep

    !> @brief For either one or two weakly-compressible fluid particle systems with stress ensor, calculate acceleration and density
    !>        rate-of-change due to cauchy stress, artificial viscosity, and mass continuity.
    !> @param self The sweeper class holding artificial viscosity constants and gravity.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS weakly compressible particles with stress tensor involved in the interactions.
    !> @param psys_rhs Ths RHS "                                                                           ".
    subroutine viscous_stress_fluid_sweep(self, pairs, psys_lhs, psys_rhs, dt)
        class(viscous_stress_fluid_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        class(eos_viscous_stress_particle_t), pointer:: fluid_lhs(:), fluid_rhs(:)
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        real(fp):: dummy_drhodt, dummy_dvxdt(ndims) ! dummy variables for when update_rhs is .false.

        ! point to lhs particlse for access to pressure member
        select type (ps => psys_lhs%particles)
        class is (eos_viscous_stress_particle_t)
            fluid_lhs => ps
        class default
            error stop "Invalid type for psys_lhs"
        end select

        ! intialize LHS acceleration and density rate-of-change arrays
        if (self%initialize) then
            do i = 1, psys_lhs%size
                fluid_lhs(i)%dvxdt(:) = 0._fp
                fluid_lhs(i)%dvxdt(ndims) = self%g
                fluid_lhs(i)%drhodt = 0._fp
            end do
        end if

        ! branch to handle logic for when psys_rhs is present as well as whether to update rhs
        if (present(psys_rhs)) then
            ! point to rhs particlse for access to pressure member
            select type (ps => psys_rhs%particles)
            class is (eos_viscous_stress_particle_t)
                fluid_rhs => ps
            class default
                error stop "Invalid type for psys_rhs"
            end select
            if (self%update_rhs) then ! sweep using both lhs and rhs, and updating both
                ! intialize RHS acceleration and density rate-of-change arrays
                if (self%initialize) then
                    do i = 1, size(fluid_rhs)
                        fluid_rhs(i)%dvxdt(:) = 0._fp
                        fluid_rhs(i)%dvxdt(ndims) = self%g
                        fluid_rhs(i)%drhodt = 0._fp
                    end do
                end if

                ! perform sweep
                do k = 1, pairs%npairs_total
                    i = pairs%pair_ij(1, k)
                    j = pairs%pair_ij(2, k)
                    call artificial_viscosity_monaghan1994( &
                        fluid_lhs(i)%x(:), fluid_rhs(j)%x(:), fluid_lhs(i)%v(:), fluid_rhs(j)%v(:), fluid_lhs(i)%rho, &
                        fluid_rhs(j)%rho, self%h, self%h, fluid_lhs(i)%c, fluid_rhs(j)%c, fluid_lhs(i)%mass, fluid_rhs(j)%mass, &
                        fluid_lhs(i)%dvxdt(:), fluid_rhs(j)%dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                        )
                    call cauchy_stress_force( &
                        fluid_lhs(i)%stress, fluid_rhs(j)%stress, fluid_lhs(i)%rho, fluid_rhs(j)%rho, fluid_lhs(i)%mass, &
                        fluid_rhs(j)%mass, fluid_lhs(i)%dvxdt(:), fluid_rhs(j)%dvxdt(:), pairs%dwdx(:, k) &
                        )
                    call continuity_density( &
                        fluid_lhs(i)%v(:), fluid_rhs(j)%v(:), fluid_lhs(i)%mass, fluid_rhs(j)%mass, fluid_lhs(i)%drhodt, &
                        fluid_rhs(j)%drhodt, pairs%dwdx(:, k) &
                        )
                    call diffusion_density( &
                        fluid_lhs(i)%rho, fluid_rhs(j)%rho, fluid_lhs(i)%x(:), fluid_rhs(j)%x(:), fluid_lhs(i)%mass, &
                        fluid_rhs(j)%mass, self%h, self%h, fluid_lhs(i)%c, fluid_rhs(j)%c, pairs%dwdx(:, k), fluid_lhs(i)%drhodt, &
                        fluid_rhs(j)%drhodt &
                        )
                end do
            else ! sweep using both lhs and rhs, but updating only lhs
                ! perform sweep
                do k = 1, pairs%npairs_total
                    i = pairs%pair_ij(1, k)
                    j = pairs%pair_ij(2, k)
                    call artificial_viscosity_monaghan1994( &
                        fluid_lhs(i)%x(:), fluid_rhs(j)%x(:), fluid_lhs(i)%v(:), fluid_rhs(j)%v(:), fluid_lhs(i)%rho, &
                        fluid_rhs(j)%rho, self%h, self%h, fluid_lhs(i)%c, fluid_rhs(j)%c, fluid_lhs(i)%mass, fluid_rhs(j)%mass, &
                        fluid_lhs(i)%dvxdt(:), dummy_dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                        )
                    call cauchy_stress_force( &
                        fluid_lhs(i)%stress, fluid_rhs(j)%stress, fluid_lhs(i)%rho, fluid_rhs(j)%rho, fluid_lhs(i)%mass, &
                        fluid_rhs(j)%mass, fluid_lhs(i)%dvxdt(:), dummy_dvxdt(:), pairs%dwdx(:, k) &
                        )
                    call continuity_density( &
                        fluid_lhs(i)%v(:), fluid_rhs(j)%v(:), fluid_lhs(i)%mass, fluid_rhs(j)%mass, &
                        fluid_lhs(i)%drhodt, dummy_drhodt, pairs%dwdx(:, k) &
                        )
                    call diffusion_density( &
                        fluid_lhs(i)%rho, fluid_rhs(j)%rho, fluid_lhs(i)%x(:), fluid_rhs(j)%x(:), fluid_lhs(i)%mass, &
                        fluid_rhs(j)%mass, self%h, self%h, fluid_lhs(i)%c, fluid_rhs(j)%c, pairs%dwdx(:, k), fluid_lhs(i)%drhodt, &
                        dummy_drhodt &
                        )
                end do
            end if

        else ! self-sweep using only psys_lhs

            ! perform sweep
            do k = 1, pairs%npairs_total
                i = pairs%pair_ij(1, k)
                j = pairs%pair_ij(2, k)
                call artificial_viscosity_monaghan1994( &
                    fluid_lhs(i)%x(:), fluid_lhs(j)%x(:), fluid_lhs(i)%v(:), fluid_lhs(j)%v(:), fluid_lhs(i)%rho, &
                    fluid_lhs(j)%rho, self%h, self%h, fluid_lhs(i)%c, fluid_lhs(j)%c, fluid_lhs(i)%mass, fluid_lhs(j)%mass, &
                    fluid_lhs(i)%dvxdt(:), fluid_lhs(j)%dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                    )
                call cauchy_stress_force( &
                    fluid_lhs(i)%stress, fluid_lhs(j)%stress, fluid_lhs(i)%rho, fluid_lhs(j)%rho, fluid_lhs(i)%mass, &
                    fluid_lhs(j)%mass, fluid_lhs(i)%dvxdt(:), fluid_lhs(j)%dvxdt(:), pairs%dwdx(:, k) &
                    )
                call continuity_density( &
                    fluid_lhs(i)%v(:), fluid_lhs(j)%v(:), fluid_lhs(i)%mass, fluid_lhs(j)%mass, &
                    fluid_lhs(i)%drhodt, fluid_lhs(j)%drhodt, pairs%dwdx(:, k) &
                    )
                call diffusion_density( &
                    fluid_lhs(i)%rho, fluid_lhs(j)%rho, fluid_lhs(i)%x(:), fluid_lhs(j)%x(:), fluid_lhs(i)%mass, &
                    fluid_lhs(j)%mass, self%h, self%h, fluid_lhs(i)%c, fluid_lhs(j)%c, pairs%dwdx(:, k), fluid_lhs(i)%drhodt, &
                    fluid_lhs(j)%drhodt &
                    )
            end do

        end if
    end subroutine viscous_stress_fluid_sweep

    !> @brief For a weakly compressible particle system with stress tensor and base particle system (aka boundary particles),
    !>        calculate the contribution of the boundary particles to the acceleration and rate of change to the real particles
    !>        using the morris boundary.
    !> @param self The sweeper class holding artificial viscosity constants and boundary surface point and unit normal vector.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS weakly compressible particles with stress tensor involved in the interactions.
    !> @param psys_rhs Ths RHS boundary particles involved in the interactions.
    subroutine viscous_stress_morris_boundary_sweep(self, pairs, psys_lhs, psys_rhs, dt)
        class(eos_viscous_stress_morris_boundary_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        class(eos_viscous_stress_particle_t), pointer:: ps_fluid(:)
        real(fp):: dummy_drhodt, dummy_dvxdt(ndims), vb(ndims), da, db

        if (.not. present(psys_rhs)) then
            error stop "Expected psys_rhs to be passed in."
        end if

        select type (ps => psys_lhs%particles)
        class is (eos_viscous_stress_particle_t)
            ps_fluid => ps
        class default
            error stop "Expected psys_lhs%particles to be eos_viscous_stress_particle_t."
        end select

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            da = dot_product(ps_fluid(i)%x(:) - self%point(:), self%normal(:))
            db = dot_product(psys_rhs%particles(j)%x(:) - self%point(:), self%normal(:))
            vb(:) = -min(3._fp, abs(db/da))*ps_fluid(i)%v(:)
            call artificial_viscosity_monaghan1994( &
                ps_fluid(i)%x(:), psys_rhs%particles(j)%x(:), ps_fluid(i)%v(:), vb(:), ps_fluid(i)%rho, &
                ps_fluid(i)%rho, self%h, self%h, ps_fluid(i)%c, ps_fluid(i)%c, ps_fluid(i)%mass, ps_fluid(i)%mass, &
                ps_fluid(i)%dvxdt(:), dummy_dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                )
            call cauchy_stress_force( &
                ps_fluid(i)%stress, ps_fluid(i)%stress, ps_fluid(i)%rho, ps_fluid(i)%rho, ps_fluid(i)%mass, ps_fluid(i)%mass, &
                ps_fluid(i)%dvxdt(:), dummy_dvxdt(:), pairs%dwdx(:, k) &
                )
            call continuity_density( &
                ps_fluid(i)%v(:), vb(:), ps_fluid(i)%mass, ps_fluid(i)%mass, &
                ps_fluid(i)%drhodt, dummy_drhodt, pairs%dwdx(:, k) &
                )
        end do

    end subroutine viscous_stress_morris_boundary_sweep

    !> @brief For two base particle systems, calculate the contribution of the boundary particles (RHS) to the strain rate of the
    !>        real particles (LHS) using the morris boundary.
    !> @param self The sweeper class holding boundary surface point and unit normal vector.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS particles involved in the interactions.
    !> @param psys_rhs Ths RHS "                                    ".
    subroutine strain_rate_morris_boundary_sweep(self, pairs, psys_lhs, psys_rhs, dt)

        use weakly_compressible_particles_m, only: ntensor_elems_voigt

        class(strain_rate_morris_boundary_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        real(fp), optional, intent(in):: dt
        class(eos_viscous_stress_particle_t), pointer:: ps_lhs(:)
        real(fp):: da, db, vb(ndims), dummy_strain_rate(ntensor_elems_voigt)
        integer:: i, j, k

        select type (ps => psys_lhs%particles)
        class is (eos_viscous_stress_particle_t)
            ps_lhs => ps
        class default
            error stop "Expected psys_lhs%particles class to be eos_viscous_stress_particle_t."
        end select

        if (.not. present(psys_rhs)) error stop "Expected psys_rhs to be passed in."

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            da = dot_product(ps_lhs(i)%x(:) - self%point(:), self%normal(:))
            db = dot_product(psys_rhs%particles(j)%x(:) - self%point(:), self%normal(:))
            vb(:) = -min(3._fp, abs(db/da))*ps_lhs(i)%v(:)
            call strain_rate( &
                ps_lhs(i)%v(:), vb(:), ps_lhs(i)%mass, ps_lhs(i)%mass, ps_lhs(i)%rho, ps_lhs(i)%rho, &
                pairs%dwdx(:, k), ps_lhs(i)%strain_rate(:), dummy_strain_rate(:) &
                )
        end do

    end subroutine strain_rate_morris_boundary_sweep

end module weakly_compressible_interactions_m
