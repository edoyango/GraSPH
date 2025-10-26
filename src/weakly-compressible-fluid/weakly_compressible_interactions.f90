!> @file weakly_compressible_interactions.f90
!> @brief Module containing subroutines for describing interactions between weakly compressible particles
!> @author Edward Yang
!> @date 2025-09-21
module weakly_compressible_interactions_m

    use grasph_constants_m, only: fp, ndims
    use grasph_particle_system_m, only: particle_system_t
    use weakly_compressible_particles_m, only: eos_particle_t, eos_ghost_particle_t
    use grasph_system_interactions_m, only: base_sweeper_t
    use grasph_pairs_m, only: particle_pairs_t
    use grasph_pair_interactions_m, only: artificial_viscosity_monaghan1994, continuity_density, isotropic_pressure_force

    implicit none

    private

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

    type, extends(base_sweeper_t):: ghost_timestep_setuper_t
        real(fp):: cutoff
        real(fp):: surface_normal(ndims), point(ndims)
    contains
        procedure:: sweep => ghost_timestep_setup_sweep
    end type

    ! define how fluid particles interact with boundary
    type, extends(fluid_sweeper_t):: morris_boundary_sweeper_t
        real(fp):: point(ndims), normal(ndims)
    contains
        procedure:: sweep => morris_boundary_sweep
    end type morris_boundary_sweeper_t

    public:: fluid_sweeper_t, ghost_timestep_setuper_t, morris_boundary_sweeper_t

contains

    !> @brief For a single set of weakly-compressible fluid particles, calculate acceleration and density rate-of-change due to
    !>        isotropic pressure, artificial viscosity, and mass continuity.
    !> @param self The sweeper class holding artificial viscosity constants and gravity.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS particles involved in the interactions.
    !> @param psys_rhs Ths RHS particles involved in the interactions. Expecting not to be passed in.
    subroutine fluid_sweep(self, pairs, psys_lhs, psys_rhs)
        class(fluid_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
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

    subroutine ghost_timestep_setup_sweep(self, pairs, psys_lhs, psys_rhs)
        class(ghost_timestep_setuper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
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

    subroutine morris_boundary_sweep(self, pairs, psys_lhs, psys_rhs)
        class(morris_boundary_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
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

end module weakly_compressible_interactions_m
