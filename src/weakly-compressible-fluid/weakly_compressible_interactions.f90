!> @file weakly_compressible_interactions.f90
!> @brief Module containing subroutines for describing interactions between weakly compressible particles
!> @author Edward Yang
!> @date 2025-09-21
module weakly_compressible_interactions

    use grasph_constants, only: fp, ndims
    use grasph_particles, only: base_particles
    use weakly_compressible_particles, only: linear_eos_particles, linear_eos_particle
    use grasph_pair_sets, only: particle_interactions, base_sweeper
    use grasph_pairs, only: particle_pairs
    use grasph_pair_interactions, only: artificial_viscosity_monaghan1994, continuity_density, isotropic_pressure_force

    implicit none

    type, extends(base_sweeper):: fluid_sweeper
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
    end type fluid_sweeper

contains

    !> @brief For a single set of weakly-compressible fluid particles, calculate acceleration and density rate-of-change due to
    !>        isotropic pressure, artificial viscosity, and mass continuity.
    !> @param self The sweeper class holding artificial viscosity constants and gravity.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS particles involved in the interactions.
    !> @param ps_rhs Ths RHS particles involved in the interactions. Expecting not to be passed in.
    subroutine fluid_sweep(self, pairs, ps_lhs, ps_rhs)
        class(fluid_sweeper), intent(in):: self
        type(particle_pairs), intent(in):: pairs
        class(base_particles), intent(inout):: ps_lhs
        class(base_particles), optional, intent(inout):: ps_rhs
        class(linear_eos_particle), pointer:: fluid_lhs(:), fluid_rhs(:)
        integer:: i, j, k
        real(fp):: dummy_drhodt, dummy_dvxdt(ndims) ! dummy variables for when update_rhs is .false.

        ! point to lhs particlse for access to pressure member
        select type (ps => ps_lhs%ps)
        class is (linear_eos_particle)
            fluid_lhs => ps
        class default
            error stop "Invalid type for ps_lhs"
        end select

        ! intialize LHS acceleration and density rate-of-change arrays
        if (self%initialize) then
            do i = 1, size(fluid_lhs)
                fluid_lhs(i)%dvxdt(:) = 0._fp
                fluid_lhs(i)%dvxdt(ndims) = self%g
                fluid_lhs(i)%drhodt = 0._fp
            end do
        end if

        ! branch to handle logic for when ps_rhs is present as well as whether to update rhs
        if (present(ps_rhs)) then
            ! point to rhs particlse for access to pressure member
            select type (ps => ps_rhs%ps)
            class is (linear_eos_particle)
                fluid_rhs => ps
            class default
                error stop "Invalid type for ps_rhs"
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

        else ! self-sweep using only ps_lhs

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

end module weakly_compressible_interactions
