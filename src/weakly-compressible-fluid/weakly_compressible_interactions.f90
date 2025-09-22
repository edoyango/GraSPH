!> @file weakly_compressible_interactions.f90
!> @brief Module containing subroutines for describing interactions between weakly compressible particles
!> @author Edward Yang
!> @date 2025-09-21
module weakly_compressible_interactions

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles
    use weakly_compressible_particles, only: linear_eos_particles
    use grasph_pair_sets, only: particle_interactions_base, base_sweeper
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

    !> @brief Describes how a single set of weakly-compressible fluid particles interact with itself.
    type, extends(particle_interactions_base):: fluid_self_interaction
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
        procedure:: sweep => fluid_self_sweep
    end type fluid_self_interaction

    type, extends(fluid_self_interaction):: fluid_fluid_interaction
    contains
        procedure:: sweep => fluid_fluid_sweep
    end type fluid_fluid_interaction

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
        class(linear_eos_particles), pointer:: fluid_lhs, fluid_rhs
        integer:: i, j, k
        real(fp):: dummy_drhodt, dummy_dvxdt(ps_lhs%ndims) ! dummy variables for when update_rhs is .false.

        ! point to lhs particlse for access to pressure member
        select type (ps => ps_lhs)
        class is (linear_eos_particles)
            fluid_lhs => ps
        class default
            error stop "Invalid type for ps_lhs"
        end select

        ! intialize LHS acceleration and density rate-of-change arrays
        do i = 1, fluid_lhs%size
            fluid_lhs%dvxdt(:, i) = 0._fp
            fluid_lhs%dvxdt(fluid_lhs%ndims, i) = self%g
            fluid_lhs%drhodt(i) = 0._fp
        end do

        ! branch to handle logic for when ps_rhs is present as well as whether to update rhs
        if (present(ps_rhs)) then
            ! point to rhs particlse for access to pressure member
            select type (ps => ps_rhs)
            class is (linear_eos_particles)
                fluid_rhs => ps
            class default
                error stop "Invalid type for ps_rhs"
            end select
            if (self%update_rhs) then ! sweep using both lhs and rhs, and updating both
                ! intialize RHS acceleration and density rate-of-change arrays
                do i = 1, fluid_rhs%size
                    fluid_rhs%dvxdt(:, i) = 0._fp
                    fluid_rhs%dvxdt(fluid_rhs%ndims, i) = self%g
                    fluid_rhs%drhodt(i) = 0._fp
                end do

                ! perform sweep
                do k = 1, pairs%npairs_total
                    i = pairs%pair_ij(1, k)
                    j = pairs%pair_ij(2, k)
                    call artificial_viscosity_monaghan1994( &
                        fluid_lhs%ndims, fluid_lhs%x(:, i), fluid_rhs%x(:, j), fluid_lhs%v(:, i), fluid_rhs%v(:, j), &
                        fluid_lhs%rho(i), fluid_rhs%rho(j), self%h, self%h, fluid_lhs%c(i), fluid_rhs%c(j), fluid_lhs%mass(i), &
                        fluid_rhs%mass(j), fluid_lhs%dvxdt(:, i), fluid_rhs%dvxdt(:, j), pairs%dwdx(:, k), self%artvisc_alpha, &
                        self%artvisc_beta &
                        )
                    call isotropic_pressure_force( &
                        fluid_lhs%ndims, fluid_lhs%p(i), fluid_rhs%p(j), fluid_lhs%rho(i), fluid_rhs%rho(j), fluid_lhs%mass(i), &
                        fluid_rhs%mass(j), fluid_lhs%dvxdt(:, i), fluid_rhs%dvxdt(:, j), pairs%dwdx(:, k) &
                        )
                    call continuity_density( &
                        fluid_lhs%ndims, fluid_lhs%v(:, i), fluid_rhs%v(:, j), fluid_lhs%mass(i), fluid_rhs%mass(j), &
                        fluid_lhs%drhodt(i), fluid_rhs%drhodt(j), pairs%dwdx(:, k) &
                        )
                end do
            else ! sweep using both lhs and rhs, but updating only lhs

                ! perform sweep
                do k = 1, pairs%npairs_total
                    i = pairs%pair_ij(1, k)
                    j = pairs%pair_ij(2, k)
                    call artificial_viscosity_monaghan1994( &
                        fluid_lhs%ndims, fluid_lhs%x(:, i), fluid_rhs%x(:, j), fluid_lhs%v(:, i), fluid_rhs%v(:, j), &
                        fluid_lhs%rho(i), fluid_rhs%rho(j), self%h, self%h, fluid_lhs%c(i), fluid_rhs%c(j), fluid_lhs%mass(i), &
                        fluid_rhs%mass(j), fluid_lhs%dvxdt(:, i), dummy_dvxdt(:), pairs%dwdx(:, k), self%artvisc_alpha, &
                        self%artvisc_beta &
                        )
                    call isotropic_pressure_force( &
                        fluid_lhs%ndims, fluid_lhs%p(i), fluid_rhs%p(j), fluid_lhs%rho(i), fluid_rhs%rho(j), fluid_lhs%mass(i), &
                        fluid_rhs%mass(j), fluid_lhs%dvxdt(:, i), dummy_dvxdt(:), pairs%dwdx(:, k) &
                        )
                    call continuity_density( &
                        fluid_lhs%ndims, fluid_lhs%v(:, i), fluid_rhs%v(:, j), fluid_lhs%mass(i), fluid_rhs%mass(j), &
                        fluid_lhs%drhodt(i), dummy_drhodt, pairs%dwdx(:, k) &
                        )
                end do

            end if

        else ! self-sweep using only ps_lhs

            ! perform sweep
            do k = 1, pairs%npairs_total
                i = pairs%pair_ij(1, k)
                j = pairs%pair_ij(2, k)
                call artificial_viscosity_monaghan1994( &
                    fluid_lhs%ndims, fluid_lhs%x(:, i), fluid_lhs%x(:, j), fluid_lhs%v(:, i), fluid_lhs%v(:, j), fluid_lhs%rho(i), &
                    fluid_lhs%rho(j), self%h, self%h, fluid_lhs%c(i), fluid_lhs%c(j), fluid_lhs%mass(i), fluid_lhs%mass(j), &
                    fluid_lhs%dvxdt(:, i), fluid_lhs%dvxdt(:, j), pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                    )
                call isotropic_pressure_force( &
                    fluid_lhs%ndims, fluid_lhs%p(i), fluid_lhs%p(j), fluid_lhs%rho(i), fluid_lhs%rho(j), fluid_lhs%mass(i), &
                    fluid_lhs%mass(j), fluid_lhs%dvxdt(:, i), fluid_lhs%dvxdt(:, j), pairs%dwdx(:, k) &
                    )
                call continuity_density( &
                    fluid_lhs%ndims, fluid_lhs%v(:, i), fluid_lhs%v(:, j), fluid_lhs%mass(i), fluid_lhs%mass(j), &
                    fluid_lhs%drhodt(i), fluid_lhs%drhodt(j), pairs%dwdx(:, k) &
                    )
            end do

        end if
    end subroutine fluid_sweep

    !> @brief For a single set of weakly-compressible fluid particles, calculate acceleration and density rate-of-change due to
    !>        isotropic pressure, artificial viscosity, and mass continuity. 2d case.
    !> @param self The particle interaction class to find particles within.
    subroutine fluid_self_sweep(self)
        class(fluid_self_interaction), intent(inout):: self
        integer:: i, j, k
        class(linear_eos_particles), pointer:: ps_fluid

        ! associate linear_eos_pointer if base pointer is associated
        select type (ps => self%ps_lhs)
        class is (linear_eos_particles)
            ps_fluid => ps
        class default
            error stop 'Particles are not "weakly_compressible_particles"'
        end select

        ! intialize acceleration and density rate-of-change arrays
        do i = 1, ps_fluid%size
            ps_fluid%dvxdt(:, i) = 0._fp
            ps_fluid%dvxdt(ps_fluid%ndims, i) = self%g
            ps_fluid%drhodt(i) = 0._fp
        end do

        ! perform sweep
        do k = 1, self%pairs%npairs_total
            i = self%pairs%pair_ij(1, k)
            j = self%pairs%pair_ij(2, k)
            call artificial_viscosity_monaghan1994( &
                ps_fluid%ndims, ps_fluid%x(:, i), ps_fluid%x(:, j), ps_fluid%v(:, i), ps_fluid%v(:, j), ps_fluid%rho(i), &
                ps_fluid%rho(j), self%h, self%h, ps_fluid%c(i), ps_fluid%c(j), ps_fluid%mass(i), &
                ps_fluid%mass(j), ps_fluid%dvxdt(:, i), ps_fluid%dvxdt(:, j), self%pairs%dwdx(:, k), &
                self%artvisc_alpha, self%artvisc_beta)
            call isotropic_pressure_force(ps_fluid%ndims, ps_fluid%p(i), ps_fluid%p(j), ps_fluid%rho(i), ps_fluid%rho(j), &
                                          ps_fluid%mass(i), ps_fluid%mass(j), ps_fluid%dvxdt(:, i), &
                                          ps_fluid%dvxdt(:, j), self%pairs%dwdx(:, k))
            call continuity_density(ps_fluid%ndims, ps_fluid%v(:, i), ps_fluid%v(:, j), ps_fluid%mass(i), ps_fluid%mass(j), &
                                    ps_fluid%drhodt(i), ps_fluid%drhodt(j), self%pairs%dwdx(:, k))
        end do

    end subroutine fluid_self_sweep

    subroutine fluid_fluid_sweep(self)
        class(fluid_fluid_interaction), intent(inout):: self
        integer:: i, j, k
        real(fp):: dvxdtj(self%ps_rhs%ndims), drhodtj
        class(linear_eos_particles), pointer:: ps_fluid1, ps_fluid2

        select type (ps => self%ps_lhs)
        class is (linear_eos_particles)
            ps_fluid1 => ps
        class default
            error stop 'LHS Particles are not "weakly_compressible_particles"'
        end select

        select type (ps => self%ps_rhs)
        class is (linear_eos_particles)
            ps_fluid2 => ps
        class default
            error stop 'RHS Particles are not "weakly_compressible_particles"'
        end select

        do k = 1, self%pairs%npairs_total
            i = self%pairs%pair_ij(1, k)
            j = self%pairs%pair_ij(2, k)
            call artificial_viscosity_monaghan1994( &
                ps_fluid1%ndims, ps_fluid1%x(:, i), ps_fluid2%x(:, j), ps_fluid1%v(:, i), ps_fluid2%v(:, j), &
                ps_fluid1%rho(i), ps_fluid2%rho(j), self%h, self%h, ps_fluid1%c(i), &
                ps_fluid2%c(j), ps_fluid1%mass(i), ps_fluid2%mass(j), ps_fluid1%dvxdt(:, i), &
                dvxdtj, self%pairs%dwdx(:, k), self%artvisc_alpha, self%artvisc_beta &
                )
            call isotropic_pressure_force(ps_fluid1%ndims, ps_fluid1%p(i), ps_fluid2%p(j), ps_fluid1%rho(i), &
                                          ps_fluid2%rho(j), ps_fluid1%mass(i), ps_fluid2%mass(j), ps_fluid1%dvxdt(:, i), &
                                          dvxdtj, self%pairs%dwdx(:, k) &
                                          )
            call continuity_density( &
                ps_fluid1%ndims, ps_fluid1%v(:, i), ps_fluid2%v(:, j), ps_fluid1%mass(i), &
                ps_fluid2%mass(j), ps_fluid1%drhodt(i), drhodtj, self%pairs%dwdx(:, k) &
                )
        end do

    end subroutine fluid_fluid_sweep

end module weakly_compressible_interactions
