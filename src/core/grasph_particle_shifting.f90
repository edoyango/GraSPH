!> @file grasph_particle_shifting.f90
!> @brief Module containing subroutines for calculating pair-wise contributions for particle-shifting.
!> @author Edward Yang
!> @date 2025-09-22
module grasph_particle_shifting

    use grasph_constants, only: fp
    use grasph_pair_sets, only: base_shifter
    use grasph_pairs, only: particle_pairs
    use grasph_particles, only: base_particles

    implicit none

    private

    !> @brief XSPH shifter class.
    type, extends(base_shifter):: xsph_shifter
        !> @brief Coefficient controlling strength of shifting.
        real(fp):: epsilon = 0.5_fp
    contains
        !> @brief Performs XSPH particle shifting, as described in Monaghan 1994.
        procedure:: shift => xsph_shift
    end type xsph_shifter

    public:: xsph_shifter

contains

    !> @brief Updates i, j particles' positions with XSPH shifting.
    !> @param ndims Number of spatial dimensions.
    !> @param xi Position of particle i [m].
    !> @param xj Position of particle j [m].
    !> @param vi Velocity of particle i [m].
    !> @param vj Velocity of particle j [m].
    !> @param rhoi Density of particle i [kg/m^3].
    !> @param rhoj Density of particle j [kg/m^3].
    !> @param massi Mass of particle i [kg].
    !> @param massj Mass of particle j [kg].
    !> @param w Kernel weight (W_ij) [dimensionless].
    !> @param dt Timestep increment [s].
    !> @param epsilon XSPH shifting coefficient [dimensionless].
    subroutine xsph_shift_ij(ndims, xi, xj, vi, vj, rhoi, rhoj, massi, massj, w, dt, epsilon)
        integer, intent(in):: ndims
        real(fp), intent(inout):: xi(ndims), xj(ndims)
        real(fp), intent(in):: vi(ndims), vj(ndims), rhoi, rhoj, massi, massj, w, dt, epsilon
        real(fp):: dv(ndims), mrho

        mrho = 0.5_fp*(rhoi + rhoj)
        dv(:) = epsilon*(vj(:) - vi(:))/mrho*w
        xi(:) = xi(:) + massj*dv(:)*dt
        xj(:) = xj(:) - massi*dv(:)*dt

    end subroutine xsph_shift_ij

    !> @brief Performs XSPH particle shifting, as described in Monaghan 1994.
    !> @param self The shifter class. Used to access epsilon.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS particles involved in the interactions.
    !> @param ps_rhs The RHS particles involved in the interactions. ps_rhs will not be passed in if not associated in the owning
    !>        particle_interactions class.
    !> @param dt The time-step increment.
    subroutine xsph_shift(self, pairs, ps_lhs, ps_rhs, dt)
        class(xsph_shifter), intent(in):: self
        type(particle_pairs), intent(in):: pairs
        class(base_particles), intent(inout):: ps_lhs
        class(base_particles), optional, intent(inout):: ps_rhs
        real(fp), intent(in):: dt
        integer:: i, j, k
        real(fp):: dummyx(ps_lhs%ndims)

        if (present(ps_rhs)) then
            do k = 1, pairs%npairs_total
                i = pairs%pair_ij(1, k)
                j = pairs%pair_ij(2, k)
                if (self%update_rhs) then
                    call xsph_shift_ij(ps_lhs%ndims, ps_lhs%x(:, i), ps_rhs%x(:, j), ps_lhs%v(:, i), ps_rhs%v(:, j), &
                                       ps_lhs%rho(i), ps_rhs%rho(j), ps_lhs%mass(i), ps_rhs%mass(j), pairs%w(k), dt, &
                                       self%epsilon)
                else
                    call xsph_shift_ij(ps_lhs%ndims, ps_lhs%x(:, i), dummyx, ps_lhs%v(:, i), ps_rhs%v(:, j), ps_lhs%rho(i), &
                                       ps_rhs%rho(j), ps_lhs%mass(i), ps_rhs%mass(j), pairs%w(k), dt, self%epsilon)
                end if
            end do
        else
            do k = 1, pairs%npairs_total
                i = pairs%pair_ij(1, k)
                j = pairs%pair_ij(2, k)
                call xsph_shift_ij(ps_lhs%ndims, ps_lhs%x(:, i), ps_lhs%x(:, j), ps_lhs%v(:, i), ps_lhs%v(:, j), ps_lhs%rho(i), &
                                   ps_lhs%rho(j), ps_lhs%mass(i), ps_lhs%mass(j), pairs%w(k), dt, self%epsilon)
            end do
        end if

    end subroutine xsph_shift

end module grasph_particle_shifting
