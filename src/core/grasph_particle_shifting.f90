!> @file grasph_particle_shifting.f90
!> @brief Module containing subroutines for calculating pair-wise contributions for particle-shifting.
!> @author Edward Yang
!> @date 2025-09-22
module grasph_particle_shifting_m

    use grasph_constants_m, only: fp, ndims, max_name_len
    use grasph_particle_m, only: base_particles_t
    use grasph_system_interactions_m, only: base_sweeper_t
    use grasph_pairs_m, only: particle_pairs_t
    use grasph_particle_system_m, only: particle_system_t

    implicit none

    private

    !> @brief XSPH shifter class.
    type, extends(base_sweeper_t):: xsph_shifter_t
        !> @brief Coefficient controlling strength of shifting.
        real(fp):: epsilon = 0.5_fp
    contains
        !> @brief Performs XSPH particle shifting for one particle system, as described in Monaghan 1994.
        procedure:: sweep_1system => xsph_shift_1system
        !> @brief Performs XSPH particle shifting for two particle systems, as described in Monaghan 1994.
        procedure:: sweep_2system => xsph_shift_2system
        !> @brief Performs XSPH particle shifting for two particle systems, as described in Monaghan 1994.
        procedure:: sweep_2system_norhsupdate => xsph_shift_2system_norhsupdate
        !> @brief Readonly name of the shifter.
        procedure, nopass:: name => xsph_shifter_name
    end type xsph_shifter_t

    public:: xsph_shifter_t

contains

    !> @brief Updates i, j particles' positions with XSPH shifting.
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
    subroutine xsph_shift_ij(xi, xj, vi, vj, rhoi, rhoj, massi, massj, w, dt, epsilon)
        real(fp), intent(inout):: xi(ndims), xj(ndims)
        real(fp), intent(in):: vi(ndims), vj(ndims), rhoi, rhoj, massi, massj, w, dt, epsilon
        real(fp):: dv(ndims), mrho

        mrho = 0.5_fp*(rhoi + rhoj)
        dv(:) = epsilon*(vj(:) - vi(:))/mrho*w
        xi(:) = xi(:) + massj*dv(:)*dt
        xj(:) = xj(:) - massi*dv(:)*dt

    end subroutine xsph_shift_ij

    !> @brief Performs XSPH particle shifting on particles within a single particle system, as described in Monaghan 1994.
    !> @param self The shifter class. Used to access epsilon.
    !> @param pairs The class storing particle pair index information.
    !> @param ps the particles involved in the interactions.
    !> @param dt The time-step increment.
    subroutine xsph_shift_1system(self, pairs, ps, dt)
        class(xsph_shifter_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), target, intent(inout):: ps
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        real(fp):: dummyx(ndims)

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call xsph_shift_ij( &
                ps%x(:, i), ps%x(:, j), ps%v(:, i), ps%v(:, j), ps%rho(i), ps%rho(j), ps%mass(i), ps%mass(j), pairs%w(k), dt, &
                self%epsilon &
                )
        end do

    end subroutine xsph_shift_1system

    !> @brief Performs XSPH particle shifting on particles within two particle systems, as described in Monaghan 1994.
    !> @param self The shifter class. Used to access epsilon.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS particles involved in the interactions.
    !> @param ps_rhs The RHS particles involved in the interactions.
    !> @param dt The time-step increment.
    subroutine xsph_shift_2system(self, pairs, ps_lhs, ps_rhs, dt)
        class(xsph_shifter_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), target, intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call xsph_shift_ij( &
                ps_lhs%x(:, i), ps_rhs%x(:, j), ps_lhs%v(:, i), ps_rhs%v(:, j), ps_lhs%rho(i), ps_rhs%rho(j), ps_lhs%mass(i), &
                ps_rhs%mass(j), pairs%w(k), dt, self%epsilon &
                )
        end do

    end subroutine xsph_shift_2system

    !> @brief Performs XSPH particle shifting on particles within two particle systems, as described in Monaghan 1994.
    !> @param self The shifter class. Used to access epsilon.
    !> @param pairs The class storing particle pair index information.
    !> @param ps_lhs the LHS particles involved in the interactions.
    !> @param ps_rhs The RHS particles involved in the interactions.
    !> @param dt The time-step increment.
    subroutine xsph_shift_2system_norhsupdate(self, pairs, ps_lhs, ps_rhs, dt)
        class(xsph_shifter_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), target, intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        real(fp):: dummyx(ndims)

        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            call xsph_shift_ij( &
                ps_lhs%x(:, i), dummyx, ps_lhs%v(:, i), ps_rhs%v(:, j), ps_lhs%rho(i), ps_rhs%rho(j), ps_lhs%mass(i), &
                ps_rhs%mass(j), pairs%w(k), dt, self%epsilon &
                )
        end do

    end subroutine xsph_shift_2system_norhsupdate

    !> @brief Read-only name of the shifter.
    !> @returns The name of the shifter.
    pure character(max_name_len) function xsph_shifter_name()
        xsph_shifter_name = "xsph_shifter_t"
    end function xsph_shifter_name

end module grasph_particle_shifting_m
