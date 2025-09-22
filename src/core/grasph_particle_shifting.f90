!> @file grasph_particle_shifting.f90
!> @brief Module containing subroutines for calculating pair-wise contributions for particle-shifting.
!> @author Edward Yang
!> @date 2025-09-22
module grasph_particle_shifting

    use grasph_constants, only: fp
    use grasph_pair_sets, only: particle_interactions_base

    implicit none

    private
    public:: xsph_shift

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
    subroutine xsph_shift(ndims, xi, xj, vi, vj, rhoi, rhoj, massi, massj, w, dt, epsilon)
        integer, intent(in):: ndims
        real(fp), intent(inout):: xi(ndims), xj(ndims)
        real(fp), intent(in):: vi(ndims), vj(ndims), rhoi, rhoj, massi, massj, w, dt, epsilon
        real(fp):: dv(ndims), mrho

        mrho = 0.5_fp*(rhoi + rhoj)
        dv(:) = epsilon*(vj(:) - vi(:))/mrho*w
        xi(:) = xi(:) + massj*dv(:)*dt
        xj(:) = xj(:) - massi*dv(:)*dt

    end subroutine xsph_shift

end module grasph_particle_shifting
