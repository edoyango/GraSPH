!> @file grasph_particle.f90
!> @brief Module containing the base particle class.
!> @author Edward Yang
!> @date 2025-10-02
module grasph_particle_m

    use grasph_constants_m, only: fp, ndims

    implicit none

    private

    !> @brief base particle type.
    type:: base_particle_t
        !> @brief the ID of the particle.
        integer:: id = 0
        !> @brief An integer indicating the "type" of the particle. Not currently used for anything and may be removed.
        integer:: type = 0
        !> @brief The particle's position.
        real(fp):: x(ndims) = 0._fp
        !> @brief Lagrangian velocity of the particle.
        real(fp):: v(ndims) = 0._fp
        !> @brief The density of the particle.
        real(fp):: rho = 0._fp
        !> @brief The mass of the particle.
        real(fp):: mass = 0._fp
        !> @brief The local speed of sound associated with the particle.
        real(fp):: c = 0._fp
        !> @brief The acceleration of the particle.
        real(fp):: dvxdt(ndims) = 0._fp
        !> @brief The density rate-of-change of the particle.
        real(fp):: drhodt = 0._fp
    end type base_particle_t

    public:: base_particle_t

end module grasph_particle_m
