!> @file grasph_particle.f90
!> @brief Module containing the base particle class.
!> @author Edward Yang
!> @date 2025-10-02
module grasph_particle_m

    use grasph_constants, only: fp, ndims

    implicit none

    private

    !> @brief base particle type.
    type:: base_particle_t
        !> @brief the ID of the particle.
        integer:: id
        !> @brief An integer indicating the "type" of the particle. Not currently used for anything and may be removed.
        integer:: type
        !> @brief The particle's position.
        real(fp):: x(ndims)
        !> @brief Lagrangian velocity of the particle.
        real(fp):: v(ndims)
        !> @brief The density of the particle.
        real(fp):: rho
        !> @brief The mass of the particle.
        real(fp):: mass
        !> @brief The local speed of sound associated with the particle.
        real(fp):: c
        !> @brief The acceleration of the particle.
        real(fp):: dvxdt(ndims)
        !> @brief The density rate-of-change of the particle.
        real(fp):: drhodt
    end type base_particle_t

    public:: base_particle_t

end module grasph_particle_m
