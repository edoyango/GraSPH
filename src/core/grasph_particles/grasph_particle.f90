!> @file grasph_particle.f90
!> @brief Module containing the base particle class.
!> @author Edward Yang
!> @date 2025-10-02
module grasph_particle_m

    use grasph_constants_m, only: fp, ndims

    implicit none

    private

    !> @brief base particle type.
    type:: base_particles_t
        !> @brief the ID of the particle.
        integer, dimension(:), allocatable:: id
        !> @brief An integer indicating the "type" of the particle. Not currently used for anything and may be removed.
        integer, dimension(:), allocatable:: type
        !> @brief The particle's position.
        real(fp), dimension(:, :), allocatable:: x
        !> @brief Lagrangian velocity of the particle.
        real(fp), dimension(:, :), allocatable:: v
        !> @brief The density of the particle.
        real(fp), dimension(:), allocatable:: rho
        !> @brief The mass of the particle.
        real(fp), dimension(:), allocatable:: mass
        !> @brief The local speed of sound associated with the particle.
        real(fp), dimension(:), allocatable:: c
        !> @brief The acceleration of the particle.
        real(fp), dimension(:, :), allocatable:: dvxdt
        !> @brief The density rate-of-change of the particle.
        real(fp), dimension(:), allocatable:: drhodt
        !> @brief The allocated size of the particle data.
        integer:: size = 0
        !> @brief Whether the particles have been allocated.
        logical, private:: allocated_ = .false.
    contains
        procedure:: init => base_particles_init
        procedure:: deallocate => base_particles_deallocate
        procedure:: allocated => particles_allocated
    end type base_particles_t

    public:: base_particles_t, base_particles_init, base_particles_deallocate

contains

    subroutine base_particles_init(self, n)
        class(base_particles_t), intent(inout):: self
        integer, intent(in):: n

        ! make sure everything is deallocated.
        call self%deallocate()

        ! allocate integer arrays and initialize with 0.
        allocate ( &
            self%id(n), &
            self%type(n), &
            source=0 &
            )

        ! allocate real arrays and initialize with 0..
        allocate ( &
            self%x(ndims, n), &
            self%v(ndims, n), &
            self%rho(n), &
            self%mass(n), &
            self%c(n), &
            self%dvxdt(ndims, n), &
            self%drhodt(n), &
            source=0._fp &
            )

        self%size = n
        self%allocated_ = .false.

    end subroutine base_particles_init

    subroutine base_particles_deallocate(self)
        class(base_particles_t), intent(inout):: self

        if (allocated(self%id)) deallocate (self%id)
        if (allocated(self%type)) deallocate (self%type)
        if (allocated(self%x)) deallocate (self%x)
        if (allocated(self%v)) deallocate (self%v)
        if (allocated(self%rho)) deallocate (self%rho)
        if (allocated(self%mass)) deallocate (self%mass)
        if (allocated(self%c)) deallocate (self%c)
        if (allocated(self%dvxdt)) deallocate (self%dvxdt)
        if (allocated(self%drhodt)) deallocate (self%drhodt)

        self%allocated_ = .false.

    end subroutine base_particles_deallocate

    !> @brief Returns wither the particle data has been allocated.
    !> @param self The particles to check allocation of.
    pure logical function particles_allocated(self)
        class(base_particles_t), intent(in):: self
        particles_allocated = self%allocated_
    end function particles_allocated

end module grasph_particle_m
