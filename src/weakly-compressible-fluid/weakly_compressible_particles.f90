!> @file weakly_compressible_particles.f90
!> @brief Module containing weakly compressible particles type and methods
!> @author Edward Yang
!> @date 2025-09-21
module weakly_compressible_particles

    use grasph_constants, only: fp
    use grasph_particles, only: base_particle_t, particle_system_t, base_state_updater

    implicit none

    private

    type, extends(base_particle_t):: eos_particle
        real(fp):: p
    end type eos_particle

    type, extends(base_state_updater):: linear_eos_state_updater
        real(fp):: rho_ref
    contains
        procedure:: update_state => linear_eos_update_state
    end type linear_eos_state_updater

    type, extends(linear_eos_state_updater):: tait_eos_state_updater
        integer:: gamma = 7
    contains
        procedure:: update_state => tait_eos_update_state
    end type tait_eos_state_updater

    public:: eos_particle, linear_eos_state_updater, tait_eos_state_updater

contains

    !> @brief The linear state equation to update stress using the particles' speed of sound (c),
    !>        density (rho), and reference density (rho_ref). Overrides particle system's state_update
    !>        subroutine.
    !> @param self The particles' pressure to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine linear_eos_update_state(self, ps, n, dt)
        class(linear_eos_state_updater), intent(in):: self
        integer, intent(in):: n
        class(base_particle_t), intent(inout):: ps(n)
        real(fp), intent(in), optional:: dt
        integer:: i
        select type (ps_eos => ps)
        class is (eos_particle)
            do i = 1, n
                ps_eos(i)%p = ps_eos(i)%c**2*(ps_eos(i)%rho - self%rho_ref)
            end do
        class default
            error stop "eos_particle required"
        end select
    end subroutine linear_eos_update_state

    !> @brief The Tait state equation to update stress using the particles' speed of sound (c),
    !>        density (rho), and reference density (rho_ref). Overrides particle system's state_update
    !>        subroutine.
    !> @param self The particles' pressure to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine tait_eos_update_state(self, ps, n, dt)
        class(tait_eos_state_updater), intent(in):: self
        integer, intent(in):: n
        class(base_particle_t), intent(inout):: ps(n)
        real(fp), intent(in), optional:: dt
        integer:: i
        select type (ps_eos => ps)
        class is (eos_particle)
            do i = 1, n
                ps_eos(i)%p = self%rho_ref*ps_eos(i)%c*ps_eos(i)%c/real(self%gamma, kind=fp)* &
                              ((ps_eos(i)%rho/self%rho_ref)**self%gamma - 1._fp)
            end do
        class default
            error stop "eos_particle required"
        end select
    end subroutine tait_eos_update_state

end module weakly_compressible_particles
