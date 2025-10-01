!> @file weakly_compressible_particles.f90
!> @brief Module containing weakly compressible particles type and methods
!> @author Edward Yang
!> @date 2025-09-21
module weakly_compressible_particles

    use grasph_constants, only: fp
    use grasph_particles, only: base_particle, base_particles, base_state_updater

    type, extends(base_particle):: linear_eos_particle
        real(fp):: p
    end type linear_eos_particle

    !> @brief particle type which adds pressure, determined from density with a linear EOS
    type, extends(base_particles):: linear_eos_particles
    contains
        !> @brief Custom intializer to initialize pressure and reference density.
        !>        Also calls base_init to initialize base data.
        procedure:: init => wcp_init
    end type linear_eos_particles

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

contains

    !> @brief Custom init function for weakly-compressible particles. Will also initialize base
    !>        particles' data.
    !> @param self The weakly-compressible particles to initialize.
    !> @param n Number of particles to allocate space for.
    !> @param d Spatial dimensions of the particles.
    !> @param name A label to give the particles. Used to label output/terminal information.
    !> @param rho_ref Reference density used in the linear EOS.
    subroutine wcp_init(self, n, name, ps_template, state_updater)
        class(linear_eos_particles), intent(inout):: self
        integer, intent(in):: n
        character(*), intent(in):: name
        class(linear_eos_particle), optional, intent(in):: ps_template
        class(base_state_updater), optional, intent(in):: state_updater
        type(linear_eos_particle):: ps_default
        if (present(ps_template)) then
            call self%base_init(n, name, ps_template, state_updater=state_updater)
            call self%register_io%register_variable(ps_template, "p", ps_template%p)
        else
            call self%base_init(n, name, ps_default, state_updater=state_updater)
            call self%register_io%register_variable(ps_default, "p", ps_default%p)
        end if
    end subroutine wcp_init

    !> @brief The linear state equation to update stress using the particles' speed of sound (c),
    !>        density (rho), and reference density (rho_ref). Overrides base_particles' state_update
    !>        subroutine.
    !> @param self The particles' pressure to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine linear_eos_update_state(self, ps, n, dt)
        class(linear_eos_state_updater), intent(in):: self
        integer, intent(in):: n
        class(base_particle), intent(inout):: ps(n)
        real(fp), intent(in), optional:: dt
        integer:: i
        select type (ps_eos => ps)
        class is (linear_eos_particle)
            do i = 1, n
                ps_eos(i)%p = ps_eos(i)%c**2*(ps_eos(i)%rho - self%rho_ref)
            end do
        class default
            error stop "linear_eos_particle required"
        end select
    end subroutine linear_eos_update_state

    !> @brief The Tait state equation to update stress using the particles' speed of sound (c),
    !>        density (rho), and reference density (rho_ref). Overrides base_particles' state_update
    !>        subroutine.
    !> @param self The particles' pressure to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine tait_eos_update_state(self, ps, n, dt)
        class(tait_eos_state_updater), intent(in):: self
        integer, intent(in):: n
        class(base_particle), intent(inout):: ps(n)
        real(fp), intent(in), optional:: dt
        integer:: i
        select type (ps_eos => ps)
        class is (linear_eos_particle)
            do i = 1, n
                ps_eos(i)%p = self%rho_ref*ps_eos(i)%c*ps_eos(i)%c/real(self%gamma, kind=fp)* &
                              ((ps_eos(i)%rho/self%rho_ref)**self%gamma - 1._fp)
            end do
        class default
            error stop "linear_eos_particle required"
        end select
    end subroutine tait_eos_update_state

end module weakly_compressible_particles
