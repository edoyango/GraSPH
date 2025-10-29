!> @file weakly_compressible_particles.f90
!> @brief Module containing weakly compressible particles type and methods
!> @author Edward Yang
!> @date 2025-09-21
module weakly_compressible_particles_m

    use grasph_constants_m, only: fp, ndims, pi
    use grasph_particle_system_m, only: base_particle_t, particle_system_t, base_state_updater_t

    implicit none

    private

    type, extends(base_particle_t):: eos_particle_t
        real(fp):: p
    end type eos_particle_t

    type, extends(eos_particle_t):: eos_ghost_particle_t
        class(eos_particle_t), pointer:: original
    end type eos_ghost_particle_t

    type, extends(base_state_updater_t):: linear_eos_state_updater_t
        real(fp):: rho_ref
    contains
        procedure:: update_state => linear_eos_update_state
    end type linear_eos_state_updater_t

    type, extends(linear_eos_state_updater_t):: tait_eos_state_updater_t
        integer:: gamma = 7
    contains
        procedure:: update_state => tait_eos_update_state
    end type tait_eos_state_updater_t

    type, extends(base_state_updater_t):: ghost_state_updater_t
        real(fp):: surface_normal(ndims)
    contains
        procedure:: update_state => ghost_state_update
    end type ghost_state_updater_t

    integer, parameter:: ntensor_elems = ndims*ndims
    integer, parameter:: ntensor_offaxis_elems = (ntensor_elems - ndims)/2
    integer, parameter:: ntensor_elems_voigt = ndims + ntensor_offaxis_elems

    type, extends(eos_particle_t):: eos_viscous_stress_particle_t
        real(fp):: strain_rate(ntensor_elems_voigt), stress(ntensor_elems_voigt)
    end type eos_viscous_stress_particle_t

    type, extends(eos_viscous_stress_particle_t):: eos_viscous_stress_ghost_particle_t
        class(eos_viscous_stress_particle_t), pointer:: original
    end type eos_viscous_stress_ghost_particle_t

    type, extends(linear_eos_state_updater_t):: dp_visco_elastic_state_updater_t
        real(fp):: friction_angle = pi/6._fp ! 30 degrees
        real(fp):: cohesion = 0._fp
    contains
        procedure:: update_state => dp_visco_elastic_state_update
    end type dp_visco_elastic_state_updater_t

    type, extends(base_state_updater_t):: eos_viscous_stress_ghost_state_updater_t
        real(fp):: surface_normal(ndims)
    contains
        procedure:: update_state => eos_viscous_stress_ghost_state_update
    end type eos_viscous_stress_ghost_state_updater_t

    public:: eos_particle_t, eos_ghost_particle_t, linear_eos_state_updater_t, tait_eos_state_updater_t, ghost_state_updater_t, &
             eos_viscous_stress_particle_t, dp_visco_elastic_state_updater_t, eos_viscous_stress_ghost_particle_t, &
             eos_viscous_stress_ghost_state_updater_t, ntensor_elems_voigt

contains

    !> @brief The linear state equation to update stress using the particles' speed of sound (c),
    !>        density (rho), and reference density (rho_ref). Overrides particle system's state_update
    !>        subroutine.
    !> @param self The particles' pressure to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine linear_eos_update_state(self, ps, n, dt)
        class(linear_eos_state_updater_t), intent(in):: self
        integer, intent(in):: n
        class(base_particle_t), intent(inout):: ps(n)
        real(fp), intent(in), optional:: dt
        integer:: i
        select type (ps_eos => ps)
        class is (eos_particle_t)
            do i = 1, n
                ps_eos(i)%p = ps_eos(i)%c**2*(ps_eos(i)%rho - self%rho_ref)
            end do
        class default
            error stop "eos_particle_t required"
        end select
    end subroutine linear_eos_update_state

    !> @brief The Tait state equation to update stress using the particles' speed of sound (c),
    !>        density (rho), and reference density (rho_ref). Overrides particle system's state_update
    !>        subroutine.
    !> @param self The particles' pressure to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine tait_eos_update_state(self, ps, n, dt)
        class(tait_eos_state_updater_t), intent(in):: self
        integer, intent(in):: n
        class(base_particle_t), intent(inout):: ps(n)
        real(fp), intent(in), optional:: dt
        integer:: i
        select type (ps_eos => ps)
        class is (eos_particle_t)
            do i = 1, n
                ps_eos(i)%p = self%rho_ref*ps_eos(i)%c*ps_eos(i)%c/real(self%gamma, kind=fp)* &
                              ((ps_eos(i)%rho/self%rho_ref)**self%gamma - 1._fp)
            end do
        class default
            error stop "eos_particle_t required"
        end select
    end subroutine tait_eos_update_state

    subroutine ghost_state_update(self, ps, n, dt)
        class(ghost_state_updater_t), intent(in):: self
        integer, intent(in):: n
        class(base_particle_t), intent(inout):: ps(n)
        real(fp), optional, intent(in):: dt
        integer:: i
        real(fp):: projection(ndims)

        select type (ps_ghost => ps)
        class is (eos_ghost_particle_t)
            do i = 1, n
                projection(:) = dot_product(ps_ghost(i)%original%v(:), self%surface_normal(:))*self%surface_normal(:)
                ps_ghost(i)%v(:) = ps_ghost(i)%original%v(:) - 2._fp*projection(:)
                ps_ghost(i)%rho = ps_ghost(i)%original%rho
                ps_ghost(i)%mass = ps_ghost(i)%original%mass
                ps_ghost(i)%p = ps_ghost(i)%original%p
                ps_ghost(i)%c = ps_ghost(i)%original%c
            end do
        class default
            error stop "Expected self%particles to be eos_ghost_particle_t."
        end select

    end subroutine ghost_state_update

    subroutine dp_visco_elastic_state_update(self, ps, n, dt)
        class(dp_visco_elastic_state_updater_t), intent(in):: self
        integer, intent(in):: n
        class(base_particle_t), intent(inout):: ps(n)
        real(fp), intent(in), optional:: dt
        integer:: i, d
        real(fp):: mag_strain_rate
        class(eos_viscous_stress_particle_t), pointer:: ps_ve(:)

        select type (ps => ps)
        class is (eos_viscous_stress_particle_t)
            ps_ve => ps
        class default
            error stop "ps is required to be eos_viscous_stress_particle_t"
        end select

        call linear_eos_update_state(self, ps_ve, n, dt)

        do i = 1, n
            mag_strain_rate = 0._fp
            do d = 1, ndims
                mag_strain_rate = mag_strain_rate + ps_ve(i)%strain_rate(d)**2
            end do
            do d = 1, ntensor_offaxis_elems
                mag_strain_rate = mag_strain_rate + 2._fp*ps_ve(i)%strain_rate(ndims + d)**2
            end do

            mag_strain_rate = max(sqrt(mag_strain_rate), tiny(1._fp))
            ps_ve(i)%stress(:) = (self%cohesion + tan(self%friction_angle)*ps_ve(i)%p)/mag_strain_rate*ps_ve(i)%strain_rate(:)
            ps_ve(i)%stress(1:ndims) = ps_ve(i)%stress(1:ndims) - ps_ve(i)%p
        end do

    end subroutine dp_visco_elastic_state_update

    subroutine eos_viscous_stress_ghost_state_update(self, ps, n, dt)
        class(eos_viscous_stress_ghost_state_updater_t), intent(in):: self
        integer, intent(in):: n
        class(base_particle_t), intent(inout):: ps(n)
        real(fp), optional, intent(in):: dt
        integer:: i
        real(fp):: projection(ndims)

        select type (ps_ghost => ps)
        class is (eos_viscous_stress_ghost_particle_t)
            do i = 1, n
                projection(:) = dot_product(ps_ghost(i)%original%v(:), self%surface_normal(:))*self%surface_normal(:)
                ps_ghost(i)%v(:) = ps_ghost(i)%original%v(:) - 2._fp*projection(:)
                ps_ghost(i)%rho = ps_ghost(i)%original%rho
                ps_ghost(i)%mass = ps_ghost(i)%original%mass
                ps_ghost(i)%p = ps_ghost(i)%original%p
                ps_ghost(i)%c = ps_ghost(i)%original%c
                ps_ghost(i)%stress(:) = ps_ghost(i)%original%stress(:)
            end do
        class default
            error stop "Expected self%particles to be eos_viscous_stress_ghost_particle_t."
        end select

    end subroutine eos_viscous_stress_ghost_state_update

end module weakly_compressible_particles_m
