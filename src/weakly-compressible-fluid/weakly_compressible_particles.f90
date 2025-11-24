!> @file weakly_compressible_particles.f90
!> @brief Module containing weakly compressible particles type and methods
!> @author Edward Yang
!> @date 2025-09-21
module weakly_compressible_particles_m

    use grasph_constants_m, only: fp, ndims, pi
    use grasph_particle_m, only: base_particles_t, base_particles_init, base_particles_deallocate
    use grasph_particle_system_m, only: particle_system_t, base_state_updater_t

    implicit none

    private

    !> @brief Weakly compressible particle type.
    type, extends(base_particles_t):: eos_particles_t
        !> @brief Pressure
        real(fp), allocatable:: p(:)
    contains
        procedure:: init => eos_particles_init
        procedure:: deallocate => eos_particles_deallocate
    end type eos_particles_t

    !> @brief Weakly compressible ghost particle type.
    type, extends(eos_particles_t):: eos_ghost_particles_t
        !> @brief Pointer to the original particles the ghost particles are based on.
        class(eos_particles_t), pointer:: ps_original => null()
        !> @brief The index in ps_original that the given ghost particle is based on.
        integer, allocatable:: idx_original(:)
    contains
        procedure:: init => eos_ghost_particles_init
        procedure:: deallocate => eos_ghost_particles_deallocate
    end type eos_ghost_particles_t

    !> @brief State updater for eos particles using linear state equation.
    type, extends(base_state_updater_t):: linear_eos_state_updater_t
        !> @brief Reference density.
        real(fp):: rho_ref = 0._fp
    contains
        !> @brief Linear equation of state update subroutine.
        procedure:: update_state => linear_eos_update_state
    end type linear_eos_state_updater_t

    !> @brief State updater for eos particles using Tait equation.
    type, extends(linear_eos_state_updater_t):: tait_eos_state_updater_t
        !> @brief Gamma constant.
        integer:: gamma = 7
    contains
        !> @brief Tait equation of state update routine.
        procedure:: update_state => tait_eos_update_state
    end type tait_eos_state_updater_t

    !> @brief Weakly compressible ghost particle state updater.
    type, extends(base_state_updater_t):: ghost_state_updater_t
        !> @brief The unit normal vector used to calculate the ghost particles' velocity (for enforcing free-slip conditions).
        real(fp):: surface_normal(ndims)
    contains
        !> @brief The ghost particle state updater.
        procedure:: update_state => ghost_state_update
    end type ghost_state_updater_t

    !> @brief Number of elements in the cauchy stress matrix.
    integer, parameter:: ntensor_elems = ndims*ndims
    !> @brief Number of off-axis elements in the cauchy stress matrix.
    integer, parameter:: ntensor_offaxis_elems = (ntensor_elems - ndims)/2
    !> @brief Number of elements in the cauchy stress matrix in Voigt notation.
    integer, parameter:: ntensor_elems_voigt = ndims + ntensor_offaxis_elems

    !> @brief Weakly compressible particle type with stress and strain rate tensors in voigt notation.
    type, extends(eos_particles_t):: eos_viscous_stress_particles_t
        !> @brief Strain rate tensor.
        real(fp), allocatable:: strain_rate(:, :)
        !> @brief Cauchy stress tensor.
        real(fp), allocatable:: stress(:, :)
    contains
        procedure:: init => eos_viscous_stress_particles_init
        procedure:: deallocate => eos_viscous_stress_particles_deallocate
    end type eos_viscous_stress_particles_t

    !> @brief Weakly compressible ghost particle type with stress and strain rate tensors in voigt notation.
    type, extends(eos_viscous_stress_particles_t):: eos_viscous_stress_ghost_particles_t
        !> @brief Pointer to the original particles the ghost particles are based on.
        class(eos_viscous_stress_particles_t), pointer:: ps_original
        !> @brief The index in ps_original that the given ghost particle is based on.
        integer, allocatable:: idx_original(:)
    contains
        procedure:: init => eos_viscous_stress_ghost_particles_init
        procedure:: deallocate => eos_viscous_stress_ghost_particles_deallocate
    end type eos_viscous_stress_ghost_particles_t

    !> @brief Stress and pressure state updater using visco-plasticity with Drucker-Prager-like yield criterion, and linear equation
    !>        of state.
    type, extends(linear_eos_state_updater_t):: dp_visco_elastic_state_updater_t
        !> @brief Friction angle for DP-like yield criterion.
        real(fp):: friction_angle = pi/6._fp ! 30 degrees
        !> @brief Cohesion for DP-like yield criterion.
        real(fp):: cohesion = 0._fp
    contains
        !> @brief Updates particles' stress using visco-plasticity with Drucker-Prager-like yield criterion, and linear equation
        !>        of state.
        procedure:: update_state => dp_visco_elastic_state_update
    end type dp_visco_elastic_state_updater_t

    !> @brief Weakly compressible ghost particle with stress tensor state updater.
    type, extends(base_state_updater_t):: eos_viscous_stress_ghost_state_updater_t
        !> @brief The unit normal vector used to calculate the ghost particles' velocity (for enforcing free-slip conditions).
        real(fp):: surface_normal(ndims)
    contains
        !> @brief The ghost particle state updater.
        procedure:: update_state => eos_viscous_stress_ghost_state_update
    end type eos_viscous_stress_ghost_state_updater_t

    public:: eos_particles_t, eos_ghost_particles_t, linear_eos_state_updater_t, tait_eos_state_updater_t, ghost_state_updater_t, &
             eos_viscous_stress_particles_t, dp_visco_elastic_state_updater_t, eos_viscous_stress_ghost_particles_t, &
             eos_viscous_stress_ghost_state_updater_t, ntensor_elems_voigt

contains

    !> @brief Deallocates all the arrays in the eos_particles_t class.
    !> @param self The class to deallocate members of.
    subroutine eos_particles_deallocate(self)
        class(eos_particles_t), intent(inout):: self
        call base_particles_deallocate(self)
        if (allocated(self%p)) deallocate (self%p)
    end subroutine eos_particles_deallocate

    !> @brief Allocates particle data arrays and initialises everything to zero.
    !> @param self The particles to initialise.
    !> @param n The number of particles to allocate.
    subroutine eos_particles_init(self, n)
        class(eos_particles_t), intent(inout):: self
        integer, intent(in):: n
        call self%deallocate()
        call base_particles_init(self, n)
        allocate (self%p(n), source=0._fp)
    end subroutine eos_particles_init

    !> @brief Deallocates all the arrays in the eos_ghost_particles_t class.
    !> @param self The class to deallocate members of.
    subroutine eos_ghost_particles_deallocate(self)
        class(eos_ghost_particles_t), intent(inout):: self
        call eos_particles_deallocate(self)
        self%ps_original => null()
        if (allocated(self%idx_original)) deallocate (self%idx_original)
    end subroutine eos_ghost_particles_deallocate

    !> @brief Allocates particle data arrays and initialises everything to zero.
    !> @param self The particles to initialise.
    !> @param n The number of particles to allocate.
    subroutine eos_ghost_particles_init(self, n)
        class(eos_ghost_particles_t), intent(inout):: self
        integer, intent(in):: n
        call self%deallocate()
        call eos_particles_init(self, n)
        allocate (self%idx_original(n), source=0)
    end subroutine eos_ghost_particles_init

    !> @brief Deallocates all the arrays in the eos_viscous_stress_particles_t class.
    !> @param self The class to deallocate members of.
    subroutine eos_viscous_stress_particles_deallocate(self)
        class(eos_viscous_stress_particles_t), intent(inout):: self
        call eos_particles_deallocate(self)
        if (allocated(self%strain_rate)) deallocate (self%strain_rate)
        if (allocated(self%stress)) deallocate (self%stress)
    end subroutine eos_viscous_stress_particles_deallocate

    !> @brief Allocates particle data arrays and initialises everything to zero.
    !> @param self The particles to initialise.
    !> @param n The number of particles to allocate.
    subroutine eos_viscous_stress_particles_init(self, n)
        class(eos_viscous_stress_particles_t), intent(inout):: self
        integer, intent(in):: n
        call eos_particles_init(self, n)
        allocate (self%strain_rate(ntensor_elems_voigt, n), source=0._fp)
        allocate (self%stress(ntensor_elems_voigt, n), source=0._fp)
    end subroutine eos_viscous_stress_particles_init

    !> @brief Deallocates all the arrays in the eos_viscous_stress_ghost_particles_t class.
    !> @param self The class to deallocate members of.
    subroutine eos_viscous_stress_ghost_particles_deallocate(self)
        class(eos_viscous_stress_ghost_particles_t), intent(inout):: self
        call eos_viscous_stress_particles_deallocate(self)
        self%ps_original => null()
        if (allocated(self%idx_original)) deallocate (self%idx_original)
    end subroutine eos_viscous_stress_ghost_particles_deallocate

    !> @brief Allocates particle data arrays and initialises everything to zero.
    !> @param self The particles to initialise.
    !> @param n The number of particles to allocate.
    subroutine eos_viscous_stress_ghost_particles_init(self, n)
        class(eos_viscous_stress_ghost_particles_t), intent(inout):: self
        integer, intent(in):: n
        call eos_viscous_stress_particles_init(self, n)
        allocate (self%idx_original(n), source=0)
    end subroutine eos_viscous_stress_ghost_particles_init

    !> @brief The linear state equation to update stress using the particles' speed of sound (c),
    !>        density (rho), and reference density (rho_ref). Overrides particle system's state_update
    !>        subroutine.
    !> @param self The state updater holding reference density constant.
    !> @param ps The particles who's pressure are to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine linear_eos_update_state(self, ps, dt)
        class(linear_eos_state_updater_t), intent(in):: self
        class(base_particles_t), target, intent(inout):: ps
        real(fp), intent(in), optional:: dt
        integer:: i
        select type (ps_eos => ps)
        class is (eos_particles_t)
            do i = 1, ps%size
                ps_eos%p(i) = ps_eos%c(i)**2*(ps_eos%rho(i) - self%rho_ref)
            end do
        class default
            error stop "eos_particles_t required"
        end select
    end subroutine linear_eos_update_state

    !> @brief The Tait state equation to update stress using the particles' speed of sound (c),
    !>        density (rho), and reference density (rho_ref). Overrides particle system's state_update
    !>        subroutine.
    !> @param self The state updater holding reference density and gamma constants.
    !> @param ps The particles who's pressure is to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine tait_eos_update_state(self, ps, dt)
        class(tait_eos_state_updater_t), intent(in):: self
        class(base_particles_t), target, intent(inout):: ps
        real(fp), intent(in), optional:: dt
        integer:: i
        select type (ps_eos => ps)
        class is (eos_particles_t)
            do i = 1, ps%size
                ps_eos%p(i) = self%rho_ref*ps_eos%c(i)**2/real(self%gamma, kind=fp)* &
                              ((ps_eos%rho(i)/self%rho_ref)**self%gamma - 1._fp)
            end do
        class default
            error stop "eos_particles_t required"
        end select
    end subroutine tait_eos_update_state

    !> @brief Updates ghost particles' state using its original particles' properties and the boundary surface unit normal vector.
    !> @param self The state updater holding boundary surface normal.
    !> @param ps The ghost particles who's state is to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine ghost_state_update(self, ps, dt)
        class(ghost_state_updater_t), intent(in):: self
        class(base_particles_t), target, intent(inout):: ps
        real(fp), optional, intent(in):: dt
        integer:: i, i_original
        real(fp):: projection(ndims), v_original(ndims)

        select type (ps_ghost => ps)
        class is (eos_ghost_particles_t)
            do i = 1, ps%size
                i_original = ps_ghost%idx_original(i)
                v_original(:) = ps_ghost%ps_original%v(:, i_original)
                projection(:) = dot_product(v_original, self%surface_normal(:))*self%surface_normal(:)
                ps_ghost%v(:, i) = v_original(:) - 2._fp*projection(:)
                ps_ghost%rho(i) = ps_ghost%ps_original%rho(i_original)
                ps_ghost%mass(i) = ps_ghost%ps_original%mass(i_original)
                ps_ghost%c(i) = ps_ghost%ps_original%c(i_original)
                ps_ghost%p(i) = ps_ghost%ps_original%p(i_original)
            end do
        class default
            error stop "Expected self%particles to be eos_ghost_particles_t."
        end select

    end subroutine ghost_state_update

    !> @brief Updates stress of weakly-compressible particles with stress using a visco-plastic stress-strain relation with DP-like
    !> @brief yield criterion and linear equation of state.
    !> @param self The state updater holding reference density, friction angle, and cohesion.
    !> @param ps The particles who's stress is to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine dp_visco_elastic_state_update(self, ps, dt)
        class(dp_visco_elastic_state_updater_t), intent(in):: self
        class(base_particles_t), target, intent(inout):: ps
        real(fp), intent(in), optional:: dt
        integer:: i, d
        real(fp):: mag_strain_rate
        class(eos_viscous_stress_particles_t), pointer:: ps_ve

        select type (ps => ps)
        class is (eos_viscous_stress_particles_t)
            ps_ve => ps
        class default
            error stop "ps is required to be eos_viscous_stress_particles_t"
        end select

        ! first calculate pressure component of stress tensor
        call linear_eos_update_state(self, ps_ve, dt)

        do i = 1, ps%size
            ! calculate second invariant of deformation rate tensor.
            mag_strain_rate = 0._fp
            do d = 1, ndims
                mag_strain_rate = mag_strain_rate + ps_ve%strain_rate(d, i)**2
            end do
            do d = 1, ntensor_offaxis_elems
                mag_strain_rate = mag_strain_rate + 2._fp*ps_ve%strain_rate(ndims + d, i)**2
            end do
            mag_strain_rate = max(sqrt(mag_strain_rate), tiny(1._fp)) ! tiny(1) to make sure non-zero

            ! viscous stress with yield criterion
            ps_ve%stress(:, i) = (self%cohesion + tan(self%friction_angle)*ps_ve%p(i))/mag_strain_rate*ps_ve%strain_rate(:, i)
            ! minus pressure along principal components.
            ps_ve%stress(1:ndims, i) = ps_ve%stress(1:ndims, i) - ps_ve%p(i)
        end do

    end subroutine dp_visco_elastic_state_update

    !> @brief Updates ghost particles' state using its original particles' properties and the boundary surface unit normal vector.
    !> @param self The state updater holding boundary surface normal.
    !> @param ps The ghost particles who's state is to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine eos_viscous_stress_ghost_state_update(self, ps, dt)
        class(eos_viscous_stress_ghost_state_updater_t), intent(in):: self
        class(base_particles_t), target, intent(inout):: ps
        real(fp), optional, intent(in):: dt
        integer:: i, i_original
        real(fp):: projection(ndims), v_original(ndims)

        select type (ps_ghost => ps)
        class is (eos_viscous_stress_ghost_particles_t)
            do i = 1, ps%size
                i_original = ps_ghost%idx_original(i)
                v_original = ps_ghost%ps_original%v(:, i_original)
                projection(:) = dot_product(v_original(:), self%surface_normal(:))*self%surface_normal(:)
                ps_ghost%v(:, i) = v_original(:) - 2._fp*projection(:)
                ps_ghost%rho(i) = ps_ghost%ps_original%rho(i_original)
                ps_ghost%mass(i) = ps_ghost%ps_original%mass(i_original)
                ps_ghost%p(i) = ps_ghost%ps_original%p(i_original)
                ps_ghost%c(i) = ps_ghost%ps_original%c(i_original)
                ps_ghost%stress(:, i) = ps_ghost%ps_original%stress(:, i_original)
            end do
        class default
            error stop "Expected self%particles to be eos_viscous_stress_ghost_particles_t."
        end select

    end subroutine eos_viscous_stress_ghost_state_update

end module weakly_compressible_particles_m
