!> @file grasph_pair_sets.f90
!> @brief Module containing types that are used to hold information about particle interactions
!> @author Edward Yang
!> @date 2025-06-09
module grasph_pair_sets

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles
    use grasph_pairs, only: particle_pairs, cell_list_search
    use grasph_kernels, only: grasph_base_kernel

    implicit none

    private

    !> @brief base interaction class to extend into new ones. The base isn't that useful on its own.
    type:: particle_interactions_base
        !> @brief Pointer to LHS particles involved in the interaction.
        !>        Note that these can't be overriden, so if LHS extended particles' properties are needed, then a new pointer is
        !>        needed in the derived type.
        class(base_particles), pointer:: ps_lhs
        !> @brief Pointer to RHS particles involved in the interaction.
        !>        Note that these can't be overriden, so if LHS extended particles' properties are needed, then a new pointer is
        !>        needed in the derived type.
        class(base_particles), pointer:: ps_rhs
        !> @brief The pairs of particles found either in ps_lhs or between ps_lhs and ps_rhs (depends on whether
        !>        ps_rhs was passed to initializer).
        type(particle_pairs):: pairs
        !> @brief Whether particle_interactions has been initialized.
        logical:: initialized = .false.
        !> @brief Whether the particle_interactions describes ps_lhs interaction with itself, or with ps_rhs.
        logical:: is_pair_set = .false.
    contains
        !> @brief A placeholder subroutine intended to update particles' state that depend on interpolated information.
        !>        E.g. Updating virtual particles' data, which requires a sweep.
        !>        Does nothing in the base type and is intended to be overridden when necessary.
        procedure:: sweep_prologue => donothing_sweep
        !> @brief A placeholder subroutine intended to calculate time-evolving data's rate-of-change.
        !>        E.g. acceleration. Does nothing in base type and is intended to be overridden when necessary.
        procedure:: sweep => donothing_sweep
        !> @brief A placeholder subroutine intended to perform any particle shifting via position or velocity adjustments.
        procedure:: shift => donothing_shift
        !> @brief If is_pair_set is .true., calculates pairs within ps_lhs, else pairs between ps_lhs and ps_rhs.
        !>        Shouldn't need to be overridden.
        procedure:: find_pairs => particle_interactions_base_find_pairs
        !> @brief Initializes base data (pointers and pairs). Intended not to be overridden and instead be used in derived types'
        !>        initializer.
        procedure:: base_init => particle_interactions_base_init
        !> @brief Initializes data (pointers and pairs). Intended to be overridden by derived types' intializer.
        procedure:: init => particle_interactions_base_init
    end type particle_interactions_base

    !> @brief A simple container class to facilitate polymorphism if extended particle_interactions.
    type:: particle_interactions_container
        !> @brief The polymorphic container to be allocated to particle_interactions_base or its derivatives.
        class(particle_interactions_base), allocatable:: pi
    end type particle_interactions_container

    public:: particle_interactions_base, particle_interactions_container

contains

    !> @brief A do-nothing placeholder subroutine. Intended to be overriden when particles have
    !>        meaningful sweeps to perform.
    !> @param self The particle interactions class which performing the sweep.
    subroutine donothing_sweep(self)
        class(particle_interactions_base), intent(inout):: self
    end subroutine donothing_sweep

    !> @brief A do-nothing placeholder subroutine. Intended to be overriden when a position/velocity adjustment is implemented.
    !>        AKA for particle shifting.
    subroutine donothing_shift(self, dt)
        class(particle_interactions_base), intent(inout):: self
        real(fp), intent(in):: dt
    end

    !> @brief The subroutine to find pairs of particles contained in particle_interactions.
    !>        Adapts to whether the particle_interactions instance is a pair set or not.
    !> @param self The particle interactions class to find pairs within.
    !> @param cutoff The interacting distance of particles.
    !> @param kernel The SPH kernel used to calculate values and kernel gradient values from.
    subroutine particle_interactions_base_find_pairs(self, cutoff, kernel)
        class(particle_interactions_base), intent(inout):: self
        real(fp), intent(in):: cutoff
        class(grasph_base_kernel), intent(in):: kernel

        if (self%is_pair_set) then
            call cell_list_search(self%ps_lhs%x, self%ps_rhs%x, self%ps_rhs%size, cutoff, kernel, self%pairs)
        else
            call cell_list_search(self%ps_lhs%x, cutoff, kernel, self%pairs)
        end if
    end subroutine particle_interactions_base_find_pairs

    !> @brief The base initializer of particle_interactions instances. Not intended to be overridden
    !>        and instead intended to be called within an extended type's initializer subroutine.
    !> @param npairs_per_particle Maximum number of particle interactions expected for each LHS particle.
    !> @param self The particle_interactions instance to initialize.
    !> @param ps_lhs The LHS particles to be attached to the instance.
    !> @param ps_rhs The RHS particles to be attached to the instance.
    subroutine particle_interactions_base_init(self, npairs_per_particle, ps_lhs, ps_rhs)
        class(particle_interactions_base), intent(out):: self
        integer, intent(in):: npairs_per_particle
        class(base_particles), target, intent(in):: ps_lhs
        class(base_particles), target, optional, intent(in):: ps_rhs
        self%ps_lhs => ps_lhs
        if (present(ps_rhs)) then
            self%ps_rhs => ps_rhs
            self%is_pair_set = .true.
        end if
        call self%pairs%init(ps_lhs%size, npairs_per_particle, ps_lhs%ndims)
    end subroutine particle_interactions_base_init

end module grasph_pair_sets
