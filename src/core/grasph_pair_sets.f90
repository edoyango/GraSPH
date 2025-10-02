!> @file grasph_pair_sets.f90
!> @brief Module containing types that are used to hold information about particle interactions
!> @author Edward Yang
!> @date 2025-06-09
module grasph_pair_sets

    use grasph_constants, only: fp
    use grasph_particles, only: particle_system_t
    use grasph_pairs, only: particle_pairs, cell_list_search
    use grasph_kernels, only: grasph_base_kernel

    implicit none

    private

    !> @brief base interaction class to extend into new ones. The base isn't that useful on its own.
    type:: particle_interactions
        !> @brief Pointer to LHS particles involved in the interaction.
        !>        Note that these can't be overriden, so if LHS extended particles' properties are needed, then a new pointer is
        !>        needed in the derived type.
        class(particle_system_t), pointer:: psys_lhs
        !> @brief Pointer to RHS particles involved in the interaction.
        !>        Note that these can't be overriden, so if LHS extended particles' properties are needed, then a new pointer is
        !>        needed in the derived type.
        class(particle_system_t), pointer:: psys_rhs
        !> @brief The pairs of particles found either in psys_lhs or between psys_lhs and psys_rhs (depends on whether
        !>        psys_rhs was passed to initializer).
        type(particle_pairs):: pairs
        !> @brief Whether particle_interactions has been initialized.
        logical:: initialized = .false.
        !> @brief Whether the particle_interactions describes psys_lhs interaction with itself, or with psys_rhs.
        logical:: is_pair_set = .false.
        !> @brief Overridable "strategy" class that performs sweep prologue.
        class(base_sweeper), allocatable:: prologue_sweeper
        !> @brief Overridable "strategy" class that performs sweep.
        class(base_sweeper), allocatable:: sweeper
        !> @brief Overridable "strategy" class that performs shift.
        class(base_shifter), allocatable:: shifter
    contains
        !> @brief Subroutine to update particles' state that depend on interpolated information. E.g. Updating virtual particles'
        !>        data, which requires a sweep. Does so by using the prologue_sweeper's sweep method.
        procedure:: do_sweep_prologue
        !> @brief Subroutine to calculate time-evolving data's rate-of-change e.g. acceleration. Does so by using the sweeper's
        !>        sweep method.
        procedure:: do_sweep
        !> @brief Subroutine to perform any particle shifting via position or velocity adjustments. Does so by using the shifter's
        !>        shift method.
        procedure:: do_shift
        !> @brief If is_pair_set is .true., calculates pairs within psys_lhs, else pairs between psys_lhs and psys_rhs.
        procedure:: find_pairs => particle_interactions_find_pairs
        !> @brief Initializes data (pointers, pairs, and strategy classes).
        procedure:: init => particle_interactions_init
    end type particle_interactions

    !> @brief Base "strategy" class whose sweep method is used to update time-evolving data's rate-of-change e.g. acceleration.
    !>        extensions of this class override the sweep to, for example, work with different particle types and implement
    !>        different physics.
    type:: base_sweeper
        !> @brief Controls whether to update the RHS particles (if they're associated).
        logical:: update_rhs = .true.
        !> @brief Controls whether the sweep initializes particles' rate-of-change data.
        logical:: initialize = .true.
    contains
        !> @brief Update particles' rate-of-change data by sweeping through particle pairs.
        procedure:: sweep => donothing_sweep
    end type base_sweeper

    !> @brief Base "strategy" class whose shift method is used to perform any particle shifting via position or velocity &
    !>        adjustments.
    type:: base_shifter
        !> @brief Controls whether to update the RHS particles (if they're associated).
        logical:: update_rhs = .true.
    contains
        !> @brief Perform particle shifting.
        procedure:: shift => donothing_shift
    end type base_shifter

    public:: particle_interactions, base_sweeper, base_shifter

contains

    subroutine do_sweep_prologue(self)
        class(particle_interactions), intent(inout):: self

        ! check that prologue_sweeper has been allocated
        if (.not. allocated(self%prologue_sweeper)) error stop "prologue sweeper not allocated in particle_interactions."

        ! pass in psys_rhs if associated
        if (associated(self%psys_rhs)) then
            call self%prologue_sweeper%sweep(self%pairs, self%psys_lhs, self%psys_rhs)
        else
            call self%prologue_sweeper%sweep(self%pairs, self%psys_lhs)
        end if

    end subroutine do_sweep_prologue

    subroutine do_sweep(self)
        class(particle_interactions), intent(inout):: self

        ! check that sweeper has been allocated
        if (.not. allocated(self%sweeper)) error stop "sweeper not allocated in particle_interactions."

        ! pass in psys_rhs if associated
        if (associated(self%psys_rhs)) then
            call self%sweeper%sweep(self%pairs, self%psys_lhs, self%psys_rhs)
        else
            call self%sweeper%sweep(self%pairs, self%psys_lhs)
        end if

    end subroutine do_sweep

    subroutine do_shift(self, dt)
        class(particle_interactions), intent(inout):: self
        real(fp), intent(in):: dt

        ! check that sweeper has been allocated
        if (.not. allocated(self%shifter)) error stop "shifter not allocated in particle_interactions."

        ! pass in psys_rhs if associated
        if (associated(self%psys_rhs)) then
            call self%shifter%shift(self%pairs, self%psys_lhs, self%psys_rhs, dt)
        else
            call self%shifter%shift(self%pairs, self%psys_lhs, dt=dt)
        end if

    end subroutine do_shift

    !> @brief A do-nothing placeholder subroutine. Intended to be overriden when particles have
    !>        meaningful sweeps to perform.
    !> @param self The sweeper class. Used to access constants.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS particles involved in the interactions.
    !> @param psys_rhs The RHS particles involved in the interactions. psys_rhs will not be passed in if not associated in the
    !>        owning particle_interactions class.
    subroutine donothing_sweep(self, pairs, psys_lhs, psys_rhs)
        class(base_sweeper), intent(in):: self
        type(particle_pairs), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
    end subroutine donothing_sweep

    !> @brief A do-nothing placeholder subroutine. Intended to be overriden when particles are to be
    !>        shifted at the end of a timestep.
    !> @param self The shifter class. Used to access constants.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS particles involved in the interactions.
    !> @param psys_rhs The RHS particles involved in the interactions. psys_rhs will not be passed in if not associated in the
    !>        owning particle_interactions class.
    !> @param dt The time-step increment.
    subroutine donothing_shift(self, pairs, psys_lhs, psys_rhs, dt)
        class(base_shifter), intent(in):: self
        type(particle_pairs), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        real(fp), intent(in):: dt
    end subroutine donothing_shift

    !> @param self The particle interactions class which performing the sweep.
    subroutine donothing_sweep_old(self)
        class(particle_interactions), intent(inout):: self
    end subroutine donothing_sweep_old

    !> @brief The subroutine to find pairs of particles contained in particle_interactions.
    !>        Adapts to whether the particle_interactions instance is a pair set or not.
    !> @param self The particle interactions class to find pairs within.
    !> @param cutoff The interacting distance of particles.
    !> @param kernel The SPH kernel used to calculate values and kernel gradient values from.
    subroutine particle_interactions_find_pairs(self, cutoff, kernel)
        class(particle_interactions), intent(inout):: self
        real(fp), intent(in):: cutoff
        class(grasph_base_kernel), intent(in):: kernel

        if (self%is_pair_set) then
            call cell_list_search(self%psys_lhs%particles, self%psys_rhs%particles, cutoff, kernel, self%pairs)
        else
            call cell_list_search(self%psys_lhs%particles, cutoff, kernel, self%pairs)
        end if
    end subroutine particle_interactions_find_pairs

    !> @brief The base initializer of particle_interactions instances. Not intended to be overridden
    !>        and instead intended to be called within an extended type's initializer subroutine.
    !> @param npairs_per_particle Maximum number of particle interactions expected for each LHS particle.
    !> @param self The particle_interactions instance to initialize.
    !> @param psys_lhs The LHS particles to be attached to the instance.
    !> @param psys_rhs The RHS particles to be attached to the instance.
    !> @param sweeper The sweeper to use in this interaction.
    !> @param shifter The shifter to use in this interaction.
    subroutine particle_interactions_init(self, npairs_per_particle, psys_lhs, psys_rhs, prologue_sweeper, sweeper, shifter)
        class(particle_interactions), intent(out):: self
        integer, intent(in):: npairs_per_particle
        class(particle_system_t), target, intent(in):: psys_lhs
        class(particle_system_t), target, optional, intent(in):: psys_rhs
        class(base_sweeper), optional, intent(in):: prologue_sweeper, sweeper
        class(base_shifter), optional, intent(in):: shifter
        type(base_sweeper):: tmp_base_sweeper
        type(base_shifter):: tmp_base_shifter
        self%psys_lhs => psys_lhs
        if (present(psys_rhs)) then
            self%psys_rhs => psys_rhs
            self%is_pair_set = .true.
        end if
        call self%pairs%init(psys_lhs%size, npairs_per_particle)

        if (present(prologue_sweeper)) then
            allocate (self%prologue_sweeper, source=prologue_sweeper)
        else
            allocate (self%prologue_sweeper, source=tmp_base_sweeper)
        end if

        if (present(sweeper)) then
            allocate (self%sweeper, source=sweeper)
        else
            allocate (self%sweeper, source=tmp_base_sweeper)
        end if

        if (present(shifter)) then
            allocate (self%shifter, source=shifter)
        else
            allocate (self%shifter, source=tmp_base_shifter)
        end if
    end subroutine particle_interactions_init

end module grasph_pair_sets
