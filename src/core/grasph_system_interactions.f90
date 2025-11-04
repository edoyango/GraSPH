!> @file grasph_system_interactions.f90
!> @brief Module containing types that are used to hold information about particle interactions
!> @author Edward Yang
!> @date 2025-06-09
module grasph_system_interactions_m

    use grasph_constants_m, only: fp
    use grasph_particle_system_m, only: particle_system_t
    use grasph_pairs_m, only: particle_pairs_t, cell_list_search
    use grasph_kernels_m, only: base_kernel_t

    implicit none

    private

    !> @brief Manages the interaction between a particle_system_t and itself or between two particle_system_t instances.
    type:: system_interaction_t
        !> @brief Pointer to LHS particle system involved in the interaction.
        class(particle_system_t), pointer:: psys_lhs
        !> @brief Pointer to RHS particles involved in the interaction.
        class(particle_system_t), pointer:: psys_rhs
        !> @brief The pairs of particles found either in psys_lhs or between psys_lhs and psys_rhs (depends on whether
        !>        psys_rhs was passed to initializer).
        type(particle_pairs_t):: pairs
        !> @brief Whether system_interaction_t has been initialized.
        logical:: initialized = .false.
        !> @brief Whether the system_interaction_t describes psys_lhs interaction with itself, or with psys_rhs.
        logical:: is_pair_set = .false.
        !> @brief Overridable "strategy" class that performs timestep setup for the interaction.
        class(base_sweeper_t), allocatable:: timestep_setuper
        !> @brief Overridable "strategy" class that performs sweep prologue.
        class(base_sweeper_t), allocatable:: prologue_sweeper
        !> @brief Overridable "strategy" class that performs sweep.
        class(base_sweeper_t), allocatable:: sweeper
        !> @brief Overridable "strategy" class that performs shift.
        class(base_sweeper_t), allocatable:: shifter
    contains
        !> @brief Called at start of every time-step for any special setup. E.g. creating ghost particles.
        procedure:: do_timestep_setup
        !> @brief Updates particles' state that depend on interpolated information. E.g. Updating virtual particles' data, which
        !>        requires a sweep. Does so by using the prologue_sweeper's sweep method.
        procedure:: do_sweep_prologue
        !> @brief Calculates time-evolving data's rate-of-change e.g. acceleration. Does so by using the sweeper's sweep method.
        procedure:: do_sweep
        !> @brief Performs any particle shifting via position or velocity adjustments. Does so by using the shifter's shift method.
        procedure:: do_shift
        !> @brief If is_pair_set is .true., calculates pairs within psys_lhs, else pairs between psys_lhs and psys_rhs.
        procedure:: find_pairs => particle_interactions_find_pairs
        !> @brief Initializes data (pointers, pairs, and strategy classes).
        procedure:: init => particle_interactions_init
    end type system_interaction_t

    !> @brief Base "strategy" class whose sweep method is used to update time-evolving data's rate-of-change e.g. acceleration.
    !>        extensions of this class override the sweep to, for example, work with different particle types and implement
    !>        different physics.
    type, abstract:: base_sweeper_t
        !> @brief Controls whether to update the RHS particles (if they're associated).
        logical:: update_rhs = .true.
        !> @brief Controls whether the sweep initializes particles' rate-of-change data.
        logical:: initialize = .true.
    contains
        procedure(sweep_1system_interface), deferred:: sweep_1system
        procedure(sweep_2system_interface), deferred:: sweep_2system
        procedure(sweep_2system_interface), deferred:: sweep_2system_norhsupdate
    end type base_sweeper_t

    abstract interface
        !> @brief Interface that describes the sweep method of the base_sweeper strategy class.
        !> @param self The sweeper class. Used to access constants.
        !> @param pairs The class storing particle pair index information.
        !> @param psys the particles involved in the interactions.
        !> @param dt The timestep size.
        subroutine sweep_1system_interface(self, pairs, psys, dt)
            import:: base_sweeper_t, particle_pairs_t, particle_system_t, fp
            class(base_sweeper_t), intent(in):: self
            type(particle_pairs_t), intent(in):: pairs
            class(particle_system_t), intent(inout):: psys
            real(fp), optional, intent(in):: dt
        end subroutine sweep_1system_interface
        !> @brief Interface that describes the sweep method of the base_sweeper strategy class.
        !> @param self The sweeper class. Used to access constants.
        !> @param pairs The class storing particle pair index information.
        !> @param psys_lhs the LHS particles involved in the interactions.
        !> @param psys_rhs The RHS particles involved in the interactions. psys_rhs will not be passed in if not associated in the
        !>        owning system_interaction_t class.
        !> @param dt The timestep size.
        subroutine sweep_2system_interface(self, pairs, psys_lhs, psys_rhs, dt)
            import:: base_sweeper_t, particle_pairs_t, particle_system_t, fp
            class(base_sweeper_t), intent(in):: self
            type(particle_pairs_t), intent(in):: pairs
            class(particle_system_t), intent(inout):: psys_lhs, psys_rhs
            real(fp), optional, intent(in):: dt
        end subroutine sweep_2system_interface
    end interface

    !> @brief A default sweeper which does nothing when the sweep method is called.
    type, extends(base_sweeper_t):: default_sweeper_t
    contains
        !> @brief A sweep that does nothing with the input particle system.
        procedure:: sweep_1system => donothing_sweep_1system
        !> @brief A sweep that does nothing with the input particle systems.
        procedure:: sweep_2system => donothing_sweep_2system
        !> @brief A sweep that does nothing with the input particle systems.
        procedure:: sweep_2system_norhsupdate => donothing_sweep_2system
    end type default_sweeper_t

    public:: system_interaction_t, base_sweeper_t

contains

    !> @brief Executes setup at start of timestep that requires information from 2 particle systems e.g. generating ghost particle
    !>        positions. Utilizes the timestep_setuper strategy member class.
    !> @param self The system interaction with the two particle systems that need to be updated.
    subroutine do_timestep_setup(self)
        class(system_interaction_t), intent(inout):: self

        ! check that timestep setuper has been allocated
        if (.not. allocated(self%timestep_setuper)) error stop "timestep setuper not allocated in system_interaction_t."

        ! pass in psys_rhs if associated
        if (associated(self%psys_rhs)) then
            if (self%timestep_setuper%update_rhs) then
                call self%timestep_setuper%sweep_2system(self%pairs, self%psys_lhs, self%psys_rhs)
            else
                call self%timestep_setuper%sweep_2system_norhsupdate(self%pairs, self%psys_lhs, self%psys_rhs)
            end if
        else
            call self%timestep_setuper%sweep_1system(self%pairs, self%psys_lhs)
        end if

    end subroutine do_timestep_setup

    !> @brief Executes sweep step prior to particle state setup. THis is for updating particle properties like virtual particles'
    !>        velocity or density, or calculating strain rate. Uses the prologue_sweeper strategy member class.
    !> @param self The system interaction to perform the sweep between.
    subroutine do_sweep_prologue(self)
        class(system_interaction_t), intent(inout):: self

        ! check that prologue_sweeper has been allocated
        if (.not. allocated(self%prologue_sweeper)) error stop "prologue sweeper not allocated in system_interaction_t."

        ! pass in psys_rhs if associated
        if (associated(self%psys_rhs)) then
            if (self%prologue_sweeper%update_rhs) then
                call self%prologue_sweeper%sweep_2system(self%pairs, self%psys_lhs, self%psys_rhs)
            else
                call self%prologue_sweeper%sweep_2system_norhsupdate(self%pairs, self%psys_lhs, self%psys_rhs)
            end if
        else
            call self%prologue_sweeper%sweep_1system(self%pairs, self%psys_lhs)
        end if

    end subroutine do_sweep_prologue

    !> @brief Executes sweep step for calculating rate of changes e.g. motion or density. Uses the sweeper strategy member class.
    !> @param self The system interaction to perform the sweep between.
    !> @param dt The time-step size.
    subroutine do_sweep(self, dt)
        class(system_interaction_t), intent(inout):: self
        real(fp), optional, intent(in):: dt

        ! check that sweeper has been allocated
        if (.not. allocated(self%sweeper)) error stop "sweeper not allocated in system_interaction_t."

        ! pass in psys_rhs if associated
        if (associated(self%psys_rhs)) then
            if (self%sweeper%update_rhs) then
                call self%sweeper%sweep_2system(self%pairs, self%psys_lhs, self%psys_rhs, dt=dt)
            else
                call self%sweeper%sweep_2system_norhsupdate(self%pairs, self%psys_lhs, self%psys_rhs, dt=dt)
            end if
        else
            call self%sweeper%sweep_1system(self%pairs, self%psys_lhs, dt=dt)
        end if

    end subroutine do_sweep

    !> @brief Executes a sweep in the particle system to calculate properties required for shifting particles' positions.
    !> @param self The system interaction to perfrom the sweep between.
    !> @param dt The timestep size.
    subroutine do_shift(self, dt)
        class(system_interaction_t), intent(inout):: self
        real(fp), intent(in):: dt

        ! check that sweeper has been allocated
        if (.not. allocated(self%shifter)) error stop "shifter not allocated in system_interaction_t."

        ! pass in psys_rhs if associated
        if (associated(self%psys_rhs)) then
            if (self%shifter%update_rhs) then
                call self%shifter%sweep_2system(self%pairs, self%psys_lhs, self%psys_rhs, dt)
            else
                call self%shifter%sweep_2system_norhsupdate(self%pairs, self%psys_lhs, self%psys_rhs, dt)
            end if
        else
            call self%shifter%sweep_1system(self%pairs, self%psys_lhs, dt=dt)
        end if

    end subroutine do_shift

    !> @brief A do-nothing placeholder subroutine.
    !> @param self The sweeper class. Used to access constants.
    !> @param pairs The class storing particle pair index information.
    !> @param psys the particles involved in the interactions.
    !> @param dt The time-step size.
    subroutine donothing_sweep_1system(self, pairs, psys, dt)
        class(default_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys
        real(fp), optional, intent(in):: dt
    end subroutine donothing_sweep_1system

    !> @brief A do-nothing placeholder subroutine. Intended to be overriden when particles have
    !>        meaningful sweeps to perform.
    !> @param self The sweeper class. Used to access constants.
    !> @param pairs The class storing particle pair index information.
    !> @param psys_lhs the LHS particles involved in the interactions.
    !> @param psys_rhs The RHS particles involved in the interactions. psys_rhs will not be passed in if not associated in the
    !>        owning system_interaction_t class.
    !> @param dt The time-step size.
    subroutine donothing_sweep_2system(self, pairs, psys_lhs, psys_rhs, dt)
        class(default_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs, psys_rhs
        real(fp), optional, intent(in):: dt
    end subroutine donothing_sweep_2system

    !> @brief The subroutine to find pairs of particles contained in system_interaction_t.
    !>        Adapts to whether the system_interaction_t instance is a pair set or not.
    !> @param self The particle interactions class to find pairs within.
    !> @param cutoff The interacting distance of particles.
    !> @param kernel The SPH kernel used to calculate values and kernel gradient values from.
    subroutine particle_interactions_find_pairs(self, cutoff, kernel)
        class(system_interaction_t), intent(inout):: self
        real(fp), intent(in):: cutoff
        class(base_kernel_t), intent(in):: kernel

        if (self%is_pair_set) then
            call cell_list_search( &
                self%psys_lhs%particles(1:self%psys_lhs%size), &
                self%psys_rhs%particles(1:self%psys_rhs%size), &
                cutoff, &
                kernel, &
                self%pairs &
                )
        else
            call cell_list_search( &
                self%psys_lhs%particles(1:self%psys_lhs%size), &
                cutoff, &
                kernel, &
                self%pairs &
                )
        end if
    end subroutine particle_interactions_find_pairs

    !> @brief The base initializer of system_interaction_t instances.
    !> @param npairs_per_particle Maximum number of particle interactions expected for each LHS particle.
    !> @param self The system_interaction_t instance to initialize.
    !> @param psys_lhs The LHS particles to be attached to the instance.
    !> @param psys_rhs The RHS particles to be attached to the instance.
    !> @param timestep_setuper The strategy class that performs any setup needed at the start of the timestep - before the pair
    !>        finding has occurred.
    !> @param prologue_sweeper The sweeper to use in the prologue sweep in this interaction.
    !> @param sweeper The sweeper to use in this interaction.
    !> @param shifter The shifter to use in this interaction.
    subroutine particle_interactions_init(self, npairs_per_particle, psys_lhs, psys_rhs, timestep_setuper, prologue_sweeper, &
                                          sweeper, shifter)
        class(system_interaction_t), intent(out):: self
        integer, intent(in):: npairs_per_particle
        class(particle_system_t), target, intent(in):: psys_lhs
        class(particle_system_t), target, optional, intent(in):: psys_rhs
        class(base_sweeper_t), optional, intent(in):: timestep_setuper, prologue_sweeper, sweeper, shifter
        type(default_sweeper_t):: tmp_base_sweeper
        self%psys_lhs => psys_lhs
        if (present(psys_rhs)) then
            self%psys_rhs => psys_rhs
            self%is_pair_set = .true.
        end if

        call self%pairs%init(psys_lhs%size, npairs_per_particle)

        if (present(timestep_setuper)) then
            allocate (self%timestep_setuper, source=timestep_setuper)
        else
            allocate (self%timestep_setuper, source=tmp_base_sweeper)
        end if

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
            allocate (self%shifter, source=tmp_base_sweeper)
        end if
    end subroutine particle_interactions_init

end module grasph_system_interactions_m
