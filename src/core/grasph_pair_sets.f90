module grasph_pair_sets

    use grasph_constants, only: fp, ndims
    use grasph_particles, only: base_particles
    use grasph_pairs, only: particle_pairs, cell_list_search
    use grasph_kernels, only: grasph_base_kernel

    implicit none

    private

    type:: particle_interactions_base
        class(base_particles), pointer:: ps_lhs, ps_rhs
        type(particle_pairs):: pairs
        logical:: initialized = .false., is_pair_set = .false.
    contains
        procedure:: sweep_prologue => donothing_sweep, sweep => donothing_sweep
        procedure:: find_pairs => particle_interactions_base_find_pairs
        ! base_init is made a separate method so overrides of init can still use it
        procedure:: base_init => particle_interactions_base_init, init => particle_interactions_base_init
    end type particle_interactions_base

    type:: particle_interactions_container
        class(particle_interactions_base), allocatable:: pi
    end type particle_interactions_container

    public:: particle_interactions_base, particle_interactions_container

contains

    subroutine donothing_sweep(self)
        class(particle_interactions_base), intent(inout):: self
    end subroutine donothing_sweep

    subroutine particle_interactions_base_find_pairs(self, cutoff, kernel)
        class(particle_interactions_base), intent(inout):: self
        real(fp), intent(in):: cutoff
        class(grasph_base_kernel), intent(in):: kernel

        if (self%is_pair_set) then
            call cell_list_search(self%ps_lhs%x, self%ps_rhs%x, self%ps_rhs%size, cutoff, kernel, self%pairs)
        else
            call cell_list_search(self%ps_lhs%x, cutoff, kernel, self%pairs)
        endif
    end subroutine particle_interactions_base_find_pairs

    subroutine particle_interactions_base_init(self, npairs_per_particle, ps_lhs, ps_rhs)
        class(particle_interactions_base), intent(out):: self
        integer, intent(in):: npairs_per_particle
        class(base_particles), target, intent(in):: ps_lhs
        class(base_particles), target, optional, intent(in):: ps_rhs
        self%ps_lhs => ps_lhs
        if (present(ps_rhs)) then
            self%ps_rhs => ps_rhs
            self%is_pair_set = .true.
        endif
        call self%pairs%init(ps_lhs%size, npairs_per_particle, ps_lhs%ndims)
    end subroutine particle_interactions_base_init

end module grasph_pair_sets
