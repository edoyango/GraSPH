module grasph_pair_sets

    use grasph_constants, only: fp, ndims
    use grasph_particles, only: base_particles
    use grasph_pairs, only: particle_pairs, cell_list_search
    use grasph_kernels, only: grasph_base_kernel

    implicit none

    private

    type, abstract:: particle_interactions_base
        type(particle_pairs):: pairs
        logical:: initialized = .false.
    contains
        procedure:: sweep
        procedure(find_pairs_interface), deferred:: find_pairs
    end type particle_interactions_base

    interface
        subroutine find_pairs_interface(self, cutoff, kernel)
            import:: particle_interactions_base, fp, grasph_base_kernel
            class(particle_interactions_base), intent(inout):: self
            real(fp), intent(in):: cutoff
            class(grasph_base_kernel), intent(in):: kernel
        end subroutine find_pairs_interface
    end interface

    public:: particle_interactions_base

contains

    subroutine sweep(self)
        class(particle_interactions_base), intent(inout):: self
    end subroutine sweep

end module grasph_pair_sets
