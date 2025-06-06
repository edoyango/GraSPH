module grasph_pair_sets

    use grasph_constants, only: fp, ndims
    use grasph_particles, only: base_particles
    use grasph_pairs, only: particle_pairs, cell_list_search
    use grasph_kernels, only: grasph_base_kernel

    implicit none

    private

    type, abstract:: interacting_particle_set
        class(base_particles), pointer:: lhs_particles => null(), rhs_particles => null()
        type(particle_pairs):: pairs
        logical:: is_pair_set = .false.
    contains
        procedure:: find_pairs => particle_set_find_pairs
        generic:: base_init => particle_pair_set_base_init, particle_single_set_base_init
        procedure:: particle_pair_set_base_init, particle_single_set_base_init
        procedure:: sweep
        procedure(pair_update_interface), deferred:: pair_update
    end type interacting_particle_set

    interface
        subroutine pair_update_interface(self, i, j)
            import:: interacting_particle_set
            class(interacting_particle_set), intent(inout):: self
            integer, intent(in):: i, j
        end subroutine pair_update_interface
    end interface

    public:: interacting_particle_set

contains

    subroutine particle_pair_set_base_init(self, lhs, rhs, npairs_per_particle)
        class(interacting_particle_set), intent(inout):: self
        class(base_particles), intent(in), target:: lhs, rhs
        integer, intent(in):: npairs_per_particle

        self%lhs_particles => lhs
        self%rhs_particles => rhs
        call self%pairs%init(lhs%size, npairs_per_particle, lhs%ndims)
        self%is_pair_set = .true.

    end subroutine particle_pair_set_base_init

    subroutine particle_single_set_base_init(self, ps, npairs_per_particle)
        class(interacting_particle_set), intent(inout):: self
        class(base_particles), intent(in), target:: ps
        integer, intent(in):: npairs_per_particle

        self%lhs_particles => ps
        self%rhs_particles => null()
        call self%pairs%init(ps%size, npairs_per_particle, ps%ndims)
        self%is_pair_set = .false.
    end subroutine particle_single_set_base_init

    subroutine particle_set_find_pairs(self, cutoff, kernel)
        class(interacting_particle_set), intent(inout):: self
        real(fp), intent(in):: cutoff
        class(grasph_base_kernel), intent(in):: kernel

        if (self%is_pair_set) then
            call cell_list_search( &
                self%lhs_particles%x, &
                self%rhs_particles%x, &
                self%rhs_particles%size, &
                cutoff, &
                kernel, &
                self%pairs &
            )
        else
            call cell_list_search( &
                self%lhs_particles%x, &
                cutoff, &
                kernel, &
                self%pairs &
            )
        endif
    end subroutine particle_set_find_pairs

    subroutine sweep(self)
        class(interacting_particle_set), intent(inout):: self
        integer:: i, jj, j
        if (self%pairs%initialized) then
            do i = 1, self%pairs%n
                do jj = self%pairs%offsets(i)+1, self%pairs%offsets(i+1)
                    j = self%pairs%rhs(jj)
                    call self%pair_update(i, j)
                enddo
            enddo
        endif
    end subroutine sweep

end module grasph_pair_sets
