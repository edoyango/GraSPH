module grasph_particles

    use grasph_constants, only: fp, ndims
    use grasph_pairs, only: particle_pairs, cell_list_search
    use grasph_kernels, only: grasph_base_kernel

    implicit none
    private

    ! declaring both strcture (to be used in an array) and structure of arrays for testing
    type, abstract:: base_particle
        integer:: id, type
        real(fp):: x(ndims), v(ndims), rho, mass
    end type base_particle

    type, abstract:: base_particles
        integer, allocatable:: id(:), type(:)
        real(fp), allocatable:: x(:, :), v(:, :), rho(:), mass(:)
        type(particle_pairs):: pairs
        logical:: initialized = .false.
        integer:: ndims = 0, size = 0
    contains
        procedure:: base_init, base_clear, find_pairs
    end type base_particles

    type, extends(base_particles):: weakly_compressible_particles
        real(fp), allocatable:: p(:)
    contains
        procedure:: init => wcp_init
    end type weakly_compressible_particles

    public:: base_particle, base_particles, weakly_compressible_particles

contains
    pure subroutine base_init(self, n, d, npairs_per_particle)
        class(base_particles), intent(inout):: self
        integer, intent(in):: n, d, npairs_per_particle
        if (self%initialized) call self%base_clear()
        allocate(self%id(n), self%type(n))
        allocate(self%x(d, n), self%v(d, n), self%rho(n), self%mass(n))
        self%initialized = .true.
        self%size = n
        self%ndims = d
        call self%pairs%init(n, npairs_per_particle, d)
    end subroutine base_init

    pure subroutine base_clear(self)
        class(base_particles), intent(inout):: self
        if (self%initialized) deallocate(self%id, self%type, self%x, self%v, self%rho, self%mass)
        self%initialized = .false.
        self%size = 0
    end subroutine base_clear

    subroutine find_pairs(self, kernel)
        class(base_particles), intent(inout):: self
        class(grasph_base_kernel), intent(in):: kernel
        call cell_list_search(self%x, kernel%cutoff, kernel, self%pairs)
    end subroutine find_pairs

    subroutine wcp_init(self, n, d, npairs_per_particle)
        class(weakly_compressible_particles), intent(inout):: self
        integer, intent(in):: n, d, npairs_per_particle
        if (self%initialized) deallocate(self%p)
        call self%base_init(n, d, npairs_per_particle)
        allocate(self%p(n))
    end subroutine wcp_init
end module grasph_particles
