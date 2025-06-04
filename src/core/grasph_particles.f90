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
        real(fp), allocatable:: x(:, :), v(:, :), rho(:), mass(:), c(:)
        real(fp), allocatable:: dvxdt(:, :), drhodt(:)
        type(particle_pairs):: pairs
        logical:: initialized = .false.
        integer:: ndims = 0, size = 0
    contains
        procedure:: base_init, base_clear, find_pairs, state_update => base_state_update
    end type base_particles

    type, extends(base_particles):: weakly_compressible_particles
        real(fp), allocatable:: p(:)
        real(fp):: rho_ref
    contains
        procedure:: init => wcp_init, state_update => linear_eos
    end type weakly_compressible_particles

    public:: base_particle, base_particles, weakly_compressible_particles

contains
    pure subroutine base_init(self, n, d, npairs_per_particle)
        class(base_particles), intent(inout):: self
        integer, intent(in):: n, d, npairs_per_particle
        if (self%initialized) call self%base_clear()
        allocate (self%id(n), self%type(n))
        allocate (self%x(d, n), self%v(d, n), self%rho(n), self%mass(n), self%c(n))
        allocate (self%dvxdt(d, n), self%drhodt(n))
        self%initialized = .true.
        self%size = n
        self%ndims = d
        call self%pairs%init(n, npairs_per_particle, d)
    end subroutine base_init

    pure subroutine base_clear(self)
        class(base_particles), intent(inout):: self
        if (self%initialized) then
            deallocate (self%id, self%type, self%x, self%v, self%rho, self%mass, self%c)
            deallocate (self%dvxdt, self%drhodt)
        endif
        self%initialized = .false.
        self%size = 0
    end subroutine base_clear

    subroutine find_pairs(self, kernel)
        class(base_particles), intent(inout):: self
        class(grasph_base_kernel), intent(in):: kernel
        call cell_list_search(self%x, kernel%cutoff, kernel, self%pairs)
    end subroutine find_pairs

    subroutine base_state_update(self, dt)
        class(base_particles), intent(inout):: self
        real(fp), intent(in), optional:: dt
        ! do nothing e.g. when using static repulsive boundaries that have no state
    end subroutine base_state_update

    subroutine wcp_init(self, n, d, npairs_per_particle, rho_ref)
        class(weakly_compressible_particles), intent(inout):: self
        integer, intent(in):: n, d, npairs_per_particle
        real(fp), intent(in):: rho_ref
        self%rho_ref = rho_ref
        if (self%initialized) deallocate (self%p)
        call self%base_init(n, d, npairs_per_particle)
        allocate (self%p(n))
    end subroutine wcp_init

    subroutine linear_eos(self, dt)
        class(weakly_compressible_particles), intent(inout):: self
        real(fp), intent(in), optional:: dt
        integer:: i
        do i = 1, self%size
            self%p(i) = self%c(i)**2*(self%rho(i) - self%rho_ref)
        enddo
    end subroutine linear_eos
end module grasph_particles
