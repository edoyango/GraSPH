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

    type:: base_particles
        integer, allocatable:: id(:), type(:)
        real(fp), allocatable:: x(:, :), v(:, :), rho(:), mass(:), c(:)
        real(fp), allocatable:: dvxdt(:, :), drhodt(:), v0(:, :), rho0(:) ! time-integration related data
        logical:: initialized = .false.
        integer:: ndims = 0, size = 0
    contains
        procedure:: base_init, base_clear
        procedure:: state_update => base_state_update
        procedure:: start_timestep => base_timestep_start, mid_timestep_update => base_midtimestep_update, &
                    full_timestep_update => base_fulltimestep_update
        procedure:: dump => base_dump, read => base_read
    end type base_particles

    type:: particles_container
        class(base_particles), allocatable:: p
    end type particles_container

    type, extends(base_particles):: weakly_compressible_particles
        real(fp), allocatable:: p(:)
        real(fp):: rho_ref
    contains
        procedure:: init => wcp_init, state_update => linear_eos, dump => wcp_dump, read => wcp_read
    end type weakly_compressible_particles

    public:: base_particle, base_particles, weakly_compressible_particles, particles_container

contains
    pure subroutine base_init(self, n, d)
        class(base_particles), intent(inout):: self
        integer, intent(in):: n, d
        if (self%initialized) call self%base_clear()
        allocate (self%id(n), self%type(n))
        allocate (self%x(d, n), self%v(d, n), self%rho(n), self%mass(n), self%c(n))
        allocate (self%dvxdt(d, n), self%drhodt(n), self%v0(d, n), self%rho0(n))
        self%initialized = .true.
        self%size = n
        self%ndims = d
    end subroutine base_init

    pure subroutine base_clear(self)
        class(base_particles), intent(inout):: self
        if (self%initialized) then
            deallocate (self%id, self%type, self%x, self%v, self%rho, self%mass, self%c)
            deallocate (self%dvxdt, self%drhodt, self%v0, self%rho0)
        endif
        self%initialized = .false.
        self%size = 0
    end subroutine base_clear

    subroutine base_state_update(self, dt)
        class(base_particles), intent(inout):: self
        real(fp), intent(in), optional:: dt
        ! do nothing e.g. when using static repulsive boundaries that have no state
    end subroutine base_state_update

    subroutine base_timestep_start(self)
        class(base_particles), intent(inout):: self
        self%v0(:, :) = self%v(:, :)
        self%rho0(:) = self%rho0(:)
    end subroutine base_timestep_start

    subroutine base_midtimestep_update(self, dt)
        class(base_particles), intent(inout):: self
        real(fp), intent(in):: dt
        self%v(:, :) = self%v(:, :) + dt*self%dvxdt(:, :)
        self%rho(:) = self%rho(:) + dt*self%drhodt(:)
    end subroutine base_midtimestep_update

    subroutine base_fulltimestep_update(self, dt, update_position)
        class(base_particles), intent(inout):: self
        real(fp), intent(in):: dt
        logical, intent(in):: update_position
        self%v(:, :) = self%v0(:, :) + dt*self%dvxdt(:, :)
        self%rho(:) = self%rho0(:) + dt*self%drhodt(:)
        if (update_position) self%x(:, :) = self%x(:, :) + dt*self%v(:, :)
    end subroutine base_fulltimestep_update

    subroutine base_dump(self, itimestep, path, prefix, comp_level)
        use h5fortran, only: hdf5_file
        class(base_particles), intent(in):: self
        integer, intent(in):: itimestep
        character(*), intent(in):: path, prefix
        integer, intent(in), optional:: comp_level
        character(*), parameter:: group = "/base"
        character(200):: filename
        integer:: ierr
        type(hdf5_file):: h5f
        character(10):: ic

        write(ic, "(I10.10)") itimestep
        filename = path // "/" // prefix // "grasph_particles_" // ic // ".h5"

        call h5f%open(filename, action="w", comp_lvl = comp_level)
        call h5f%write("/n", self%size)
        call h5f%write("/ndims", self%ndims)
        call h5f%write(group // "/id", self%id)
        call h5f%write(group // "/type", self%type)
        call h5f%write(group // "/x", self%x)
        call h5f%write(group // "/v", self%v)
        call h5f%write(group // "/rho", self%rho)
        call h5f%write(group // "/mass", self%mass)
        call h5f%write(group // "/c", self%c)
        call h5f%write(group // "/dvxdt", self%dvxdt)
        call h5f%write(group // "/drhodt", self%drhodt)
        call h5f%write(group // "/v0", self%v0)
        call h5f%write(group // "/rho0", self%rho0)
        call h5f%close()

    end subroutine base_dump

    subroutine base_read(self, itimestep, path, prefix)
        use h5fortran, only: hdf5_file
        class(base_particles), intent(out):: self
        integer, intent(in):: itimestep
        character(*), intent(in):: path, prefix
        character(*), parameter:: group = "/base"
        character(200):: filename
        integer:: ierr, d, n
        type(hdf5_file):: h5f
        character(10):: ic

        write(ic, "(I10.10)") itimestep
        filename = path // "/" // prefix // "grasph_particles_" // ic // ".h5"

        call h5f%open(filename, action="r")
        call h5f%read("/n", n)
        call h5f%read("/ndims", d)
        call self%base_init(n, d)
        call h5f%read(group // "/id", self%id)
        call h5f%read(group // "/type", self%type)
        call h5f%read(group // "/x", self%x)
        call h5f%read(group // "/v", self%v)
        call h5f%read(group // "/rho", self%rho)
        call h5f%read(group // "/mass", self%mass)
        call h5f%read(group // "/c", self%c)
        call h5f%read(group // "/dvxdt", self%dvxdt)
        call h5f%read(group // "/drhodt", self%drhodt)
        call h5f%read(group // "/v0", self%v0)
        call h5f%read(group // "/rho0", self%rho0)
        call h5f%close()

    end subroutine base_read

    subroutine wcp_init(self, n, d, rho_ref)
        class(weakly_compressible_particles), intent(inout):: self
        integer, intent(in):: n, d
        real(fp), intent(in):: rho_ref
        self%rho_ref = rho_ref
        if (self%initialized) deallocate (self%p)
        call self%base_init(n, d)
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

    subroutine wcp_dump(self, itimestep, path, prefix, comp_level)
        use h5fortran, only: hdf5_file
        class(weakly_compressible_particles), intent(in):: self
        integer, intent(in):: itimestep
        character(*), intent(in):: path, prefix
        integer, intent(in), optional:: comp_level
        character(*), parameter:: group = "/weakly_compressible"
        character(200):: filename
        integer:: ierr
        type(hdf5_file):: h5f
        character(10):: ic

        write(ic, "(I10.10)") itimestep
        filename = path // "/" // prefix // "grasph_particles_" // ic // ".h5"

        call base_dump(self, itimestep, path, prefix, comp_level)
        call h5f%open(filename, action="a")
        call h5f%write(group // "/p", self%p)
        call h5f%close()

    end subroutine wcp_dump

    subroutine wcp_read(self, itimestep, path, prefix)
        use h5fortran, only: hdf5_file
        class(weakly_compressible_particles), intent(out):: self
        integer, intent(in):: itimestep
        character(*), intent(in):: path, prefix
        character(*), parameter:: group = "/weakly_compressible"
        character(200):: filename
        integer:: ierr, d, n
        type(hdf5_file):: h5f
        character(10):: ic
        call base_read(self, itimestep, path, prefix)
        
        write(ic, "(I10.10)") itimestep
        filename = path // "/" // prefix // "grasph_particles_" // ic // ".h5"

        call h5f%open(filename, action="r")
        call h5f%read(group // "/p", self%p)
        call h5f%close()
    end subroutine wcp_read
end module grasph_particles
