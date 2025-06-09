!> @file grasph_particles.f90
!> @brief Module containing core particles derived types and methods
!> @author Edward Yang
!> @date 2025-06-09
module grasph_particles

    use grasph_constants, only: fp
    use grasph_pairs, only: particle_pairs, cell_list_search
    use grasph_kernels, only: grasph_base_kernel

    implicit none
    private

    ! declaring both strcture (to be used in an array) and structure of arrays for testing
    ! type, abstract:: base_particle
    !     integer:: id, type
    !     real(fp):: x(ndims), v(ndims), rho, mass
    ! end type base_particle

    !> @brief The core particles derived type
    type:: base_particles
        integer, allocatable:: id(:), type(:)
        real(fp), allocatable:: x(:, :), v(:, :), rho(:), mass(:), c(:)
        real(fp), allocatable:: dvxdt(:, :), drhodt(:), v0(:, :), rho0(:) ! time-integration related data
        logical:: initialized = .false., to_print_summary = .true.
        integer:: ndims = 0, size = 0
        character(100):: name ! used in naming groups in output hdf5 file
    contains
        procedure:: base_init, base_clear
        procedure:: state_update => base_state_update
        procedure:: start_timestep => base_timestep_start, mid_timestep_update => base_midtimestep_update, &
                    full_timestep_update => base_fulltimestep_update
        procedure:: dump => base_dump, read => base_read
        procedure:: generate_summary => base_generate_summary
    end type base_particles

    !> @brief particles container class for setting up simulation
    type:: particles_container
        class(base_particles), allocatable:: p
    end type particles_container

    !> @brief particle type which adds pressure, determined from density with a linear EOS
    type, extends(base_particles):: weakly_compressible_particles
        real(fp), allocatable:: p(:)
        real(fp):: rho_ref
    contains
        procedure:: init => wcp_init, state_update => linear_eos, dump => wcp_dump, read => wcp_read
    end type weakly_compressible_particles

    public:: base_particles, weakly_compressible_particles, particles_container

contains
    !> @brief Initializes base_particles' internal arrays.
    !> @param self The particles to initialize.
    !> @param n The number of particles.
    !> @param d The dimension of the problem (1-3).
    !> @param name A label to give the particles. Used to label output/terminal information.
    subroutine base_init(self, n, d, name)
        class(base_particles), intent(inout):: self
        integer, intent(in):: n, d
        character(*), intent(in):: name
        if (self%initialized) call self%base_clear()
        allocate (self%id(n), self%type(n))
        allocate (self%x(d, n), self%v(d, n), self%rho(n), self%mass(n), self%c(n), source=0._fp)
        allocate (self%dvxdt(d, n), self%drhodt(n), self%v0(d, n), self%rho0(n), source=0._fp)
        self%initialized = .true.
        self%size = n
        self%ndims = d
        self%name = name
    end subroutine base_init

    !> @brief Deallocates internal arrays of self and sets state to uninitialized
    !> @param self The particles to clear
    subroutine base_clear(self)
        class(base_particles), intent(out):: self
    end subroutine base_clear

    !> @brief A do-nothing placeholder subroutine used in time-integration. Extend this with particles' internal state update code e.g. updating pressure, stress.
    !> @param self Particles whose state is to be updated.
    !> @param dt A time-increment which may be used to update particles' state.
    subroutine base_state_update(self, dt)
        class(base_particles), intent(inout):: self
        real(fp), intent(in), optional:: dt
        ! do nothing e.g. when using static repulsive boundaries that have no state
    end subroutine base_state_update

    !> @brief Subroutine to perform setup at start of every time-step.
    !> @param self Particles to setup.
    subroutine base_timestep_start(self)
        class(base_particles), intent(inout):: self
        ! save current velocity/density
        self%v0(:, :) = self%v(:, :)
        self%rho0(:) = self%rho(:)
    end subroutine base_timestep_start

    !> @brief Subroutine to update time-evolving variables.
    !> @param self Particles to evolve.
    !> @param dt Time increment.
    subroutine base_midtimestep_update(self, dt)
        class(base_particles), intent(inout):: self
        real(fp), intent(in):: dt
        self%v(:, :) = self%v(:, :) + dt*self%dvxdt(:, :)
        self%rho(:) = self%rho(:) + dt*self%drhodt(:)
    end subroutine base_midtimestep_update

    !> @brief Subroutine to update time-evolving variables using start-of-timestep data.
    !> @param self Particles to evolve.
    !> @param dt Time increment.
    !> @param update_position Whether to update position
    subroutine base_fulltimestep_update(self, dt, update_position)
        class(base_particles), intent(inout):: self
        real(fp), intent(in):: dt
        logical, intent(in):: update_position
        self%v(:, :) = self%v0(:, :) + dt*self%dvxdt(:, :)
        self%rho(:) = self%rho0(:) + dt*self%drhodt(:)
        if (update_position) self%x(:, :) = self%x(:, :) + dt*self%v(:, :)
    end subroutine base_fulltimestep_update

    !> @brief Writes base particle data to HDF5 file.
    !> @param self The particles to write.
    !> @param itimestep The timestep to use to encode in the filename.
    !> @param path The output directory.
    !> @param prefix The prefix to give to the output filenames.
    !> @param comp_level The level of gzip compression to use.
    subroutine base_dump(self, itimestep, path, prefix_in, comp_level)
        use h5fortran, only: hdf5_file
        class(base_particles), intent(in):: self
        integer, intent(in):: itimestep
        character(*), intent(in):: path
        character(*), intent(in), optional:: prefix_in
        integer, intent(in), optional:: comp_level
        character(*), parameter:: group = "base/"
        character(200):: filename_prefix, file_path, this_group
        integer:: ierr
        type(hdf5_file):: h5f
        character(10):: ic

        if (present(prefix_in)) then
            filename_prefix = prefix_in
        else
            filename_prefix = "grasph_particles"
        endif

        write(ic, "(I10.10)") itimestep
        file_path = path // "/" // trim(filename_prefix) // "_" // ic // ".h5"
        this_group = "/" // trim(self%name) // "/" // group

        call h5f%open(file_path, action="a", comp_lvl = comp_level)
        call h5f%write("/" // trim(self%name) // "/n", self%size)
        call h5f%write("/" // trim(self%name) // "/ndims", self%ndims)
        call h5f%write(trim(this_group) // "id", self%id)
        call h5f%write(trim(this_group) // "type", self%type)
        call h5f%write(trim(this_group) // "x", self%x)
        call h5f%write(trim(this_group) // "v", self%v)
        call h5f%write(trim(this_group) // "rho", self%rho)
        call h5f%write(trim(this_group) // "mass", self%mass)
        call h5f%write(trim(this_group) // "c", self%c)
        call h5f%write(trim(this_group) // "dvxdt", self%dvxdt)
        call h5f%write(trim(this_group) // "drhodt", self%drhodt)
        call h5f%write(trim(this_group) // "v0", self%v0)
        call h5f%write(trim(this_group) // "rho0", self%rho0)
        call h5f%close()

    end subroutine base_dump

    !> @brief Reads base particle data from HDF5 file.
    !> @param self The particles to read data into.
    !> @param file_path The path to the file to read.
    !> @param name The name to of particles to read and assign to the read particles.
    subroutine base_read(self, file_path, name)
        use h5fortran, only: hdf5_file
        class(base_particles), intent(out):: self
        character(*), intent(in):: name, file_path
        character(*), parameter:: group = "base/"
        character(200):: this_group
        integer:: d, n
        type(hdf5_file):: h5f

        this_group = "/" // trim(name) // "/" // group

        call h5f%open(file_path, action="r")
        call h5f%read("/" // trim(name) // "/n", n)
        call h5f%read("/" // trim(name) // "/ndims", d)
        call self%base_init(n, d, name)
        call h5f%read(trim(this_group) // "id", self%id)
        call h5f%read(trim(this_group) // "type", self%type)
        call h5f%read(trim(this_group) // "x", self%x)
        call h5f%read(trim(this_group) // "v", self%v)
        call h5f%read(trim(this_group) // "rho", self%rho)
        call h5f%read(trim(this_group) // "mass", self%mass)
        call h5f%read(trim(this_group) // "c", self%c)
        call h5f%read(trim(this_group) // "dvxdt", self%dvxdt)
        call h5f%read(trim(this_group) // "drhodt", self%drhodt)
        call h5f%read(trim(this_group) // "v0", self%v0)
        call h5f%read(trim(this_group) // "rho0", self%rho0)
        call h5f%close()

    end subroutine base_read

    subroutine base_generate_summary(self, out_str)

        class(base_particles), intent(in):: self
        character(len=:), allocatable, intent(out):: out_str
        real(fp):: val
        integer:: offset, d, i
        integer, parameter:: line_length = 60, nlines = 5
        character(*), parameter:: format_str = "(4x, A, f12.5, A, I10)"

        allocate(character(nlines*line_length)::out_str)
        offset = 0

        ! save max accel
        i = maxloc(sum(self%dvxdt(:, :)**2, dim=1), dim=1)
        val = sqrt(sum(self%dvxdt(:, i)**2))
        write(out_str(offset+1:offset+line_length), format_str) "  max(|dvdt|) of ", val, " at particle ", i
        offset = offset + line_length
        write(out_str(offset:offset), "(A1)") new_line("a")

        ! save max vel
        i = maxloc(sum(self%v(:, :)**2, dim=1), dim=1)
        val = sqrt(sum(self%v(:, i)**2))
        write(out_str(offset+1:offset+line_length), format_str) "     max(|v|) of ", val, " at particle ", i
        offset = offset + line_length
        write(out_str(offset:offset), "(A1)") new_line("a")

        ! save min rho
        i = minloc(self%rho, dim=1)
        write(out_str(offset+1:offset+line_length), format_str) "     min(rho) of ", self%rho(i), " at particle ", i
        offset = offset + line_length
        write(out_str(offset:offset), "(A1)") new_line("a")

        ! save max rho
        i = maxloc(self%rho, dim=1)
        write(out_str(offset+1:offset+line_length), format_str) "     max(rho) of ", self%rho(i), " at particle ", i
        offset = offset + line_length
        write(out_str(offset:offset), "(A1)") new_line("a")
        i = maxloc(abs(self%drhodt), dim=1)
        write(out_str(offset+1:offset+line_length), format_str) "max(|drhodt|) of ", self%drhodt(i), " at particle ", i
        

    end subroutine base_generate_summary

    subroutine wcp_init(self, n, d, name, rho_ref)
        class(weakly_compressible_particles), intent(inout):: self
        integer, intent(in):: n, d
        real(fp), intent(in):: rho_ref
        character(*), intent(in):: name
        self%rho_ref = rho_ref
        if (self%initialized) deallocate (self%p)
        call self%base_init(n, d, name)
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

    subroutine wcp_dump(self, itimestep, path, prefix_in, comp_level)
        use h5fortran, only: hdf5_file
        class(weakly_compressible_particles), intent(in):: self
        integer, intent(in):: itimestep
        character(*), intent(in):: path
        character(*), intent(in), optional:: prefix_in
        integer, intent(in), optional:: comp_level
        character(*), parameter:: group = "weakly_compressible/"
        character(200):: filename_prefix, file_path, this_group
        integer:: ierr
        type(hdf5_file):: h5f
        character(10):: ic

        if (present(prefix_in)) then
            filename_prefix = prefix_in
        else
            filename_prefix = "grasph_particles"
        endif

        write(ic, "(I10.10)") itimestep
        file_path = path // "/" // trim(filename_prefix) // "_" // ic // ".h5"
        this_group = "/" // trim(self%name) // "/" // group

        call base_dump(self, itimestep, path, prefix_in, comp_level)
        call h5f%open(file_path, action="a")
        call h5f%write(trim(this_group) // "p", self%p)
        call h5f%close()

    end subroutine wcp_dump

    subroutine wcp_read(self, file_path, name)
        use h5fortran, only: hdf5_file
        class(weakly_compressible_particles), intent(out):: self
        character(*), intent(in):: file_path, name
        character(*), parameter:: group = "weakly_compressible/"
        character(200):: filename, this_group
        integer:: ierr, d, n
        type(hdf5_file):: h5f
        character(10):: ic
        call base_read(self, file_path, name)
        
        this_group = "/" // trim(name) // "/" // group

        allocate(self%p(self%size))

        call h5f%open(file_path, action="r")
        call h5f%read(trim(this_group) // "p", self%p)
        call h5f%close()
    end subroutine wcp_read
end module grasph_particles
