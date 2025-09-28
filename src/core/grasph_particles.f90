!> @file grasph_particles.f90
!> @brief Module containing core particles derived types and methods
!> @author Edward Yang
!> @date 2025-06-09
module grasph_particles

    use grasph_constants, only: fp
    use grasph_common, only: array_pointer_container

    implicit none
    private

    ! declaring both strcture (to be used in an array) and structure of arrays for testing
    ! type, abstract:: base_particle
    !     integer:: id, type
    !     real(fp):: x(ndims), v(ndims), rho, mass
    ! end type base_particle

    !> @brief Max allowable registrations. Will be removed for a more dynamic approach.
    integer, parameter:: max_registrations = 20

    !> @brief Register for variables and their derivatives to be updated in time-integrations.
    type:: data_register
        !> @brief Running total number of registrations for this register.
        integer:: nregistrations = 0
        !> @brief Names of the variables being registered. First index is the name of the variable to be registered. The second is
        !>        the name of its derivative.
        character(50):: names(2, max_registrations)
        !> @brief Pointers to the array variables being registered. First index is the variable to be registered. The second is its
        !>        derivative.
        type(array_pointer_container):: data(2, max_registrations)
    contains
        procedure:: register_data_1d, register_data_2d
        generic:: register_data => register_data_1d, register_data_2d
    end type data_register

    !> @brief The core particles derived type
    type:: base_particles
        !> @brief The id of the particles.
        integer, allocatable:: id(:)
        !> @brief The type of the particles.
        integer, allocatable:: type(:)
        !> @brief The particles' position.
        real(fp), allocatable:: x(:, :)
        !> @brief The particles' velocity.
        real(fp), allocatable:: v(:, :)
        !> @brief The particles' density.
        real(fp), allocatable:: rho(:)
        !> @brief The particles' mass.
        real(fp), allocatable:: mass(:)
        !> @brief The particles' local speed of sound.
        real(fp), allocatable:: c(:)
        !> @brief The particles' acceleration
        real(fp), allocatable:: dvxdt(:, :)
        !> @brief The particles' density rate-of-change.
        real(fp), allocatable:: drhodt(:)
        !> @brief The particles' velocity at the start of a time-step.
        real(fp), allocatable:: v0(:, :)
        !> @brief The particles' density at the start of a time-step
        real(fp), allocatable:: rho0(:)
        !> @brief Whether the particles have been initialized.
        logical:: initialized = .false.
        !> @brief Whether to update the particles' properties
        logical:: evolve = .true.
        !> @brief Whether to print information when generate_summary is called.
        logical:: to_print_summary = .true.
        !> @brief Number of spatial dimensions
        integer:: ndims = 0
        !> @brief Number of particles.
        integer:: size = 0
        !> @brief Name used in naming groups in output hdf5 file
        character(100):: name
        !> @brief Register for variables to be updated only at full-timestep e.g. position (x).
        type(data_register):: register_x
        !> @brief Register for variables to be updated at mid- and full-timestep e.g. velocity (v) and density (rho).
        type(data_register):: register_v
    contains
        !> @brief The initializer for the base class. Intended to be called in extended types' initializer method.
        procedure:: base_init
        !> @brief A manual destructor to clean up.
        procedure:: base_clear
        !> @brief A method intended to be overriden when extended particles' state needs to be updated during time-integration.
        !>        Is called after any time-evolution has occurred, but before any sweeps are supposed to happen.
        !>        Does nothing in the base_particles instance.
        procedure:: state_update => base_state_update
        !> @brief A method intended to be overriden when extended particles' have extra data that needs to be saved in the output
        !>        files. Different derived types should store their data in different groups
        !>        e.g. /\<name>/base_particles/..., and /\<name>/derived_particles/....
        procedure:: dump => base_dump
        !> @brief A method intended to be overriden when extended particles' have extra data that needs to be read from the output
        !>        files.
        procedure:: read => base_read
        !> @brief A method that can be overriden to generate a string with particles' summary data to print during time-integration.
        procedure:: generate_summary => base_generate_summary
    end type base_particles

    !> @brief particles container class for setting up simulation
    type:: particles_container
        !> @brief The polymorphic container to be allocated to base_particles or its derivatives.
        class(base_particles), allocatable:: p
    end type particles_container

    public:: base_particles, particles_container, base_dump, base_read, max_registrations

contains

    !> @brief Registers 1D variables e.g. density (rho).
    !> @param self The register to register the variable into.
    !> @param var The variable to register.
    !> @param var_name The name of the variable being registered.
    !> @param deriv The derivative variable of the variable to register.
    !> @param deriv_name The name of the derivative variable being registered.
    subroutine register_data_1d(self, var, var_name, deriv, deriv_name)
        class(data_register), intent(inout):: self
        real(fp), target, intent(in):: var(:), deriv(:)
        character(*), intent(in):: var_name, deriv_name

        if (self%nregistrations == max_registrations) error stop "can't register anymore"
        self%nregistrations = self%nregistrations + 1
        self%names(1, self%nregistrations) = var_name
        self%names(2, self%nregistrations) = deriv_name
        self%data(1, self%nregistrations)%p(1:1, 1:size(var)) => var
        self%data(2, self%nregistrations)%p(1:1, 1:size(deriv)) => deriv

    end subroutine register_data_1d

    !> @brief Registers 2D variables e.g. position (x).
    !> @param self The register to register the variable into.
    !> @param var The variable to register.
    !> @param var_name The name of the variable being registered.
    !> @param deriv The derivative variable of the variable to register.
    !> @param deriv_name The name of the derivative variable being registered.
    subroutine register_data_2d(self, var, var_name, deriv, deriv_name)
        class(data_register), intent(inout):: self
        real(fp), target, intent(in):: var(:, :), deriv(:, :)
        character(*), intent(in):: var_name, deriv_name

        if (self%nregistrations == max_registrations) error stop "can't register anymore"
        self%nregistrations = self%nregistrations + 1
        self%names(1, self%nregistrations) = var_name
        self%names(2, self%nregistrations) = deriv_name
        self%data(1, self%nregistrations)%p => var
        self%data(2, self%nregistrations)%p => deriv

    end subroutine register_data_2d

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

    !> @brief Writes base particle data to HDF5 file.
    !> @param self The particles to write.
    !> @param itimestep The timestep to add to the filename.
    !> @param path The output directory.
    !> @param prefix_in The prefix to give to the output filenames.
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
        end if

        write (ic, "(I10.10)") itimestep
        file_path = path//"/"//trim(filename_prefix)//"_"//ic//".h5"
        this_group = "/"//trim(self%name)//"/"//group

        call h5f%open(file_path, action="a", comp_lvl=comp_level)
        call h5f%write("/"//trim(self%name)//"/n", self%size)
        call h5f%write("/"//trim(self%name)//"/ndims", self%ndims)
        call h5f%write(trim(this_group)//"id", self%id)
        call h5f%write(trim(this_group)//"type", self%type)
        call h5f%write(trim(this_group)//"x", self%x)
        call h5f%write(trim(this_group)//"v", self%v)
        call h5f%write(trim(this_group)//"rho", self%rho)
        call h5f%write(trim(this_group)//"mass", self%mass)
        call h5f%write(trim(this_group)//"c", self%c)
        call h5f%write(trim(this_group)//"dvxdt", self%dvxdt)
        call h5f%write(trim(this_group)//"drhodt", self%drhodt)
        call h5f%write(trim(this_group)//"v0", self%v0)
        call h5f%write(trim(this_group)//"rho0", self%rho0)
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

        this_group = "/"//trim(name)//"/"//group

        call h5f%open(file_path, action="r")
        call h5f%read("/"//trim(name)//"/n", n)
        call h5f%read("/"//trim(name)//"/ndims", d)
        call self%base_init(n, d, name)
        call h5f%read(trim(this_group)//"id", self%id)
        call h5f%read(trim(this_group)//"type", self%type)
        call h5f%read(trim(this_group)//"x", self%x)
        call h5f%read(trim(this_group)//"v", self%v)
        call h5f%read(trim(this_group)//"rho", self%rho)
        call h5f%read(trim(this_group)//"mass", self%mass)
        call h5f%read(trim(this_group)//"c", self%c)
        call h5f%read(trim(this_group)//"dvxdt", self%dvxdt)
        call h5f%read(trim(this_group)//"drhodt", self%drhodt)
        call h5f%read(trim(this_group)//"v0", self%v0)
        call h5f%read(trim(this_group)//"rho0", self%rho0)
        call h5f%close()

    end subroutine base_read

    !> @brief Controls particles' summary stats.
    !> @param self The particles to print the stats of.
    !> @param out_str The string which contains the summary stats and any formatting.
    subroutine base_generate_summary(self, out_str)

        class(base_particles), intent(in):: self
        character(len=:), allocatable, intent(out):: out_str
        real(fp):: val
        integer:: offset, d, i
        integer, parameter:: line_length = 60, nlines = 5
        character(*), parameter:: format_str = "(4x, A, f12.5, A, I10)"

        allocate (character(nlines*line_length)::out_str)
        offset = 0

        ! save max accel
        i = maxloc(sum(self%dvxdt(:, :)**2, dim=1), dim=1)
        val = sqrt(sum(self%dvxdt(:, i)**2))
        write (out_str(offset + 1:offset + line_length), format_str) "  max(|dvdt|) of ", val, " at particle ", i
        offset = offset + line_length
        write (out_str(offset:offset), "(A1)") new_line("a")

        ! save max vel
        i = maxloc(sum(self%v(:, :)**2, dim=1), dim=1)
        val = sqrt(sum(self%v(:, i)**2))
        write (out_str(offset + 1:offset + line_length), format_str) "     max(|v|) of ", val, " at particle ", i
        offset = offset + line_length
        write (out_str(offset:offset), "(A1)") new_line("a")

        ! save min rho
        i = minloc(self%rho, dim=1)
        write (out_str(offset + 1:offset + line_length), format_str) "     min(rho) of ", self%rho(i), " at particle ", i
        offset = offset + line_length
        write (out_str(offset:offset), "(A1)") new_line("a")

        ! save max rho
        i = maxloc(self%rho, dim=1)
        write (out_str(offset + 1:offset + line_length), format_str) "     max(rho) of ", self%rho(i), " at particle ", i
        offset = offset + line_length
        write (out_str(offset:offset), "(A1)") new_line("a")
        i = maxloc(abs(self%drhodt), dim=1)
        write (out_str(offset + 1:offset + line_length), format_str) "max(|drhodt|) of ", self%drhodt(i), " at particle ", i

    end subroutine base_generate_summary

end module grasph_particles
