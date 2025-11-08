!> @file grasph_particle_system.f90
!> @brief Module containing core particles derived types and methods
!> @author Edward Yang
!> @date 2025-06-09
module grasph_particle_system_m

    use iso_fortran_env, only: error_unit
    use grasph_constants_m, only: fp, ndims
    use grasph_particle_m, only: base_particle_t
    use grasph_common_m, only: array_pointer_container_t
    use grasph_register_m, only: variable_register_t, variable_deriv_register_t

    implicit none
    private

    !> @brief State updater which does nothing when update_state method is called.
    type:: base_state_updater_t
    contains
        !> @brief A do-nothing state updater.
        procedure:: update_state => base_update_state
    end type base_state_updater_t

    type:: state_updater_container_t
        class(base_state_updater_t), allocatable:: updater
    end type state_updater_container_t

    !> @brief Manages a group of particles that behave similarly.
    type:: particle_system_t
        !> @brief The particles that comprise the system.
        class(base_particle_t), allocatable:: particles(:)
        !> @brief Whether the system have been initialized.
        logical:: initialized = .false.
        !> @brief Whether to print information when generate_summary is called.
        logical:: to_print_summary = .true.
        !> @brief Number of spatial dimensions
        integer:: ndims = ndims
        !> @brief Number of particles.
        integer:: size = 0
        !> @brief Name used in naming groups in output hdf5 file
        character(100):: name
        !> @brief Allocatable "strategy" class that performs particles' first state update.
        class(state_updater_container_t), allocatable:: state_updaters(:)
        !> @brief Register for variables to be updated only at full-timestep e.g. position (x).
        type(variable_deriv_register_t):: register_x
        !> @brief Register for variables to be updated at mid- and full-timestep e.g. velocity (v) and density (rho).
        type(variable_deriv_register_t):: register_v
        !> @brief Register for variables to be written/read.
        type(variable_register_t):: register_io
    contains
        !> @brief The initializer for the base class.
        procedure:: init => base_init
        !> @brief A manual destructor to clean up.
        procedure:: clear => base_clear
        !> @brief A method intended to be overriden when extended particles' state needs to be updated during time-integration.
        !>        Is called after any time-evolution has occurred, but before any sweeps are supposed to happen.
        procedure:: do_state_update
        !> @brief A method intended to be overriden when extended particles' have extra data that needs to be saved in the output
        !>        files. Different derived types should store their data in different groups
        !>        e.g. /\<name>/particle_system_t/..., and /\<name>/derived_particles/....
        procedure:: dump => base_dump
        !> @brief A method intended to be overriden when extended particles' have extra data that needs to be read from the output
        !>        files.
        procedure:: read => base_read
        !> @brief A method that can be overriden to generate a string with particles' summary data to print during time-integration.
        procedure:: generate_summary => base_generate_summary
        !> @brief Function that returns the next allocated index in self%particles. If there isn't enough space, self%particles
        !>        is resized.
        procedure:: safe_size_plus_1
    end type particle_system_t

    public:: base_particle_t, particle_system_t, base_state_updater_t, state_updater_container_t

contains

    !> @brief Initializes particle_system_t' internal arrays.
    !> @param self The particle system to initialize.
    !> @param n The number of particles.
    !> @param name A label to give the particles. Used to label output/terminal information.
    !> @param particle_template Template to use to define the type of particles.
    !> @param state_updater_1 The state updater strategy class used to update the particles' time-independent state.
    !>        Is called before the prologue sweep.
    !> @param state_updater_2 The state updater strategy class used to update the particles' time-independent state.
    !>        Is called before the main sweep.
    subroutine base_init(self, n, name, particle_template, state_updaters)
        class(particle_system_t), intent(inout):: self
        integer, intent(in):: n
        character(*), intent(in):: name
        class(base_particle_t), optional, intent(in):: particle_template
        class(state_updater_container_t), optional, intent(in):: state_updaters(:)

        if (self%initialized) call self%clear()
        if (present(particle_template)) then
            allocate (self%particles(n), source=particle_template)
        else
            allocate (self%particles(n))
        end if

        if (present(state_updaters)) then
            allocate (self%state_updaters, source=state_updaters)
        end if

        self%initialized = .true.
        self%size = n
        self%ndims = ndims
        self%name = name

        ! add default IO registrations
        call self%register_io%register_variable(self%particles(1), "x", self%particles(1)%x)
        call self%register_io%register_variable(self%particles(1), "v", self%particles(1)%v)
        call self%register_io%register_variable(self%particles(1), "rho", self%particles(1)%rho)
        call self%register_io%register_variable(self%particles(1), "mass", self%particles(1)%mass)
        call self%register_io%register_variable(self%particles(1), "c", self%particles(1)%c)
        call self%register_io%register_variable(self%particles(1), "dvxdt", self%particles(1)%dvxdt)
        call self%register_io%register_variable(self%particles(1), "drhodt", self%particles(1)%drhodt)

    end subroutine base_init

    !> @brief Deallocates internal arrays of self and sets state to uninitialized
    !> @param self The particles to clear
    subroutine base_clear(self)
        class(particle_system_t), intent(out):: self
    end subroutine base_clear

    subroutine do_state_update(self, i, dt)
        class(particle_system_t), intent(inout):: self
        integer, intent(in):: i
        real(fp), optional, intent(in):: dt
        integer:: ii

        call self%state_updaters(i)%updater%update_state(self%particles, self%size, dt)

    end subroutine do_state_update

    !> @brief A do-nothing placeholder subroutine used in time-integration. Extend this with particles' internal state update code e.g. updating pressure, stress.
    !> @param self Particles whose state is to be updated.
    !> @param ps Particles whose state is to be updated.
    !> @param n Number of particles in ps.
    !> @param dt A time-increment which may be used to update particles' state.
    subroutine base_update_state(self, ps, n, dt)
        class(base_state_updater_t), intent(in):: self
        integer, intent(in):: n
        class(base_particle_t), intent(inout):: ps(n)
        real(fp), intent(in), optional:: dt
        ! do nothing e.g. when using static repulsive boundaries that have no state
    end subroutine base_update_state

    !> @brief Method that adds 1 to self%size and returns the result, but ensures that self%particles at the resultant index is
    !>        allocated. However, the particle at the index is in an undefined state and should be updated manually.
    !> @param self THe particle_system_t to add one to the size of.
    integer function safe_size_plus_1(self)
        class(particle_system_t), intent(inout):: self
        class(base_particle_t), allocatable:: tmp_particle(:)

        if (.not. (allocated(self%particles)) .or. size(self%particles) == 0) &
            error stop "Cannot add to unallocated or zero-sized particles."

        if (self%size == size(self%particles)) then
            ! allocate tmp_particle to ensure it's same type and size as self%particles
            allocate (tmp_particle, mold=self%particles)
            self%particles = [self%particles, tmp_particle] ! this doubles the space in self%particles
        end if

        self%size = self%size + 1
        safe_size_plus_1 = self%size

    end function safe_size_plus_1

    !> @brief Writes a system's particle data to HDF5 file.
    !> @param self The particle system to write.
    !> @param itimestep The timestep to add to the filename.
    !> @param path The output directory.
    !> @param prefix_in The prefix to give to the output filenames.
    !> @param comp_level The level of gzip compression to use.
    subroutine base_dump(self, itimestep, path, prefix_in, comp_level)
        use h5fortran, only: hdf5_file
        class(particle_system_t), intent(in):: self
        integer, intent(in):: itimestep
        character(*), intent(in):: path
        character(*), intent(in), optional:: prefix_in
        integer, intent(in), optional:: comp_level
        character(200):: filename_prefix, file_path, this_group
        integer:: i, v
        type(hdf5_file):: h5f
        character(10):: ic
        integer, allocatable:: tmp_int(:)
        real(fp), allocatable:: tmp_real(:, :)
        real(fp), pointer:: ptr(:)

        if (present(prefix_in)) then
            filename_prefix = prefix_in
        else
            filename_prefix = "grasph_particles"
        end if

        write (ic, "(I10.10)") itimestep
        file_path = path//"/"//trim(filename_prefix)//"_"//ic//".h5"
        this_group = "/"//trim(self%name)//"/"

        call h5f%open(file_path, action="a", comp_lvl=comp_level)
        call h5f%write("/"//trim(self%name)//"/n", self%size)
        call h5f%write("/"//trim(self%name)//"/ndims", self%ndims)
        allocate (tmp_int(self%size))
        do i = 1, self%size
            tmp_int(i) = self%particles(i)%id
        end do
        call h5f%write(trim(this_group)//"id", tmp_int)
        do i = 1, self%size
            tmp_int(i) = self%particles(i)%type
        end do
        call h5f%write(trim(this_group)//"type", tmp_int)
        do v = 1, self%register_io%nregistrations
            if (allocated(tmp_real)) deallocate (tmp_real)
            if (self%register_io%dims(v) == 1) then
                allocate (tmp_real(self%size, 1))
                do i = 1, self%size
                    call self%register_io%get_variable(self%particles(i), v, ptr)
                    tmp_real(i, 1) = ptr(1)
                end do
                call h5f%write(trim(this_group)//trim(self%register_io%names(v)), tmp_real(:, 1))
            else
                allocate (tmp_real(self%register_io%dims(v), self%size))
                do i = 1, self%size
                    call self%register_io%get_variable(self%particles(i), v, ptr)
                    tmp_real(:, i) = ptr(:)
                end do
                call h5f%write(trim(this_group)//trim(self%register_io%names(v)), tmp_real)
            end if
        end do
        call h5f%close()

    end subroutine base_dump

    !> @brief Reads registered particle data from HDF5 file into the pre-initialized particle system.
    !> @param self The pre-initialized particle system to read data into.
    !> @param file_path The path to the file to read.
    !> @param name The name to of particles to read and assign to the read particles.
    subroutine base_read(self, file_path, name)
        use h5fortran, only: hdf5_file, hsize_t
        class(particle_system_t), intent(inout):: self
        character(*), intent(in):: name, file_path
        character(200):: this_group
        character(250):: arr_path
        integer:: d, n, i, v, nrank
        type(hdf5_file):: h5f
        integer, allocatable:: tmp_int(:)
        integer(hsize_t), allocatable:: dims(:)
        real(fp), allocatable:: tmp_real(:, :)
        real(fp), pointer:: ptr(:)
        character(2):: nc_dim_arr, nc_dim_h5

        this_group = "/"//trim(name)//"/"

        call h5f%open(file_path, action="r")
        call h5f%read("/"//trim(name)//"/n", n)
        call h5f%read("/"//trim(name)//"/ndims", d)
        if (d /= ndims) error stop "Input HDF5 file dimensions don't match code dimensions."
        allocate (tmp_int(n))
        call h5f%read(trim(this_group)//"id", tmp_int)
        do i = 1, n
            self%particles(i)%id = tmp_int(i)
        end do
        call h5f%read(trim(this_group)//"type", tmp_int)
        do i = 1, n
            self%particles(i)%type = tmp_int(i)
        end do
        deallocate (tmp_int)

        ! iterate over registered variables
        do v = 1, self%register_io%nregistrations
            ! name of array in hdf5 file
            arr_path = trim(this_group)//trim(self%register_io%names(v))

            ! inspect rank and shape of array in hdf5 file
            nrank = h5f%ndim(arr_path)
            call h5f%shape(arr_path, dims)

            ! handle storing data base on array rank
            if (nrank == 1) then
                if (self%register_io%dims(v) /= 1) &
                    error stop "Expected rank 1 array for "//trim(arr_path)//"in input HDF5 file, "//file_path//"."
                allocate (tmp_real(n, 1))
                call h5f%read(arr_path, tmp_real(:, 1))
                do i = 1, n
                    call self%register_io%get_variable(self%particles(i), v, ptr)
                    ptr(1) = tmp_real(i, 1)
                end do
            elseif (nrank == 2) then
                if (self%register_io%dims(v) /= dims(1)) then
                    write (nc_dim_arr, "(I2)") self%register_io%dims(v)
                    write (nc_dim_h5, "(I2)") dims(1)
                    error stop "Expected dim 1 of "//trim(arr_path)//"in input HDF5 file to be "//trim(nc_dim_arr)// &
                        ", but found "//trim(nc_dim_h5)//"."
                end if
                allocate (tmp_real(self%register_io%dims(v), n))
                call h5f%read(arr_path, tmp_real)
                do i = 1, n
                    call self%register_io%get_variable(self%particles(i), v, ptr)
                    ptr(:) = tmp_real(:, i)
                end do
            else
                error stop "HDF5 array must be either rank 1 or 2."
            end if
            deallocate (tmp_real)
        end do
        call h5f%close()

    end subroutine base_read

    !> @brief Controls particles' summary stats printed.
    !> @param self The particles system to print the stats of.
    !> @param out_str The string which contains the summary stats and any formatting.
    subroutine base_generate_summary(self, out_str)

        class(particle_system_t), intent(in):: self
        character(len=:), allocatable, intent(out):: out_str
        integer:: offset, i
        integer, parameter:: line_length = 60, nlines = 5
        character(*), parameter:: format_str = "(4x, A, f12.5, A, I10)"
        real(fp):: minv, maxv, v
        integer:: mini, maxi

        allocate (character(nlines*line_length)::out_str)
        offset = 0

        ! save max accel
        maxv = sum(self%particles(1)%dvxdt(:)**2)
        maxi = 1
        do i = 2, self%size
            v = sum(self%particles(i)%dvxdt(:)**2)
            if (v > maxv) then
                maxv = v
                maxi = i
            end if
        end do
        write (out_str(offset + 1:offset + line_length), format_str) "  max(|dvdt|) of ", sqrt(maxv), " at particle ", maxi
        offset = offset + line_length
        write (out_str(offset:offset), "(A1)") new_line("a")

        ! save max vel
        maxv = sum(self%particles(1)%v(:)**2)
        maxi = 1
        do i = 2, self%size
            v = sum(self%particles(i)%v(:)**2)
            if (v > maxv) then
                maxv = v
                maxi = i
            end if
        end do
        write (out_str(offset + 1:offset + line_length), format_str) "     max(|v|) of ", sqrt(maxv), " at particle ", maxi
        offset = offset + line_length
        write (out_str(offset:offset), "(A1)") new_line("a")

        ! save min/max rho
        maxv = self%particles(1)%rho
        minv = self%particles(1)%rho
        maxi = 1
        mini = 1
        do i = 2, self%size
            if (self%particles(i)%rho > maxv) then
                maxv = self%particles(i)%rho
                maxi = i
            end if
            if (self%particles(i)%rho < minv) then
                minv = self%particles(i)%rho
                mini = i
            end if
        end do
        write (out_str(offset + 1:offset + line_length), format_str) "     min(rho) of ", minv, " at particle ", mini
        offset = offset + line_length
        write (out_str(offset:offset), "(A1)") new_line("a")

        ! save max rho
        write (out_str(offset + 1:offset + line_length), format_str) "     max(rho) of ", maxv, " at particle ", maxi
        offset = offset + line_length
        write (out_str(offset:offset), "(A1)") new_line("a")

        maxv = abs(self%particles(1)%drhodt)
        maxi = 1
        do i = 2, self%size
            v = abs(self%particles(i)%drhodt)
            if (v > maxv) then
                maxv = v
                maxi = i
            end if
        end do
        write (out_str(offset + 1:offset + line_length), format_str) "max(|drhodt|) of ", self%particles(maxi)%drhodt, &
            " at particle ", maxi

    end subroutine base_generate_summary

end module grasph_particle_system_m
