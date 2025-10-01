!> @file grasph_particles.f90
!> @brief Module containing core particles derived types and methods
!> @author Edward Yang
!> @date 2025-06-09
module grasph_particles

    use iso_fortran_env, only: error_unit
    use grasph_constants, only: fp, ndims
    use grasph_common, only: array_pointer_container
    use iso_c_binding, only: c_intptr_t, c_f_pointer, c_ptr, c_loc

    implicit none
    private

    type:: base_particle
        integer:: id
        integer:: type
        real(fp):: x(ndims)
        real(fp):: v(ndims)
        real(fp):: rho
        real(fp):: mass
        real(fp):: c
        real(fp):: dvxdt(ndims)
        real(fp):: drhodt
    end type base_particle

    !> @brief Max allowable registrations. Will be removed for a more dynamic approach.
    integer, parameter:: max_registrations = 20

    type variable_register
        integer:: nregistrations = 0
        integer:: dims(max_registrations)
        character(20):: names(max_registrations)
        integer(c_intptr_t):: offsets(max_registrations)
    contains
        procedure:: register_variable_vector, register_variable_scalar
        generic:: register_variable => register_variable_vector, register_variable_scalar
        procedure:: get_variable
    end type variable_register

    type, extends(variable_register):: variable_deriv_register
        integer(c_intptr_t):: deriv_offsets(max_registrations)
    contains
        procedure:: register_variable_deriv_vector, register_variable_deriv_scalar
        generic:: register => register_variable_deriv_vector, register_variable_deriv_scalar
        procedure:: get
    end type variable_deriv_register

    !> @brief The core particles derived type
    type:: base_particles
        class(base_particle), allocatable:: ps(:)
        !> @brief Whether the particles have been initialized.
        logical:: initialized = .false.
        !> @brief Whether to print information when generate_summary is called.
        logical:: to_print_summary = .true.
        !> @brief Number of spatial dimensions
        integer:: ndims = ndims
        !> @brief Number of particles.
        integer:: size = 0
        !> @brief Name used in naming groups in output hdf5 file
        character(100):: name
        !> @brief Register for variables to be updated only at full-timestep e.g. position (x).
        type(variable_deriv_register):: register_x
        !> @brief Register for variables to be updated at mid- and full-timestep e.g. velocity (v) and density (rho).
        type(variable_deriv_register):: register_v
        !> @brief Register for variables to be written/read.
        type(variable_register):: register_io
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

    public:: base_particle, base_particles, particles_container, max_registrations

contains

    subroutine register_variable_deriv_vector(self, base, name, member, member_deriv)
        class(variable_deriv_register), intent(inout):: self
        class(base_particle), target, intent(in):: base
        character(*), intent(in):: name
        real(fp), target, intent(in):: member(:), member_deriv(:)

        if (size(member) == 0) error stop "Cannot register 0-size member variable."
        if (size(member_deriv) == 0) error stop "Cannot register 0-size member derivative variable."
        if (size(member) /= size(member_deriv)) error stop "member and member_deriv are not same size."

        call register_variable_deriv_scalar(self, base, name, member(1), member_deriv(1))

        self%dims(self%nregistrations) = size(member)

    end subroutine register_variable_deriv_vector

    subroutine register_variable_deriv_scalar(self, base, name, member, member_deriv)
        class(variable_deriv_register), intent(inout):: self
        class(base_particle), target, intent(in):: base
        character(*), intent(in):: name
        real(fp), target, intent(in):: member, member_deriv
        integer(c_intptr_t):: base_addr

        if (self%nregistrations == max_registrations) error stop "Exceeded maximum variable registrations."

        base_addr = transfer(c_loc(base%id), base_addr)

        self%nregistrations = self%nregistrations + 1
        self%names(self%nregistrations) = name
        self%offsets(self%nregistrations) = get_offset_(base_addr, member)
        self%deriv_offsets(self%nregistrations) = get_offset_(base_addr, member_deriv)
        self%dims(self%nregistrations) = 1

        if (self%offsets(self%nregistrations) >= sizeof(base)) error stop "Member is not a subset of base."
        if (self%deriv_offsets(self%nregistrations) >= sizeof(base)) error stop "member_deriv is not a subset of base."

    end subroutine register_variable_deriv_scalar

    elemental integer(c_intptr_t) function get_offset_(base_addr, member)
        integer(c_intptr_t), intent(in):: base_addr
        real(fp), target, intent(in):: member

        get_offset_ = transfer(c_loc(member), base_addr) - base_addr

    end function get_offset_

    subroutine get(self, base, idx, ptr, ptr_deriv)
        class(variable_deriv_register), intent(in):: self
        class(base_particle), target, intent(in):: base
        integer, intent(in):: idx
        real(fp), pointer, intent(out):: ptr(:), ptr_deriv(:)
        integer(c_intptr_t):: base_addr

        base_addr = transfer(c_loc(base%id), base_addr)

        call ptr_from_offset_(base_addr, self%offsets(idx), self%dims(idx), ptr)
        call ptr_from_offset_(base_addr, self%deriv_offsets(idx), self%dims(idx), ptr_deriv)

    end subroutine get

    subroutine ptr_from_offset_(base_addr, offset, dims, ptr)
        integer(c_intptr_t), intent(in):: base_addr, offset
        integer, intent(in):: dims
        real(fp), pointer, intent(out):: ptr(:)
        integer(c_intptr_t):: member_offset
        type(c_ptr):: member_c_ptr

        ! calculate address of member
        member_offset = base_addr + offset
        ! convert address to c_ptr
        member_c_ptr = transfer(member_offset, member_c_ptr)
        ! convert c_ptr to Fortran pointer
        call c_f_pointer(member_c_ptr, ptr, [dims])

    end subroutine ptr_from_offset_

    subroutine register_variable_vector(self, base, name, member)
        class(variable_register), intent(inout):: self
        class(base_particle), target, intent(in):: base
        character(*), intent(in):: name
        real(fp), target, intent(in):: member(:)

        if (size(member) == 0) error stop "Cannot register 0-size member variable."

        call register_variable_scalar(self, base, name, member(1))

        self%dims(self%nregistrations) = size(member)

    end subroutine register_variable_vector

    subroutine register_variable_scalar(self, base, name, member)
        class(variable_register), intent(inout):: self
        class(base_particle), target, intent(in):: base
        character(*), intent(in):: name
        real(fp), target, intent(in):: member
        integer(c_intptr_t):: base_addr

        if (self%nregistrations == max_registrations) error stop "Exceeded maximum variable registrations."

        base_addr = transfer(c_loc(base%id), base_addr)

        self%nregistrations = self%nregistrations + 1
        self%names(self%nregistrations) = name
        self%offsets(self%nregistrations) = get_offset_(base_addr, member)
        self%dims(self%nregistrations) = 1

        if (self%offsets(self%nregistrations) >= sizeof(base)) error stop "Member is not a subset of base."

    end subroutine register_variable_scalar

    subroutine get_variable(self, base, idx, ptr)
        class(variable_register), intent(in):: self
        class(base_particle), target, intent(in):: base
        integer, intent(in):: idx
        real(fp), pointer, intent(out):: ptr(:)
        integer(c_intptr_t):: base_addr

        base_addr = transfer(c_loc(base%id), base_addr)

        call ptr_from_offset_(base_addr, self%offsets(idx), self%dims(idx), ptr)

    end subroutine get_variable

    !> @brief Initializes base_particles' internal arrays.
    !> @param self The particles to initialize.
    !> @param n The number of particles.
    !> @param name A label to give the particles. Used to label output/terminal information.
    subroutine base_init(self, n, name, ps_template)
        class(base_particles), intent(inout):: self
        integer, intent(in):: n
        character(*), intent(in):: name
        class(base_particle), optional, intent(in):: ps_template
        if (self%initialized) call self%base_clear()
        if (present(ps_template)) then
            allocate (self%ps(n), source=ps_template)
        else
            allocate (self%ps(n))
        end if

        self%initialized = .true.
        self%size = n
        self%ndims = ndims
        self%name = name

        call self%register_io%register_variable(self%ps(1), "x", self%ps(1)%x)
        call self%register_io%register_variable(self%ps(1), "v", self%ps(1)%v)
        call self%register_io%register_variable(self%ps(1), "rho", self%ps(1)%rho)
        call self%register_io%register_variable(self%ps(1), "mass", self%ps(1)%mass)
        call self%register_io%register_variable(self%ps(1), "c", self%ps(1)%c)
        call self%register_io%register_variable(self%ps(1), "dvxdt", self%ps(1)%dvxdt)
        call self%register_io%register_variable(self%ps(1), "drhodt", self%ps(1)%drhodt)
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
        integer:: ierr, i, v
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
        this_group = "/"//trim(self%name)//"/"//group

        call h5f%open(file_path, action="a", comp_lvl=comp_level)
        call h5f%write("/"//trim(self%name)//"/n", self%size)
        call h5f%write("/"//trim(self%name)//"/ndims", self%ndims)
        allocate (tmp_int(self%size))
        do i = 1, self%size
            tmp_int(i) = self%ps(i)%id
        end do
        call h5f%write(trim(this_group)//"id", tmp_int)
        do i = 1, self%size
            tmp_int(i) = self%ps(i)%type
        end do
        call h5f%write(trim(this_group)//"type", tmp_int)
        do v = 1, self%register_io%nregistrations
            if (allocated(tmp_real)) deallocate (tmp_real)
            if (self%register_io%dims(v) == 1) then
                allocate (tmp_real(self%size, 1))
                do i = 1, self%size
                    call self%register_io%get_variable(self%ps(i), v, ptr)
                    tmp_real(i, 1) = ptr(1)
                end do
                call h5f%write(trim(this_group)//trim(self%register_io%names(v)), tmp_real(:, 1))
            else
                allocate (tmp_real(self%register_io%dims(v), self%size))
                do i = 1, self%size
                    call self%register_io%get_variable(self%ps(i), v, ptr)
                    tmp_real(:, i) = ptr(:)
                end do
                call h5f%write(trim(this_group)//trim(self%register_io%names(v)), tmp_real)
            end if
        end do
        call h5f%close()

    end subroutine base_dump

    !> @brief Reads base particle data from HDF5 file.
    !> @param self The particles to read data into.
    !> @param file_path The path to the file to read.
    !> @param name The name to of particles to read and assign to the read particles.
    subroutine base_read(self, file_path, name, ps_template)
        use h5fortran, only: hdf5_file, hsize_t
        class(base_particles), intent(inout):: self
        character(*), intent(in):: name, file_path
        class(base_particle), optional, intent(in):: ps_template
        character(*), parameter:: group = "base/"
        character(200):: this_group
        character(250):: arr_path
        integer:: d, n, i, v, nrank
        type(hdf5_file):: h5f
        integer, allocatable:: tmp_int(:)
        integer(hsize_t), allocatable:: dims(:)
        real(fp), allocatable:: tmp_real(:, :)
        real(fp), pointer:: ptr(:)
        character(2):: nc_dim_arr, nc_dim_h5

        this_group = "/"//trim(name)//"/"//group

        call h5f%open(file_path, action="r")
        call h5f%read("/"//trim(name)//"/n", n)
        call h5f%read("/"//trim(name)//"/ndims", d)
        if (d /= ndims) error stop "Input HDF5 file dimensions don't match code dimensions."
        allocate (tmp_int(n))
        call h5f%read(trim(this_group)//"id", tmp_int)
        do i = 1, n
            self%ps(i)%id = tmp_int(i)
        end do
        call h5f%read(trim(this_group)//"type", tmp_int)
        do i = 1, n
            self%ps(i)%type = tmp_int(i)
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
                    call self%register_io%get_variable(self%ps(i), v, ptr)
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
                    call self%register_io%get_variable(self%ps(i), v, ptr)
                    ptr(:) = tmp_real(:, i)
                end do
            else
                error stop "HDF5 array must be either rank 1 or 2."
            end if
            deallocate (tmp_real)
        end do
        call h5f%close()

    end subroutine base_read

    !> @brief Controls particles' summary stats.
    !> @param self The particles to print the stats of.
    !> @param out_str The string which contains the summary stats and any formatting.
    subroutine base_generate_summary(self, out_str)

        class(base_particles), intent(in):: self
        character(len=:), allocatable, intent(out):: out_str
        integer:: offset, d, i
        integer, parameter:: line_length = 60, nlines = 5
        character(*), parameter:: format_str = "(4x, A, f12.5, A, I10)"
        real(fp):: minv, maxv, v
        integer:: mini, maxi

        allocate (character(nlines*line_length)::out_str)
        offset = 0

        ! save max accel
        maxv = sum(self%ps(1)%dvxdt(:)**2)
        maxi = 1
        do i = 2, self%size
            v = sum(self%ps(i)%dvxdt(:)**2)
            if (v > maxv) then
                maxv = v
                maxi = i
            end if
        end do
        write (out_str(offset + 1:offset + line_length), format_str) "  max(|dvdt|) of ", sqrt(maxv), " at particle ", maxi
        offset = offset + line_length
        write (out_str(offset:offset), "(A1)") new_line("a")

        ! save max vel
        maxv = sum(self%ps(1)%v(:)**2)
        maxi = 1
        do i = 2, self%size
            v = sum(self%ps(i)%v(:)**2)
            if (v > maxv) then
                maxv = v
                maxi = i
            end if
        end do
        write (out_str(offset + 1:offset + line_length), format_str) "     max(|v|) of ", sqrt(maxv), " at particle ", maxi
        offset = offset + line_length
        write (out_str(offset:offset), "(A1)") new_line("a")

        ! save min/max rho
        maxv = self%ps(1)%rho
        minv = self%ps(1)%rho
        maxi = 1
        mini = 1
        do i = 2, self%size
            if (self%ps(i)%rho > maxv) then
                maxv = self%ps(i)%rho
                maxi = i
            end if
            if (self%ps(i)%rho < minv) then
                minv = self%ps(i)%rho
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

        maxv = abs(self%ps(1)%drhodt)
        maxi = 1
        do i = 2, self%size
            v = abs(self%ps(i)%drhodt)
            if (v > maxv) then
                maxv = v
                maxi = i
            end if
        end do
        write (out_str(offset + 1:offset + line_length), format_str) "max(|drhodt|) of ", self%ps(maxi)%drhodt, " at particle ", &
            maxi

    end subroutine base_generate_summary

end module grasph_particles
