!> @file weakly_compressible_particles.f90
!> @brief Module containing weakly compressible particles type and methods
!> @author Edward Yang
!> @date 2025-09-21
module weakly_compressible_particles

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles, base_dump, base_read

    !> @brief particle type which adds pressure, determined from density with a linear EOS
    type, extends(base_particles):: linear_eos_particles
        !> @brief Particles' isotropic pressure.
        real(fp), allocatable:: p(:)
        !> @brief Reference density to be used to calculate pressure in the linear EOS.
        real(fp):: rho_ref
    contains
        !> @brief Custom intializer to initialize pressure and reference density.
        !>        Also calls base_init to initialize base data.
        procedure:: init => wcp_init
        !> @brief Linear equation of state which overrides the do-nothing base state-update subroutine.
        procedure:: state_update => linear_eos
        !> @brief Overrides base output dump to include pressure data.
        procedure:: dump => wcp_dump
        !> @brief Overrides base input read to include pressure data.
        procedure:: read => wcp_read
    end type linear_eos_particles

    !> @brief particle type which adds pressure, determined from density with a Tait EOS
    type, extends(linear_eos_particles):: tait_eos_particles
    contains
        !> @brief Tait equation of state which overrides the do-nothing base state-update subroutine.
        procedure:: state_update => tait_eos
    end type tait_eos_particles

contains

    !> @brief Custom init function for weakly-compressible particles. Will also initialize base
    !>        particles' data.
    !> @param self The weakly-compressible particles to initialize.
    !> @param n Number of particles to allocate space for.
    !> @param d Spatial dimensions of the particles.
    !> @param name A label to give the particles. Used to label output/terminal information.
    !> @param rho_ref Reference density used in the linear EOS.
    subroutine wcp_init(self, n, d, name, rho_ref)
        class(linear_eos_particles), intent(inout):: self
        integer, intent(in):: n, d
        real(fp), intent(in):: rho_ref
        character(*), intent(in):: name
        self%rho_ref = rho_ref
        if (self%initialized) deallocate (self%p)
        call self%base_init(n, d, name)
        allocate (self%p(n))
    end subroutine wcp_init

    !> @brief Custom output subroutine to include relevant weakly-compressible data. Overrides
    !>        base particles' dump method.
    !> @param self The weakly-compressible particles to write.
    !> @param itimestep The timestep to add to the filename.
    !> @param path The output directory.
    !> @param prefix_in The prefix to give to the output filenames.
    !> @param comp_level The level of gzip compression to use.
    subroutine wcp_dump(self, itimestep, path, prefix_in, comp_level)
        use h5fortran, only: hdf5_file
        class(linear_eos_particles), intent(in):: self
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
        end if

        write (ic, "(I10.10)") itimestep
        file_path = path//"/"//trim(filename_prefix)//"_"//ic//".h5"
        this_group = "/"//trim(self%name)//"/"//group

        call base_dump(self, itimestep, path, prefix_in, comp_level)
        call h5f%open(file_path, action="a")
        call h5f%write(trim(this_group)//"p", self%p)
        call h5f%close()

    end subroutine wcp_dump

    !> @brief Reads weakly-compressible particle data from HDF5 file.
    !> @param self The weakly-compressible particles to read data into.
    !> @param file_path The path to the file to read.
    !> @param name The name to of particles to read and assign to the read particles.
    subroutine wcp_read(self, file_path, name)
        use h5fortran, only: hdf5_file
        class(linear_eos_particles), intent(out):: self
        character(*), intent(in):: file_path, name
        character(*), parameter:: group = "weakly_compressible/"
        character(200):: filename, this_group
        integer:: ierr, d, n
        type(hdf5_file):: h5f
        character(10):: ic
        call base_read(self, file_path, name)

        this_group = "/"//trim(name)//"/"//group

        allocate (self%p(self%size))

        call h5f%open(file_path, action="r")
        call h5f%read(trim(this_group)//"p", self%p)
        call h5f%close()
    end subroutine wcp_read

    !> @brief The linear state equation to update stress using the particles' speed of sound (c),
    !>        density (rho), and reference density (rho_ref). Overrides base_particles' state_update
    !>        subroutine.
    !> @param self The particles' pressure to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine linear_eos(self, dt)
        class(linear_eos_particles), intent(inout):: self
        real(fp), intent(in), optional:: dt
        integer:: i
        do i = 1, self%size
            self%p(i) = self%c(i)**2*(self%rho(i) - self%rho_ref)
        end do
    end subroutine linear_eos

    !> @brief The Tait state equation to update stress using the particles' speed of sound (c),
    !>        density (rho), and reference density (rho_ref). Overrides base_particles' state_update
    !>        subroutine.
    !> @param self The particles' pressure to be updated.
    !> @param dt The input time-increment (unused - included to match the overriden method).
    subroutine tait_eos(self, dt)
        class(tait_eos_particles), intent(inout):: self
        real(fp), intent(in), optional:: dt
        integer:: i
        integer, parameter:: gamma = 7
        do i = 1, self%size
            self%p(i) = self%rho_ref*self%c(i)*self%c(i)/real(gamma, kind=fp)*((self%rho(i)/self%rho_ref)**gamma - 1._fp)
        end do
    end subroutine tait_eos

end module weakly_compressible_particles
