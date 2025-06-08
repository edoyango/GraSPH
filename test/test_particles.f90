module test_particles

    use grasph_constants, only: fp
    use grasph_kernels, only: grasph_cubic_bspline_kernel
    use grasph_pairs, only: particle_pairs
    use grasph_particles, only: wc_particles => weakly_compressible_particles, base_particles
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

    ! data for testing
    real(fp), parameter:: dx = 0.25_fp

contains

    type(test_list) function tests()

        tests = test_list([ &
            test("test_particles_init", test_particles_init), &
            test("test_linear_eos_wc_particles", test_linear_eos_wc_particles) &
        ])

    end function tests

    subroutine test_particles_init()

        type(wc_particles):: ps

        call ps%init(16, 2, "test", 0._fp)

        ! check name assigned correctly
        call check(ps%name == "test", "Particle set name not initialized to 'test'")

        ! check member values set correctly
        call check(ps%initialized, "Particle initilization logical not set to .true.")
        call check(is_equal(ps%ndims, 2), "Particle ndims not set correctly")
        call check(is_equal(ps%size, 16), "Particle size not set correctly")

        ! check arrays are allocated and sized correctly
        call check(allocated(ps%x), "Particle x not allocated")
        call check(is_equal(size(ps%x, 1), 2), "Particle x dim 1 incorrect")
        call check(is_equal(size(ps%x, 2), 16), "Particle x dim 2 incorrect")
        call check(allocated(ps%v), "Particle v not allocated")
        call check(is_equal(size(ps%v, 1), 2), "Particle v dim 1 incorrect")
        call check(is_equal(size(ps%v, 2), 16), "Particle v dim 2 incorrect")
        call check(allocated(ps%rho), "Particle rho not allocated")
        call check(is_equal(size(ps%rho), 16), "Particle rho size incorrect")
        call check(allocated(ps%mass), "Particle mass not allocated")
        call check(is_equal(size(ps%mass), 16), "Particle mass size incorrect")

        ! check extended type array(s)
        call check(allocated(ps%p), "Particle p not allocated")
        call check(is_equal(size(ps%p), 16), "Particle p size incorrect")

        ! 3d
        call ps%init(27, 3, "test", 0._fp)

        ! check member values set correctly
        call check(ps%initialized, "Particle initilization logical not set to .true.")
        call check(is_equal(ps%ndims, 3), "Particle ndims not set correctly")
        call check(is_equal(ps%size, 27), "Particle size not set correctly")

        ! check arrays are allocated and sized correctly
        call check(allocated(ps%x), "Particle x not allocated")
        call check(is_equal(size(ps%x, 1), 3), "Particle x dim 1 incorrect")
        call check(is_equal(size(ps%x, 2), 27), "Particle x dim 2 incorrect")
        call check(allocated(ps%v), "Particle v not allocated")
        call check(is_equal(size(ps%v, 1), 3), "Particle v dim 1 incorrect")
        call check(is_equal(size(ps%v, 2), 27), "Particle v dim 2 incorrect")
        call check(allocated(ps%rho), "Particle rho not allocated")
        call check(is_equal(size(ps%rho), 27), "Particle rho size incorrect")
        call check(allocated(ps%mass), "Particle mass not allocated")
        call check(is_equal(size(ps%mass), 27), "Particle mass size incorrect")

        ! check extended type array(s)
        call check(allocated(ps%p), "Particle p not allocated")
        call check(is_equal(size(ps%p), 27), "Particle p size incorrect")
        
    end subroutine test_particles_init

    subroutine state_update_test1(self, dt)
        class(base_particles), intent(inout):: self
        real(fp), intent(in):: dt
        integer:: i
        do i = 1, self%size
            self%rho(i) = self%rho(i)*2._fp
        enddo
    end subroutine state_update_test1

    subroutine state_update_test2(self, dt)
        class(wc_particles), intent(inout):: self
        real(fp), intent(in):: dt
        integer:: i
        do i = 1, self%size
            self%rho(i) = self%rho(i)*2._fp
        enddo
    end subroutine state_update_test2

    subroutine test_linear_eos_wc_particles()

        type(wc_particles):: ps1
        integer:: i
        character:: ic
        
        call ps1%init(5, 2, "test", 1._fp)

        do i= 1, 5
            ps1%rho(i) = real(i, kind=fp)
            ps1%c(i) = 2._fp
        enddo

        call ps1%state_update()

        do i = 1, 5
            write(ic, "(I1)") i
            call check( &
                is_close(ps1%p(i), 4._fp*real(i-1, kind=fp)), &
                "State update function not applied correctly to particle " // ic &
            )
        enddo

    end subroutine test_linear_eos_wc_particles

end module test_particles

program run_tests

    use test_particles, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests