module test_particles

    use grasph_constants, only: fp, ndims
    use grasph_kernels, only: grasph_cubic_bspline_kernel
    use grasph_pairs, only: particle_pairs
    use grasph_particles, only: base_particles
    use weakly_compressible_particles, only: wc_particles => linear_eos_particles, linear_eos_particle, linear_eos_state_updater
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
        type(linear_eos_particle):: ps_template

        call ps%init(n=16, name="test", ps_template=ps_template)

        ! check name assigned correctly
        call check(ps%name == "test", "Particle set name not initialized to 'test'")

        ! check member values set correctly
        call check(ps%initialized, "Particle initilization logical not set to .true.")
        call check(is_equal(ps%ndims, ndims), "Particle ndims not set correctly")
        call check(is_equal(ps%size, 16), "Particle size not set correctly")

        ! check arrays are allocated and sized correctly
        call check(allocated(ps%ps), "Particle array not allocated")
        call check(is_equal(size(ps%ps), 16), "Particle array size incorrect")

    end subroutine test_particles_init

    subroutine test_linear_eos_wc_particles()

        type(wc_particles):: ps1
        type(linear_eos_particle):: ps_template
        type(linear_eos_state_updater):: state_updater
        integer:: i
        character:: ic

        state_updater%rho_ref = 1._fp
        call ps1%init(n=5, name="test", ps_template=ps_template, state_updater=state_updater)

        do i = 1, 5
            ps1%ps(i)%rho = real(i, kind=fp)
            ps1%ps(i)%c = 2._fp
        end do

        call ps1%do_state_update()

        do i = 1, 5
            write (ic, "(I1)") i
            select type (ps => ps1%ps)
            class is (linear_eos_particle)
                call check( &
                    is_close(ps(i)%p, 4._fp*real(i - 1, kind=fp)), &
                    "State update function not applied correctly to particle "//ic &
                    )
            end select
        end do

    end subroutine test_linear_eos_wc_particles

end module test_particles

program run_tests

    use test_particles, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
