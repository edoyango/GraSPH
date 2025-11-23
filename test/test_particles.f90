module test_particles

    use grasph_constants_m, only: fp, ndims
    use grasph_particle_m, only: base_particle_t
    use grasph_particle_system_m, only: particle_system_t, state_updater_container_t
    use weakly_compressible_particles_m, only: linear_eos_state_updater_t
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
                          test("test_linear_eos_wc_particles", test_linear_eos_wc_particles), &
                          test("test_update_del_methods", test_update_del_methods) &
                          ])

    end function tests

    subroutine test_particles_init()

        type(particle_system_t):: psys

        call psys%init(n=16, name="test")

        ! check name assigned correctly
        call check(psys%name == "test", "Particle set name not initialized to 'test'")

        ! check member values set correctly
        call check(psys%initialized, "Particle initilization logical not set to .true.")
        call check(is_equal(psys%ndims, ndims), "Particle ndims not set correctly")
        call check(is_equal(psys%size, 16), "Particle size not set correctly")

        ! check arrays are allocated and sized correctly
        call check(allocated(psys%particles), "Particle array not allocated")
        call check(is_equal(size(psys%particles), 16), "Particle array size incorrect")

    end subroutine test_particles_init

    subroutine test_linear_eos_wc_particles()

        type(particle_system_t):: psys
        type(state_updater_container_t):: updaters(1)
        integer:: i
        character:: ic

        allocate (linear_eos_state_updater_t::updaters(1)%updater)
        select type (updater => updaters(1)%updater)
        type is (linear_eos_state_updater_t)
            updater%rho_ref = 1._fp
        end select
        call psys%init(n=5, name="test", state_updaters=updaters)

        do i = 1, 5
            psys%particles(i)%rho = real(i, kind=fp)
            psys%particles(i)%c = 2._fp
        end do

        call psys%do_state_update(1)

        do i = 1, 5
            write (ic, "(I1)") i
            call check( &
                is_close(psys%particles(i)%p, 4._fp*real(i - 1, kind=fp)), &
                "State update function not applied correctly to particle "//ic &
                )
        end do

    end subroutine test_linear_eos_wc_particles

    subroutine test_update_del_methods()

        type(particle_system_t):: psys
        integer:: i, d
        character:: dc

        call psys%init(1, "test-particles")

        call check(is_equal(psys%safe_size_plus_1(), 2), "safe_size_plus_1 didn't return size + 1.")
        call check(is_equal(psys%size, 2), "safe_size_plus_1 didn't update size of psys.")
        call check(is_equal(size(psys%particles), 2), "safe_size_plus_1 didn't allocate correct space for psys%particles.")

        call check(is_equal(psys%safe_size_plus_1(), 3), "safe_size_plus_1 didn't return size + 1.")
        call check(is_equal(psys%size, 3), "safe_size_plus_1 didn't update size of psys.")
        call check(is_equal(size(psys%particles), 4), "safe_size_plus_1 didn't allocate correct space for psys%particles.")

    end subroutine test_update_del_methods

end module test_particles

program run_tests

    use test_particles, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
