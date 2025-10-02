module test_time_integration

    use grasph_constants, only: fp
    use grasph_particles, only: particle_system_t
    use weakly_compressible_particles, only: eos_particle, linear_eos_state_updater
    use grasph_kernels, only: grasph_cubic_bspline_kernel
    use grasph_system_interactions_m, only: system_interaction_t
    use grasph_time_integration, only: leap_frog_time_integration
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

contains

    type(test_list) function tests()

        tests = test_list([ &
                          test("test_LF_1particle", test_LF_1particle_nointeractions) &
                          ])

    end function tests

    subroutine test_LF_1particle_nointeractions()
        type(system_interaction_t):: wcp_interaction_pairs(1)
        type(particle_system_t):: psys(1)
        type(grasph_cubic_bspline_kernel):: kernel
        type(eos_particle):: ps_template
        type(linear_eos_state_updater):: state_updater

        state_updater%rho_ref = 1000._fp
        call psys(1)%base_init(n=1, name="test", particle_template=ps_template, state_updater=state_updater)
        ! call ps%register_x%register_data(ps%x, "x", ps%v, "v")
        call psys(1)%register_v%register(psys(1)%particles(1), "v", psys(1)%particles(1)%v, psys(1)%particles(1)%dvxdt)
        call psys(1)%register_v%register(psys(1)%particles(1), "rho", psys(1)%particles(1)%rho, psys(1)%particles(1)%drhodt)
        psys(1)%particles(1)%x(1) = 1._fp
        psys(1)%particles(1)%x(2) = 2._fp
        psys(1)%particles(1)%v(1) = 3._fp
        psys(1)%particles(1)%v(2) = 4._fp
        psys(1)%particles(1)%dvxdt(1) = 10._fp
        psys(1)%particles(1)%dvxdt(2) = 20._fp
        psys(1)%particles(1)%rho = 1000._fp
        psys(1)%particles(1)%drhodt = 1000._fp
        psys(1)%particles(1)%c = 2._fp

        call wcp_interaction_pairs(1)%init(1, psys(1))

        call kernel%init(2, 1._fp)

        call leap_frog_time_integration(1, 1, 1, psys, wcp_interaction_pairs, 1._fp, kernel, "/tmp", "test-", 4)

        call check(is_close(psys(1)%particles(1)%dvxdt(1), 10._fp), "Incorrect value for dvxdt(1)") ! should be unchanged
        call check(is_close(psys(1)%particles(1)%dvxdt(2), 20._fp), "Incorrect value for dvxdt(2)") ! should be unchanged
        call check(is_close(psys(1)%particles(1)%v(1), 8._fp), "Incorrect value for v(1)") ! should be 3 + (1/2)*10
        call check(is_close(psys(1)%particles(1)%v(2), 14._fp), "Incorrect value for v(2)") ! should be 4 + (1/2)*20
        call check(is_close(psys(1)%particles(1)%x(1), 1._fp), "Incorrect value for x(1)") ! should be unchanged
        call check(is_close(psys(1)%particles(1)%x(2), 2._fp), "Incorrect value for x(2)") ! should be unchanged
        call check(is_close(psys(1)%particles(1)%drhodt, 1000._fp), "Incorrect value for drhodt") ! should be unchanged
        call check(is_close(psys(1)%particles(1)%rho, 1500._fp), "Incorrect value for rho") ! should be 1000 + (1/2)*1000
        select type (ps => psys(1)%particles)
        type is (eos_particle)
            call check(is_close(ps(1)%p, 1000._fp), "Incorrect value for p") ! should be 2**2*((1000 + 0.5*(1/2)*1000) - 1000)
        end select

        ! add the x-registration
        call psys(1)%register_x%register(psys(1)%particles(1), "x", psys(1)%particles(1)%x, psys(1)%particles(1)%v)

        call leap_frog_time_integration(1, 1, 1, psys, wcp_interaction_pairs, 1._fp, kernel, "/tmp", "test-", 4)

        call check(is_close(psys(1)%particles(1)%dvxdt(1), 10._fp), "Incorrect value for dvxdt(1)") ! should be unchanged
        call check(is_close(psys(1)%particles(1)%dvxdt(2), 20._fp), "Incorrect value for dvxdt(2)") ! should be unchanged
        call check(is_close(psys(1)%particles(1)%v(1), 13._fp), "Incorrect value for v(1)") ! should be 8 + (1/2)*10
        call check(is_close(psys(1)%particles(1)%v(2), 24._fp), "Incorrect value for v(2)") ! should be 14 + (1/2)*20
        call check(is_close(psys(1)%particles(1)%x(1), 7.5_fp), "Incorrect value for x(1)") ! should be 1 + (1/2)*(8 + (1/2)*10)
        call check(is_close(psys(1)%particles(1)%x(2), 14._fp), "Incorrect value for x(2)") ! should be 2 + (1/2)*(14 + (1/2)*20)
        call check(is_close(psys(1)%particles(1)%drhodt, 1000._fp), "Incorrect value for drhodt") ! should be unchanged
        call check(is_close(psys(1)%particles(1)%rho, 2000._fp), "Incorrect value for rho") ! should be 1500 + (1/2)*1000
        select type (ps => psys(1)%particles)
        type is (eos_particle)
            call check(is_close(ps(1)%p, 3000._fp), "Incorrect value for p") ! should be 2**2*((1500 + 0.5*(1/2)*1000) - 1000)
        end select

    end subroutine test_LF_1particle_nointeractions

end module test_time_integration

program run_tests

    use test_time_integration, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
