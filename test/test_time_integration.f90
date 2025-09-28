module test_time_integration

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles, particles_container
    use weakly_compressible_particles, only: wcp => linear_eos_particles
    use grasph_kernels, only: grasph_cubic_bspline_kernel
    use grasph_pair_sets, only: particle_interactions
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
        type(particle_interactions):: wcp_interaction_pairs(1)
        type(particles_container):: wcp_sets(1)
        type(grasph_cubic_bspline_kernel):: kernel

        allocate (wcp::wcp_sets(1)%p)

        select type (ps => wcp_setS(1)%p)
        type is (wcp)
            call ps%init(n=1, d=2, name="test", rho_ref=1000._fp)
            ! call ps%register_x%register_data(ps%x, "x", ps%v, "v")
            call ps%register_v%register_data(ps%v, "v", ps%dvxdt, "dvxdt")
            call ps%register_v%register_data(ps%rho, "rho", ps%drhodt, "drhodt")
            ps%x(1, 1) = 1._fp
            ps%x(2, 1) = 2._fp
            ps%v(1, 1) = 3._fp
            ps%v(2, 1) = 4._fp
            ps%dvxdt(1, 1) = 10._fp
            ps%dvxdt(2, 1) = 20._fp
            ps%rho(1) = 1000._fp
            ps%drhodt(1) = 1000._fp
            ps%c(1) = 2._fp

        end select

        call wcp_interaction_pairs(1)%init(1, wcp_sets(1)%p)

        call kernel%init(2, 1._fp)

        call leap_frog_time_integration(1, 1, 1, wcp_sets, wcp_interaction_pairs, 1._fp, kernel, "/tmp", "test-", 4)

        call check(is_close(wcp_sets(1)%p%dvxdt(1, 1), 10._fp), "Incorrect value for dvxdt(1)") ! should be unchanged
        call check(is_close(wcp_sets(1)%p%dvxdt(2, 1), 20._fp), "Incorrect value for dvxdt(2)") ! should be unchanged
        call check(is_close(wcp_sets(1)%p%v(1, 1), 8._fp), "Incorrect value for v(1)") ! should be 3 + (1/2)*10
        call check(is_close(wcp_sets(1)%p%v(2, 1), 14._fp), "Incorrect value for v(2)") ! should be 4 + (1/2)*20
        call check(is_close(wcp_sets(1)%p%x(1, 1), 1._fp), "Incorrect value for x(1)") ! should be unchanged
        call check(is_close(wcp_sets(1)%p%x(2, 1), 2._fp), "Incorrect value for x(2)") ! should be unchanged
        call check(is_close(wcp_sets(1)%p%drhodt(1), 1000._fp), "Incorrect value for drhodt") ! should be unchanged
        call check(is_close(wcp_sets(1)%p%rho(1), 1500._fp), "Incorrect value for rho") ! should be 1000 + (1/2)*1000
        select type (ps => wcp_sets(1)%p)
        type is (wcp)
            call check(is_close(ps%p(1), 1000._fp), "Incorrect value for p") ! should be 2**2*((1000 + 0.5*(1/2)*1000) - 1000)
        end select

        ! add the x-registration
        call wcp_sets(1)%p%register_x%register_data(wcp_sets(1)%p%x, "x", wcp_sets(1)%p%v, "v")

        call leap_frog_time_integration(1, 1, 1, wcp_sets, wcp_interaction_pairs, 1._fp, kernel, "/tmp", "test-", 4)

        call check(is_close(wcp_sets(1)%p%dvxdt(1, 1), 10._fp), "Incorrect value for dvxdt(1)") ! should be unchanged
        call check(is_close(wcp_sets(1)%p%dvxdt(2, 1), 20._fp), "Incorrect value for dvxdt(2)") ! should be unchanged
        call check(is_close(wcp_sets(1)%p%v(1, 1), 13._fp), "Incorrect value for v(1)") ! should be 8 + (1/2)*10
        call check(is_close(wcp_sets(1)%p%v(2, 1), 24._fp), "Incorrect value for v(2)") ! should be 14 + (1/2)*20
        call check(is_close(wcp_sets(1)%p%x(1, 1), 7.5_fp), "Incorrect value for x(1)") ! should be 1 + (1/2)*(8 + (1/2)*10)
        call check(is_close(wcp_sets(1)%p%x(2, 1), 14._fp), "Incorrect value for x(2)") ! should be 2 + (1/2)*(14 + (1/2)*20)
        call check(is_close(wcp_sets(1)%p%drhodt(1), 1000._fp), "Incorrect value for drhodt") ! should be unchanged
        call check(is_close(wcp_sets(1)%p%rho(1), 2000._fp), "Incorrect value for rho") ! should be 1500 + (1/2)*1000
        select type (ps => wcp_sets(1)%p)
        type is (wcp)
            call check(is_close(ps%p(1), 3000._fp), "Incorrect value for p") ! should be 2**2*((1500 + 0.5*(1/2)*1000) - 1000)
        end select

    end subroutine test_LF_1particle_nointeractions

end module test_time_integration

program run_tests

    use test_time_integration, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
