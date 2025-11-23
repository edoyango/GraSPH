module test_time_integration

    use grasph_constants_m, only: fp
    use grasph_particle_system_m, only: particle_system_t, state_updater_container_t
    use weakly_compressible_particles_m, only: eos_particles_t, linear_eos_state_updater_t
    use grasph_kernels_m, only: cubic_bspline_kernel_t
    use grasph_system_interactions_m, only: system_interaction_t, sweeper_container_t, default_sweeper_t
    use grasph_time_integration_m, only: leap_frog_time_integration
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
        type(cubic_bspline_kernel_t):: kernel
        type(eos_particles_t):: ps_template
        type(linear_eos_state_updater_t), pointer:: state_updater
        type(sweeper_container_t):: sweepers(1)
        type(state_updater_container_t):: state_updaters(1)

        allocate (default_sweeper_t::sweepers(1)%sweeper)
        allocate (linear_eos_state_updater_t::state_updaters(1)%updater)
        select type (updater => state_updaters(1)%updater)
        type is (linear_eos_state_updater_t)
            state_updater => updater
        end select
        state_updater%rho_ref = 1000._fp
        call psys(1)%init(n=1, name="test", particle_template=ps_template, state_updaters=state_updaters)
        ! call ps%register_x%register_data(ps%x, "x", ps%v, "v")
        call psys(1)%register_v%register("v", psys(1)%particles%v, psys(1)%particles%dvxdt)
        call psys(1)%register_v%register("rho", psys(1)%particles%rho, psys(1)%particles%drhodt)
        psys(1)%particles%x(1, 1) = 1._fp
        psys(1)%particles%x(2, 1) = 2._fp
        psys(1)%particles%v(1, 1) = 3._fp
        psys(1)%particles%v(2, 1) = 4._fp
        psys(1)%particles%dvxdt(1, 1) = 10._fp
        psys(1)%particles%dvxdt(2, 1) = 20._fp
        psys(1)%particles%rho(1) = 1000._fp
        psys(1)%particles%drhodt(1) = 1000._fp
        psys(1)%particles%c(1) = 2._fp

        call wcp_interaction_pairs(1)%init(1, psys(1), sweepers=sweepers)

        call kernel%init(2, 1._fp)

        call leap_frog_time_integration(1, 1, 1, psys, wcp_interaction_pairs, 1._fp, kernel, "/tmp", "test-", 4)

        call check(is_close(psys(1)%particles%dvxdt(1, 1), 10._fp), "Incorrect value for dvxdt(1)") ! should be unchanged
        call check(is_close(psys(1)%particles%dvxdt(2, 1), 20._fp), "Incorrect value for dvxdt(2)") ! should be unchanged
        call check(is_close(psys(1)%particles%v(1, 1), 8._fp), "Incorrect value for v(1)") ! should be 3 + (1/2)*10
        call check(is_close(psys(1)%particles%v(2, 1), 14._fp), "Incorrect value for v(2)") ! should be 4 + (1/2)*20
        call check(is_close(psys(1)%particles%x(1, 1), 1._fp), "Incorrect value for x(1)") ! should be unchanged
        call check(is_close(psys(1)%particles%x(2, 1), 2._fp), "Incorrect value for x(2)") ! should be unchanged
        call check(is_close(psys(1)%particles%drhodt(1), 1000._fp), "Incorrect value for drhodt") ! should be unchanged
        call check(is_close(psys(1)%particles%rho(1), 1500._fp), "Incorrect value for rho") ! should be 1000 + (1/2)*1000
        select type (ps => psys(1)%particles)
        type is (eos_particles_t)
            call check(is_close(ps%p(1), 1000._fp), "Incorrect value for p") ! should be 2**2*((1000 + 0.5*(1/2)*1000) - 1000)
        end select

        ! add the x-registration
        call psys(1)%register_x%register("x", psys(1)%particles%x, psys(1)%particles%v)

        call leap_frog_time_integration(1, 1, 1, psys, wcp_interaction_pairs, 1._fp, kernel, "/tmp", "test-", 4)

        call check(is_close(psys(1)%particles%dvxdt(1, 1), 10._fp), "Incorrect value for dvxdt(1)") ! should be unchanged
        call check(is_close(psys(1)%particles%dvxdt(2, 1), 20._fp), "Incorrect value for dvxdt(2)") ! should be unchanged
        call check(is_close(psys(1)%particles%v(1, 1), 13._fp), "Incorrect value for v(1)") ! should be 8 + (1/2)*10
        call check(is_close(psys(1)%particles%v(2, 1), 24._fp), "Incorrect value for v(2)") ! should be 14 + (1/2)*20
        call check(is_close(psys(1)%particles%x(1, 1), 7.5_fp), "Incorrect value for x(1)") ! should be 1 + (1/2)*(8 + (1/2)*10)
        call check(is_close(psys(1)%particles%x(2, 1), 14._fp), "Incorrect value for x(2)") ! should be 2 + (1/2)*(14 + (1/2)*20)
        call check(is_close(psys(1)%particles%drhodt(1), 1000._fp), "Incorrect value for drhodt") ! should be unchanged
        call check(is_close(psys(1)%particles%rho(1), 2000._fp), "Incorrect value for rho") ! should be 1500 + (1/2)*1000
        select type (ps => psys(1)%particles)
        type is (eos_particles_t)
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
