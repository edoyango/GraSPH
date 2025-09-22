module test_time_integration

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles, particles_container
    use weakly_compressible_particles, only: wcp => linear_eos_particles
    use grasph_kernels, only: grasph_base_kernel, grasph_cubic_bspline_kernel
    use grasph_pairs, only: cell_list_search
    use grasph_pair_sets, only: particle_interactions_base, particle_interactions_container
    use grasph_time_integration, only: leap_frog_time_integration
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

    type, extends(particle_interactions_base):: wcp_interactions
    end type wcp_interactions

    type, extends(particle_interactions_base):: wcp_virt_interactions
    end type wcp_virt_interactions

contains

    type(test_list) function tests()

        tests = test_list([ &
                          test("test_compile", test_compile) &
                          ])

    end function tests

    subroutine test_compile()
        type(particle_interactions_container):: wcp_interaction_pairs(2)
        type(particles_container):: wcp_sets(2)
        type(grasph_cubic_bspline_kernel):: kernel

        allocate (wcp::wcp_sets(1)%p)
        allocate (wcp::wcp_sets(2)%p)

        allocate (wcp_interactions::wcp_interaction_pairs(1)%pi)
        allocate (wcp_virt_interactions::wcp_interaction_pairs(2)%pi)

        select type (ps => wcp_setS(1)%p)
        type is (wcp)
            call ps%init(n=1, d=2, name="test", rho_ref=1000._fp)
        end select
        select type (ps => wcp_sets(2)%p)
        type is (wcp)
            call ps%init(n=1, d=2, name="test", rho_ref=1000._fp)
        end select

        call wcp_interaction_pairs(1)%pi%base_init(1, wcp_sets(1)%p)
        call wcp_interaction_pairs(2)%pi%base_init(1, wcp_sets(1)%p, wcp_sets(2)%p)

        call kernel%init(2, 1._fp)

        call leap_frog_time_integration(2, 1, 1, wcp_sets, wcp_interaction_pairs, 1._fp, kernel, "/tmp", "test-", 4)

    end subroutine test_compile

end module test_time_integration

program run_tests

    use test_time_integration, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
