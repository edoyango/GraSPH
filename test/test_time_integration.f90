module test_time_integration

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles, wcp => weakly_compressible_particles, particles_container
    use grasph_kernels, only: grasph_base_kernel, grasph_cubic_bspline_kernel
    use grasph_pairs, only: cell_list_search
    use grasph_pair_sets, only: particle_interactions_base, particle_interactions_container
    use grasph_time_integration, only: leap_frog_time_integration
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

    type, extends(particle_interactions_base):: wcp_interactions
        class(wcp), pointer:: ps
    contains
        procedure:: init => wcp_interactions_init
    end type wcp_interactions

    type, extends(particle_interactions_base):: wcp_virt_interactions
        class(wcp), pointer:: ps_real, ps_virt
    contains
        procedure:: init => wcp_virt_interactions_init
    end type wcp_virt_interactions

contains

    type(test_list) function tests()

        tests = test_list([ &
            test("test_compile", test_compile) &
        ])

    end function tests

    subroutine wcp_interactions_init(self, npairs_per_particle, ps_lhs, ps_rhs)
        class(wcp_interactions), intent(out):: self
        class(base_particles), target, intent(in):: ps_lhs
        class(base_particles), target, intent(in), optional:: ps_rhs
        integer, intent(in):: npairs_per_particle
        select type (ps => ps_lhs)
        type is (wcp)
            self%ps => ps
        end select
        call self%base_init(npairs_per_particle, ps_lhs)
    end subroutine wcp_interactions_init

    subroutine wcp_virt_interactions_init(self, npairs_per_particle, ps_lhs, ps_rhs)
        class(wcp_virt_interactions), intent(out):: self
        class(base_particles), target, intent(in):: ps_lhs
        class(base_particles), target, optional, intent(in):: ps_rhs
        integer, intent(in):: npairs_per_particle
        select type (ps => ps_lhs)
        type is (wcp)
            self%ps_real => ps
        end select
        select type (ps => ps_rhs)
        type is (wcp)
            self%ps_virt => ps
        end select
        call self%base_init(npairs_per_particle, ps_lhs, ps_rhs)
    end subroutine wcp_virt_interactions_init

    subroutine test_compile()
        type(particle_interactions_container):: wcp_interaction_pairs(2)
        type(particles_container):: wcp_sets(2)
        type(grasph_cubic_bspline_kernel):: kernel

        allocate(wcp::wcp_sets(1)%p)
        allocate(wcp::wcp_sets(2)%p)

        allocate(wcp_interactions::wcp_interaction_pairs(1)%pi)
        allocate(wcp_virt_interactions::wcp_interaction_pairs(2)%pi)

        select type (ps => wcp_setS(1)%p)
        type is (wcp)
            call ps%init(n=1, d=2, name="test", rho_ref=1000._fp)
        end select
        select type (ps => wcp_sets(2)%p)
        type is (wcp)
            call ps%init(n=1, d=2, name="test", rho_ref=1000._fp)
        end select

        select type (psi => wcp_interaction_pairs(1)%pi)
        type is (wcp_interactions)
            call psi%init(1, wcp_sets(1)%p)
        end select
        select type (psi => wcp_interaction_pairs(2)%pi)
        type is (wcp_virt_interactions)
            call psi%init(1, wcp_sets(1)%p, wcp_sets(2)%p)
        end select

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