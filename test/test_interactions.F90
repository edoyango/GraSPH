module test_interactions

    use grasph_constants_m, only: fp, max_name_len
    use grasph_kernels_m, only: base_kernel_t, cubic_bspline_kernel_t
    use grasph_pairs_m, only: particle_pairs_t, cell_list_search
    use grasph_particle_system_m, only: particle_system_t, base_particles_t
    use weakly_compressible_particles_m, only: eos_particles_t
    use grasph_system_interactions_m, only: system_interaction_t, base_sweeper_t, sweeper_container_t
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

    ! data for testing
    real(fp), parameter:: dx = 0.25_fp

    type, extends(base_sweeper_t):: example_real_virt_sweeper_t
    contains
        procedure:: sweep_1system => example_real_virt_sweep_1system
        procedure:: sweep_2system => example_real_virt_sweep_2system
        procedure:: sweep_2system_norhsupdate => example_real_virt_sweep_2system
        procedure, nopass:: name => example_sweeper_name
    end type example_real_virt_sweeper_t

contains

    type(test_list) function tests()

        tests = test_list([ &
                          test("test_set_pair_setup", test_set_pair_setup), &
                          test("test_find_self_pairs", test_find_self_pairs) &
                          ])

    end function tests

    pure character(max_name_len) function example_sweeper_name()
        example_sweeper_name = "example_sweeper_t"
    end function example_sweeper_name

    subroutine example_real_virt_sweep_1system(self, pairs, ps, dt)

        class(example_real_virt_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps
        real(fp), optional, intent(in):: dt

        error stop "Cannot perform sweep with only 1 particle system. Ensure that both ps_lhs and ps_rhs are associated."

    end subroutine example_real_virt_sweep_1system

    subroutine example_real_virt_sweep_2system(self, pairs, ps_lhs, ps_rhs, dt)
        class(example_real_virt_sweeper_t), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(base_particles_t), intent(inout):: ps_lhs, ps_rhs
        real(fp), optional, intent(in):: dt
        integer:: i, j, k
        class(eos_particles_t), pointer:: ps_real, ps_virt

        ! assign pointers to ps_lhs/rhs for access to pressure
        select type (ps => ps_lhs)
        class is (eos_particles_t)
            ps_real => ps
        class default
            error stop "Invalid class for ps_lhs"
        end select

        select type (ps => ps_rhs)
        class is (eos_particles_t)
            ps_virt => ps
        class default
            error stop "Invalid class for ps_rhs"
        end select

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            ps_real%p(i) = ps_real%p(i) + ps_virt%p(j)
        end do

    end subroutine example_real_virt_sweep_2system

    subroutine test_set_pair_setup()

        type(system_interaction_t):: real_virt_set
        type(sweeper_container_t):: sweepers(1)
        type(particle_system_t), target:: psys_real, psys_virt
        type(cubic_bspline_kernel_t):: kernel
        integer:: ii, j, i
        character:: ic
        integer, parameter:: nd = 2, nxr = 2, nr = nxr**nd, nxv = 3, nv = nxv**nd
        class(eos_particles_t), pointer:: ps_lhs, ps_rhs
        type(eos_particles_t):: ps_template

        allocate (example_real_virt_sweeper_t::sweepers(1)%sweeper)

#ifndef THREED
        call psys_real%init(n=nr, name="test", particle_template=ps_template)
        call psys_virt%init(n=nv, name="test", particle_template=ps_template)
        select type (ps => psys_real%particles)
        class is (eos_particles_t)
            ps_lhs => ps
        class default
            error stop "Expected eos_particles_t for psys_real%particles."
        end select
        do concurrent(i=0:nxr - 1, j=0:nxr - 1)
            ii = i*nxr + j + 1
            ps_lhs%x(1, ii) = (i + 0.5_fp)*dx
            ps_lhs%x(2, ii) = (j + 0.5_fp)*dx
            ps_lhs%p(ii) = real(ii, kind=fp)
        end do
        select type (ps => psys_virt%particles)
        class is (eos_particles_t)
            ps_rhs => ps
        class default
            error stop "Expected eos_particles_t for virt%particles."
        end select
        do concurrent(i=0:nxv - 1, j=0:nxv - 1)
            ii = i*nxv + j + 1
            ps_rhs%x(1, ii) = i*dx
            ps_rhs%x(2, ii) = j*dx
            ps_rhs%p(ii) = real(ii, kind=fp)
        end do

        ! manual init
        call real_virt_set%init(nv, psys_real, psys_virt, sweepers=sweepers)
        call real_virt_set%find_pairs(1._fp, kernel)
        call real_virt_set%do_sweep(1)

        do i = 1, 4
            write (ic, "(I1)") i
            call check( &
                is_close(ps_lhs%p(i), real(i + (nv*(nv + 1)/2), kind=fp)), &
                "Incorrect updated pressure of particle "//ic//" using simple sweep" &
                )
        end do

        ! reset pressures for next test
        do i = 1, nr
            ps_lhs%p(i) = real(i, kind=fp)
        end do
        do i = 1, nv
            ps_rhs%p(i) = real(i, kind=fp)
        end do

        call real_virt_set%find_pairs(0.75_fp*dx, kernel)
        call real_virt_set%do_sweep(1)

        call check( &
            is_close(ps_lhs%p(1), 13._fp), &
            "Incorrect pressure calculated for particle 1 during second sweep" &
            )
        call check( &
            is_close(ps_lhs%p(2), 18._fp), &
            "Incorrect pressure calculated for particle 2 during second sweep" &
            )
        call check( &
            is_close(ps_lhs%p(3), 27._fp), &
            "Incorrect pressure calculated for particle 3 during second sweep" &
            )
        call check( &
            is_close(ps_lhs%p(4), 32._fp), &
            "Incorrect pressure calculated for particle 4 during second sweep" &
            )
#endif

    end subroutine test_set_pair_setup

    subroutine test_find_self_pairs()

        type(particle_system_t):: psys
        type(cubic_bspline_kernel_t):: kernel
        type(system_interaction_t):: ps_set
        integer:: i, j, k, ii
        class(eos_particles_t), pointer:: ps_real(:)

#ifdef THREED

        select type (ps => psys%particles)
        class is (eos_particles_t)
            ps_real => ps
        class default
            error stop "Expected eos_particles_t for psys%particles"
        end select

        call psys%init(27, "test", 0._fp)

        do concurrent(i=0:2, j=0:2, k=0:2)
            ii = i*9 + j*3 + k + 1
            ps_real(ii)%x(1) = (i + 0.5_fp)*dx
            ps_real(ii)%x(2) = (j + 0.5_fp)*dx
            ps_real(ii)%x(3) = (k + 0.5_fp)*dx
            ps_real(ii)%p = real(ii, kind=fp)
        end do

        call kernel%init(3, 0.9_fp*dx)

        call ps_set%init(27, psys)

        call ps_set%find_pairs(kernel%cutoff, kernel)

        ! basic check as correctness checks are in test_pair_finding
        call check( &
            is_equal(ps_set%pairs%npairs_total, 158), &
            "Particles pair finding got wrong number of pairs" &
            )

#endif

    end subroutine test_find_self_pairs

end module test_interactions

program run_tests

    use test_interactions, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
