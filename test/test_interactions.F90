module test_interactions

    use grasph_constants_m, only: fp
    use grasph_kernels_m, only: base_kernel_t, cubic_bspline_kernel_t
    use grasph_pairs_m, only: particle_pairs_t, cell_list_search
    use grasph_particle_system_m, only: particle_system_t
    use weakly_compressible_particles, only: eos_particle
    use grasph_system_interactions_m, only: system_interaction_t, base_sweeper
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

    ! data for testing
    real(fp), parameter:: dx = 0.25_fp

    type, extends(base_sweeper):: example_real_virt_sweeper
    contains
        procedure:: sweep => example_real_virt_sweep
    end type example_real_virt_sweeper

contains

    type(test_list) function tests()

        tests = test_list([ &
                          test("test_set_pair_setup", test_set_pair_setup), &
                          test("test_find_self_pairs", test_find_self_pairs) &
                          ])

    end function tests

    subroutine example_real_virt_sweep(self, pairs, psys_lhs, psys_rhs)
        class(example_real_virt_sweeper), intent(in):: self
        type(particle_pairs_t), intent(in):: pairs
        class(particle_system_t), intent(inout):: psys_lhs
        class(particle_system_t), optional, intent(inout):: psys_rhs
        integer:: i, j, k
        class(eos_particle), pointer:: ps_real(:), ps_virt(:)

        ! assign pointers to ps_lhs/rhs for access to pressure
        select type (ps => psys_lhs%particles)
        class is (eos_particle)
            ps_real => ps
        class default
            error stop "Invalid class for psys_lhs"
        end select

        select type (ps => psys_rhs%particles)
        class is (eos_particle)
            ps_virt => ps
        class default
            error stop "Invalid class for psys_rhs"
        end select

        ! perform sweep
        do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            ps_real(i)%p = ps_real(i)%p + ps_virt(j)%p
        end do

    end subroutine example_real_virt_sweep

    subroutine test_set_pair_setup()

        type(system_interaction_t):: real_virt_set
        type(example_real_virt_sweeper):: rv_sweeper
        type(particle_system_t), target:: psys_real, psys_virt
        type(cubic_bspline_kernel_t):: kernel
        integer:: ii, j, i
        character:: ic
        integer, parameter:: nd = 2, nxr = 2, nr = nxr**nd, nxv = 3, nv = nxv**nd
        class(eos_particle), pointer:: ps_lhs(:), ps_rhs(:)
        type(eos_particle):: ps_template

#ifndef THREED
        call psys_real%base_init(n=nr, name="test", particle_template=ps_template)
        call psys_virt%base_init(n=nv, name="test", particle_template=ps_template)
        select type (ps => psys_real%particles)
        class is (eos_particle)
            ps_lhs => ps
        class default
            error stop "Expected eos_particle for psys_real%particles."
        end select
        do concurrent(i=0:nxr - 1, j=0:nxr - 1)
            ii = i*nxr + j + 1
            ps_lhs(ii)%x(1) = (i + 0.5_fp)*dx
            ps_lhs(ii)%x(2) = (j + 0.5_fp)*dx
            ps_lhs(ii)%p = real(ii, kind=fp)
        end do
        select type (ps => psys_virt%particles)
        class is (eos_particle)
            ps_rhs => ps
        class default
            error stop "Expected eos_particle for virt%particles."
        end select
        do concurrent(i=0:nxv - 1, j=0:nxv - 1)
            ii = i*nxv + j + 1
            ps_rhs(ii)%x(1) = i*dx
            ps_rhs(ii)%x(2) = j*dx
            ps_rhs(ii)%p = real(ii, kind=fp)
        end do

        ! manual init
        call real_virt_set%init(nv, psys_real, psys_virt, sweeper=rv_sweeper)
        call real_virt_set%find_pairs(1._fp, kernel)
        call real_virt_set%do_sweep()

        do i = 1, 4
            write (ic, "(I1)") i
            call check( &
                is_close(ps_lhs(i)%p, real(i + (nv*(nv + 1)/2), kind=fp)), &
                "Incorrect updated pressure of particle "//ic//" using simple sweep" &
                )
        end do

        ! reset pressures for next test
        do i = 1, nr
            ps_lhs(i)%p = real(i, kind=fp)
        end do
        do i = 1, nv
            ps_rhs(i)%p = real(i, kind=fp)
        end do

        call real_virt_set%find_pairs(0.75_fp*dx, kernel)
        call real_virt_set%do_sweep()

        call check( &
            is_close(ps_lhs(1)%p, 13._fp), &
            "Incorrect pressure calculated for particle 1 during second sweep" &
            )
        call check( &
            is_close(ps_lhs(2)%p, 18._fp), &
            "Incorrect pressure calculated for particle 2 during second sweep" &
            )
        call check( &
            is_close(ps_lhs(3)%p, 27._fp), &
            "Incorrect pressure calculated for particle 3 during second sweep" &
            )
        call check( &
            is_close(ps_lhs(4)%p, 32._fp), &
            "Incorrect pressure calculated for particle 4 during second sweep" &
            )
#endif

    end subroutine test_set_pair_setup

    subroutine test_find_self_pairs()

        type(particle_system_t):: psys
        type(cubic_bspline_kernel_t):: kernel
        type(system_interaction_t):: ps_set
        integer:: i, j, k, ii
        class(eos_particle), pointer:: ps_real(:)

#ifdef THREED

        select type (ps => psys%particles)
        class is (eos_particle)
            ps_real => ps
        class default
            error stop "Expected eos_particle for psys%particles"
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
