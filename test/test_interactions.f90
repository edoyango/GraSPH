module test_interactions

    use grasph_constants, only: fp
    use grasph_kernels, only: grasph_base_kernel, grasph_cubic_bspline_kernel
    use grasph_pairs, only: particle_pairs, cell_list_search
    use grasph_particles, only: base_particles, wc_particles => weakly_compressible_particles
    use grasph_pair_sets, only: particle_interactions_base
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

    ! data for testing
    real(fp), parameter:: dx = 0.25_fp

    type, extends(particle_interactions_base):: example_real_virt_set
        class(wc_particles), pointer:: lhs_wcp => null(), rhs_wcp => null()
    contains
        procedure:: init => example_real_virt_set_init
        procedure:: sweep => example_real_virt_sweep
    end type example_real_virt_set

    type, extends(particle_interactions_base):: example_self_set
        class(wc_particles), pointer:: wcp => null()
    contains
        procedure:: init => example_self_set_init
    end type example_self_set

contains

    type(test_list) function tests()

        tests = test_list([ &
            test("test_set_pair_setup", test_set_pair_setup), &
            test("test_find_self_pairs", test_find_self_pairs) &
        ])

    end function tests

    subroutine example_real_virt_set_init(self, npairs_per_particle, ps_lhs, ps_rhs)
        class(example_real_virt_set), intent(out):: self
        class(base_particles), target, intent(in):: ps_lhs ! base_particles needed to ensure matching interface with overriden init
        class(base_particles), target, intent(in), optional:: ps_rhs
        integer, intent(in):: npairs_per_particle
        ! base_init first as it wipes out self (intent(out))
        call self%base_init(npairs_per_particle, ps_lhs, ps_rhs)
        ! select type to make sure pointer and input align
        select type (ps => ps_lhs)
        class is (wc_particles)
            self%lhs_wcp => ps
        class default
            error stop "Invalid class for ps_lhs"
        end select
        select type (ps => ps_rhs)
        class is (wc_particles)
            self%rhs_wcp => ps
        class default
            error stop "Invalid class for ps_lhs"
        end select
    end subroutine example_real_virt_set_init

    subroutine example_real_virt_sweep(self)
        class(example_real_virt_set), intent(inout):: self
        integer:: i, jj, j
        do i = 1, self%pairs%n
            do jj = self%pairs%offsets(i)+1, self%pairs%offsets(i+1)
                j = self%pairs%rhs(jj)
                self%lhs_wcp%p(i) = self%lhs_wcp%p(i) + self%rhs_wcp%p(j)
            enddo
        enddo
    end subroutine example_real_virt_sweep

    subroutine example_self_set_init(self, npairs_per_particle, ps_lhs, ps_rhs)
        class(example_self_set), intent(out):: self
        class(base_particles), target, intent(in):: ps_lhs
        class(base_particles), target, optional, intent(in):: ps_rhs
        integer, intent(in):: npairs_per_particle
        call self%base_init(npairs_per_particle, ps_lhs)
        select type (ps => ps_lhs)
        class is (wc_particles)
            self%wcp => ps
        end select
    end subroutine example_self_set_init

    subroutine test_set_pair_setup()

        type(example_real_virt_set):: real_virt_set
        type(wc_particles), target:: realp, virtp
        type(grasph_cubic_bspline_kernel):: kernel
        integer:: ii, j, i
        character:: ic
        integer, parameter:: nd = 2, nxr = 2, nr = nxr**nd, nxv = 3, nv = nxv**nd

        call realp%init(nr, 2, 1._fp)
        call virtp%init(nv, 2, 1._fp)
        do concurrent (i=0:nxr-1, j=0:nxr-1)
            ii = i*nxr + j + 1
            realp%x(1, ii) = (i+0.5_fp)*dx
            realp%x(2, ii) = (j+0.5_fp)*dx
            realp%p(ii) = real(ii, kind=fp)
        enddo
        do concurrent (i=0:nxv-1, j=0:nxv-1)
            ii = i*nxv + j + 1
            virtp%x(1, ii) = i*dx
            virtp%x(2, ii) = j*dx
            virtp%p(ii) = real(ii, kind=fp)
        enddo

        ! manual init
        call real_virt_set%init(nv, realp, virtp)
        call real_virt_set%find_pairs(1._fp, kernel)
        call real_virt_set%sweep()

        do i = 1, 4
            write(ic, "(I1)") i
            call check( &
                is_close(realp%p(i), real(i+(nv*(nv+1)/2), kind=fp)), &
                "Incorrect updated pressure of particle " // ic // " using simple sweep" &
            )
        enddo

        ! reset pressures for next test
        realp%p(:) = [(real(i, kind=fp), i = 1, nr)]
        virtp%p(:) = [(real(i, kind=fp), i = 1, nv)]

        call real_virt_set%find_pairs(0.75_fp*dx, kernel)
        call real_virt_set%sweep()

        call check( &
            is_close(real_virt_set%lhs_wcp%p(1), 13._fp), &
            "Incorrect pressure calculated for particle 1 during second sweep" &
        )
        call check( &
            is_close(real_virt_set%lhs_wcp%p(2), 18._fp), &
            "Incorrect pressure calculated for particle 2 during second sweep" &
        )
        call check( &
            is_close(real_virt_set%lhs_wcp%p(3), 27._fp), &
            "Incorrect pressure calculated for particle 3 during second sweep" &
        )
        call check( &
            is_close(real_virt_set%lhs_wcp%p(4), 32._fp), &
            "Incorrect pressure calculated for particle 4 during second sweep" &
        )
        
    end subroutine test_set_pair_setup

    subroutine test_find_self_pairs()

        type(wc_particles):: ps
        type(grasph_cubic_bspline_kernel):: kernel
        type(example_self_set):: ps_set
        integer:: i, j, k, ii

        call ps%init(27, 3, 0._fp)

        do concurrent (i=0:2, j=0:2, k=0:2)
            ii = i*9+j*3+k+1
            ps%x(1, ii) = (i+0.5_fp)*dx
            ps%x(2, ii) = (j+0.5_fp)*dx
            ps%x(3, ii) = (k+0.5_fp)*dx
            ps%p(ii) = real(ii, kind=fp)
        enddo

        call kernel%init(3, 0.9_fp*dx)
        
        call ps_set%init(27, ps)

        call ps_set%find_pairs(kernel%cutoff, kernel)

        ! basic check as correctness checks are in test_pair_finding
        call check( &
            is_equal(ps_set%pairs%npairs_total, 158), &
            "Particles pair finding got wrong number of pairs" &
        )

    end subroutine test_find_self_pairs

end module test_interactions

program run_tests

    use test_interactions, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests