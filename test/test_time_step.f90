module test_time_step

    use grasph_constants, only: fp
    use grasph_kernels, only: grasph_base_kernel, grasph_cubic_bspline_kernel
    use grasph_pairs, only: particle_pairs, cell_list_search
    use grasph_particles, only: wc_particles => weakly_compressible_particles
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
        procedure:: find_pairs => example_real_virt_find_pairs
        procedure:: init => example_real_virt_set_init
        procedure:: sweep => example_real_virt_sweep
    end type example_real_virt_set

contains

    subroutine example_real_virt_find_pairs(self, cutoff, kernel)
        class(example_real_virt_set), intent(inout):: self
        real(fp), intent(in):: cutoff
        class(grasph_base_kernel), intent(in):: kernel
        call cell_list_search(self%lhs_wcp%x, self%rhs_wcp%x, self%rhs_wcp%size, cutoff, kernel, self%pairs)
    end subroutine example_real_virt_find_pairs

    subroutine example_real_virt_set_init(self, lhs_wcp, rhs_wcp, npairs_per_particle)
        class(example_real_virt_set), intent(inout):: self
        class(wc_particles), target, intent(in):: lhs_wcp, rhs_wcp
        integer, intent(in):: npairs_per_particle
        self%lhs_wcp => lhs_wcp
        self%rhs_wcp => rhs_wcp
        call self%pairs%init(lhs_wcp%size, npairs_per_particle, lhs_wcp%ndims)
        self%initialized = .true.
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

    type(test_list) function tests()

        tests = test_list([ &
            test("test_set_pair_setup", test_set_pair_setup) &
        ])

    end function tests

    subroutine test_set_pair_setup()

        type(example_real_virt_set):: real_virt_set
        type(wc_particles), target:: realp, virtp
        type(grasph_cubic_bspline_kernel):: kernel
        integer:: ii, j, i
        character:: ic
        integer, parameter:: nd = 2, nxr = 2, nr = nxr**nd, nxv = 3, nv = nxv**nd

        call realp%init(nr, 2, nv, 1._fp)
        call virtp%init(nv, 2, 0, 1._fp)
        do i = 0, nxr-1
            do j = 0, nxr-1
                ii = i*nxr + j + 1
                realp%x(1, ii) = (i+0.5_fp)*dx
                realp%x(2, ii) = (j+0.5_fp)*dx
                realp%p(ii) = real(ii, kind=fp)
            enddo
        enddo
        do i = 0, nxv-1
            do j = 0, nxv-1
                ii = i*nxv + j + 1
                virtp%x(1, ii) = i*dx
                virtp%x(2, ii) = j*dx
                virtp%p(ii) = real(ii, kind=fp)
            enddo
        enddo

        ! manual init
        call real_virt_set%init(realp, virtp, nv)
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

end module test_time_step

program run_tests

    use test_time_step, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests