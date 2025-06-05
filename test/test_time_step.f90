module test_time_step

    use grasph_constants, only: fp
    use grasph_kernels, only: grasph_cubic_bspline_kernel
    use grasph_pairs, only: particle_pairs
    use grasph_particles, only: wc_particles => weakly_compressible_particles
    use grasph_pair_sets, only: interacting_particle_set
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

    ! data for testing
    real(fp), parameter:: dx = 0.25_fp

    type, extends(interacting_particle_set):: example_real_virt_set
        class(wc_particles), pointer:: lhs_wcp => null(), rhs_wcp => null()
    contains
        procedure:: sweep_callback => test_sweep_callback
        procedure:: init => example_real_virt_set_init
    end type example_real_virt_set

    type, extends(example_real_virt_set):: example_real_virt_set2
    contains
        procedure:: sweep_callback => test_sweep_callback2
    end type example_real_virt_set2

contains

    type(test_list) function tests()

        tests = test_list([ &
            test("test_set_pair_setup", test_set_pair_setup) &
        ])

    end function tests

    subroutine example_real_virt_set_init(self, lhs, rhs, npairs_per_particle)
        class(example_real_virt_set), intent(inout):: self
        class(wc_particles), intent(in), target:: lhs, rhs
        integer, intent(in):: npairs_per_particle
        call self%base_init(lhs, rhs, npairs_per_particle)
        self%lhs_wcp => lhs
        self%rhs_wcp => rhs
    end subroutine example_real_virt_set_init

    subroutine test_sweep_callback(self, dt)
        class(example_real_virt_set), intent(inout):: self
        real(fp), intent(in), optional:: dt
        integer:: i, j
        do i = 1, self%lhs_wcp%size
            do j = 1, self%rhs_particles%size
                self%lhs_wcp%p(i) = self%lhs_wcp%p(i) + self%rhs_wcp%p(j)
            enddo
        enddo
    end subroutine test_sweep_callback

    subroutine test_sweep_callback2(self, dt)
        class(example_real_virt_set2), intent(inout):: self
        real(fp), intent(in), optional:: dt
        integer:: i, jj, j
        do i = 1, self%pairs%n
            do jj = self%pairs%offsets(i)+1, self%pairs%offsets(i+1)
                j = self%pairs%rhs(jj)
                self%lhs_wcp%p(i) = self%lhs_wcp%p(i) + self%rhs_wcp%p(j)
            enddo
        enddo
    end subroutine test_sweep_callback2

    subroutine test_set_pair_setup()

        type(example_real_virt_set):: real_virt_set
        type(example_real_virt_set2):: real_virt_set2
        type(wc_particles), target:: realp, virtp
        type(grasph_cubic_bspline_kernel):: kernel
        integer:: ii, j, i
        character:: ic
        integer, parameter:: nd = 2, nxr = 2, nr = nxr**nd, nxv = 3, nv = nxv**nd

        call realp%init(nr, 2, 4, 1._fp)
        call virtp%init(nv, 2, 0, 1._fp)
        do i = 1, nr
            realp%p(i) = real(i, kind=fp)
        enddo
        do i = 1, nv
            virtp%p(i) = real(i, kind=fp)
        enddo

        ! manual init
        call real_virt_set%init(realp, virtp, nv)
        call real_virt_set%sweep_callback()

        do i = 1, 4
            write(ic, "(I1)") i
            call check( &
                is_close(realp%p(i), real(i+(nv*(nv+1)/2), kind=fp)), &
                "Incorrect updated pressure of particle " // ic // " using simple sweep" &
            )
        enddo

        call real_virt_set2%init(realp, virtp, 8)
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
            enddo
        enddo
        call real_virt_set2%find_pairs(0.75_fp*dx, kernel)

        call check( &
            is_equal(real_virt_set2%pairs%npairs_total, 16), &
            "Incorrect pairs calculated between real and virtual particles" &
        )

        call real_virt_set2%sweep_callback()

        call check( &
            is_close(real_virt_set2%lhs_wcp%p(1), 13._fp), &
            "Incorrect pressure calculated for particle 1 during second sweep" &
        )
        call check( &
            is_close(real_virt_set2%lhs_wcp%p(2), 18._fp), &
            "Incorrect pressure calculated for particle 2 during second sweep" &
        )
        call check( &
            is_close(real_virt_set2%lhs_wcp%p(3), 27._fp), &
            "Incorrect pressure calculated for particle 3 during second sweep" &
        )
        call check( &
            is_close(real_virt_set2%lhs_wcp%p(4), 32._fp), &
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