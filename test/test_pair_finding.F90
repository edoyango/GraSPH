module test_pair_finding

    use grasph_constants_m, only: fp
    use grasph_particle_system_m, only: base_particle_t
    use grasph_kernels_m, only: cubic_bspline_kernel_t
    use grasph_pairs_m, only: particle_pairs_t, dsearch, cell_list_search
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

    ! data for testing
    real(fp), parameter:: dx = 0.25_fp
#ifdef THREED
    type(base_particle_t):: ps(27), ps_other(64)
    integer:: pairs3d_1(2, 351), pairs3d_2(2, 158), pairs3d_other(2, 216)
#else
    type(base_particle_t):: ps(16), ps_other(25)
    integer:: pairs2d_1(2, 120), pairs2d_2(2, 42), pairs2d_other(2, 64)
#endif

contains

    type(test_list) function tests()

        call test_setup()

        tests = test_list([ &
                          test("test_dsearch", test_dsearch), &
                          test("test_dsearch_other", test_dsearch_other), &
                          test("test_cell_list", test_cell_list), &
                          test("test_cell_list_other", test_cell_list_other) &
                          ])

    end function tests

    subroutine test_setup()

        integer:: i, j, k, ii, n

#ifdef THREED

        do concurrent(i=0:2, j=0:2, k=0:2)
            ii = i*9 + j*3 + k + 1
            ps(ii)%x(1) = (i + 0.5_fp)*dx
            ps(ii)%x(2) = (j + 0.5_fp)*dx
            ps(ii)%x(3) = (k + 0.5_fp)*dx
        end do

        do concurrent(i=0:3, j=0:3, k=0:3)
            ii = i*16 + j*4 + k + 1
            ps_other(ii)%x(1) = i*dx
            ps_other(ii)%x(2) = j*dx
            ps_other(ii)%x(3) = k*dx
        end do

        n = 0
        do i = 1, 26
            do j = i + 1, 27
                n = n + 1
                pairs3d_1(1, n) = i
                pairs3d_1(2, n) = j
            end do
        end do

        n = 0
        do i = 1, 26
            do j = i + 1, 27
                if (sum((x3d(:, i) - x3d(:, j))**2) < (dx*1.8_fp)**2) then
                    n = n + 1
                    pairs3d_2(1, n) = i
                    pairs3d_2(2, n) = j
                end if
            end do
        end do

        n = 0
        do i = 1, 27
            do j = 1, 64
                if (sum((x3d(:, i) - x3d_other(:, j))**2) < (sqrt(0.8_fp)*dx)**2) then
                    n = n + 1
                    pairs3d_other(1, n) = i
                    pairs3d_other(2, n) = j
                end if
            end do
        end do

#else
        ! particles on grid such that x, y ∈ (0, 1)
        do concurrent(i=0:3, j=0:3)
            ii = i*4 + j + 1
            ps(ii)%x(1) = (i + 0.5_fp)*dx
            ps(ii)%x(2) = (j + 0.5_fp)*dx
        end do

        ! "                         " x, y ∈ [0, 1]
        do concurrent(i=0:4, j=0:4)
            ii = i*5 + j + 1
            ps_other(ii)%x(1) = i*dx
            ps_other(ii)%x(2) = j*dx
        end do

        n = 0
        do i = 1, 15
            do j = i + 1, 16
                n = n + 1
                pairs2d_1(1, n) = i
                pairs2d_1(2, n) = j
            end do
        end do

        n = 0
        do i = 1, 15
            do j = i + 1, 16
                if (sum((ps(i)%x(:) - ps(j)%x(:))**2) < (dx*1.8_fp)**2) then
                    n = n + 1
                    pairs2d_2(1, n) = i
                    pairs2d_2(2, n) = j
                end if
            end do
        end do

        ! calculate pairs for finding pairs between two sets of particles
        n = 0
        do i = 1, 16
            do j = 1, 25
                if (sum((ps(i)%x(:) - ps_other(j)%x(:))**2) < (sqrt(0.6_fp)*dx)**2) then
                    n = n + 1
                    pairs2d_other(1, n) = i
                    pairs2d_other(2, n) = j
                end if
            end do
        end do

#endif

    end subroutine test_setup

    subroutine check_pairs(pairs, ps_lhs, ps_rhs, correct_pairs, ncorrect_pairs, case_string, kernel)

        type(particle_pairs_t), intent(in):: pairs
        type(base_particle_t), intent(in):: ps_lhs(:), ps_rhs(:)
        integer, intent(in):: ncorrect_pairs, correct_pairs(2, ncorrect_pairs)
        character(*), intent(in):: case_string
        type(cubic_bspline_kernel_t), intent(in):: kernel
        integer:: n, i, j, k, jj, ii
        character(3):: c, ic, jc
        real(fp):: w, dwdx(kernel%d)

        call check( &
            is_equal(pairs%npairs_total, ncorrect_pairs), &
            "Not all pairs found!" &
            )

        n = 0
        sweep_pairs_to_check: do k = 1, pairs%npairs_total
            i = pairs%pair_ij(1, k)
            j = pairs%pair_ij(2, k)
            do ii = 1, ncorrect_pairs
                if ((correct_pairs(1, ii) == i .and. correct_pairs(2, ii) == j) .or. &
                    (correct_pairs(1, ii) == j .and. correct_pairs(2, ii) == i)) then
                    ! just to make sure the check index is updated properly
                    call check(.true., "")
                    ! checking the kernel values in pairs have been updated properly
                    call kernel%values(ps_lhs(i)%x(:) - ps_rhs(j)%x(:), w, dwdx)
                    call check( &
                        is_close(w, pairs%w(k)), &
                        "Case: "//case_string//". Incorrectly updated w value: "//ic//", j :"//jc &
                        )
                    cycle sweep_pairs_to_check
                end if
            end do
            call check(.false., "Case: "//case_string//". Couldn't find pair - i: "//ic//", j: "//jc)
        end do sweep_pairs_to_check

    end subroutine check_pairs

    subroutine test_dsearch()

        type(particle_pairs_t):: pairs
        type(cubic_bspline_kernel_t):: kernel

#ifdef THREED

        call kernel%init(3, 1._fp)
        call pairs%init(27, 26)

        ! 3d - cutoff selected for all particles to be paired with eachother
        call dsearch(ps, sqrt(3._fp), kernel, pairs)
        call check_pairs(pairs, ps, ps, pairs3d_1, 351, "brute-force (3d - all pairs)", kernel)

        ! 3d - cutoff selected so pairs are adjacent (incl. diagonal) on grid
        call dsearch(ps, dx*1.8_fp, kernel, pairs)
        call check_pairs(pairs, ps, ps, pairs3d_2, 158, "brute-force (3d - adj pairs)", kernel)

#else

        call kernel%init(2, 1._fp)
        call pairs%init(16, 15)

        ! 2d - cutoff selected for all particles to be paired with eachother
        call dsearch(ps, sqrt(2._fp), kernel, pairs)
        call check_pairs(pairs, ps, ps, pairs2d_1, 120, "brute-force (2d - all pairs)", kernel)

        ! 2d -cutoff selected so pairs are adjacent (incl. diagonal) on grid
        call dsearch(ps, dx*1.8_fp, kernel, pairs)
        call check_pairs(pairs, ps, ps, pairs2d_2, 42, "brute-force (2d - adj pairs)", kernel)
#endif

    end subroutine test_dsearch

    subroutine test_dsearch_other()

        type(particle_pairs_t):: pairs
        type(cubic_bspline_kernel_t):: kernel
#ifdef THREED
        call kernel%init(3, 1._fp)
        call pairs%init(27, 8)

        call dsearch(ps, ps_other, dx, kernel, pairs)
        call check_pairs(pairs, ps, ps_other, pairs3d_other, 216, "brute-force (3d - 2sets)", kernel)
#else
        call kernel%init(2, 1._fp)
        call pairs%init(16, 4)

        call dsearch(ps, ps_other, dx, kernel, pairs)
        call check_pairs(pairs, ps, ps_other, pairs2d_other, 64, "brute-force (2d - 2sets)", kernel)
#endif
    end subroutine test_dsearch_other

    subroutine test_cell_list()

        type(particle_pairs_t):: pairs
        type(cubic_bspline_kernel_t):: kernel
#ifdef THREED
        call kernel%init(3, 1._fp)
        call pairs%init(27, 27)

        ! 3d - cutoff selected for all particles to be paired with eachother
        call cell_list_search(ps, sqrt(3._fp), kernel, pairs)
        call check_pairs(pairs, ps, ps, pairs3d_1, 351, "cell-lists (3d - all pairs)", kernel)

        ! 3d - cutoff selected so pairs are adjacent (incl. diagonal) on grid
        call cell_list_search(ps, dx*1.8_fp, kernel, pairs)
        call check_pairs(pairs, ps, ps, pairs3d_2, 158, "cell-lists (3d - adj pairs)", kernel)
#else
        call kernel%init(2, 1._fp)
        call pairs%init(16, 16)

        ! 2d - cutoff selected for all particles to be paired with eachother
        call cell_list_search(ps, sqrt(2._fp), kernel, pairs)
        call check_pairs(pairs, ps, ps, pairs2d_1, 120, "cell-lists (2d - all pairs)", kernel)

        ! 2d -cutoff selected so pairs are adjacent (incl. diagonal) on grid
        call cell_list_search(ps, dx*1.8_fp, kernel, pairs)
        call check_pairs(pairs, ps, ps, pairs2d_2, 42, "cell-lists (3d - adj pairs)", kernel)
#endif
    end subroutine test_cell_list

    subroutine test_cell_list_other()

        type(particle_pairs_t):: pairs
        type(cubic_bspline_kernel_t):: kernel
#ifdef THREED
        call kernel%init(3, 1._fp)
        call pairs%init(27, 8)

        call cell_list_search(ps, ps_other, dx, kernel, pairs)
        call check_pairs(pairs, ps, ps_other, pairs3d_other, 216, "cell-lists (3d - 2sets)", kernel)
#else
        call kernel%init(2, 1._fp)
        call pairs%init(16, 4)

        call cell_list_search(ps, ps_other, sqrt(0.6_fp)*dx, kernel, pairs)
        call check_pairs(pairs, ps, ps_other, pairs2d_other, 64, "cell-lists (2d - 2sets)", kernel)
#endif

    end subroutine test_cell_list_other

end module test_pair_finding

program run_tests

    use test_pair_finding, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
