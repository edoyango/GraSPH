module test_pair_finding

    use grasph_constants, only: fp
    use grasph_kernels, only: grasph_cubic_bspline_kernel
    use grasph_pairs, only: particle_pairs, dsearch, cell_list_search
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

    ! data for testing
    real(fp), parameter:: dx = 0.25_fp
    real(fp):: x2d(2, 16), x3d(3, 27), x2d_other(2, 25), x3d_other(3, 64)
    integer:: pairs2d_1(2, 120), pairs2d_2(2, 42), pairs3d_1(2, 351), pairs3d_2(2, 158)
    integer:: pairs2d_other(2, 64), pairs3d_other(2, 216)

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

        ! particles on grid such that x, y ∈ (0, 1)
        do concurrent(i=0:3, j=0:3)
            ii = i*4 + j + 1
            x2d(1, ii) = (i + 0.5_fp)*dx
            x2d(2, ii) = (j + 0.5_fp)*dx
        end do

        ! "                         " x, y ∈ [0, 1]
        do concurrent(i=0:4, j=0:4)
            ii = i*5 + j + 1
            x2d_other(1, ii) = i*dx
            x2d_other(2, ii) = j*dx
        end do

        do concurrent(i=0:2, j=0:2, k=0:2)
            ii = i*9 + j*3 + k + 1
            x3d(1, ii) = (i + 0.5_fp)*dx
            x3d(2, ii) = (j + 0.5_fp)*dx
            x3d(3, ii) = (k + 0.5_fp)*dx
        end do

        do concurrent(i=0:3, j=0:3, k=0:3)
            ii = i*16 + j*4 + k + 1
            x3d_other(1, ii) = i*dx
            x3d_other(2, ii) = j*dx
            x3d_other(3, ii) = k*dx
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
                if (sum((x2d(:, i) - x2d(:, j))**2) < (dx*1.8_fp)**2) then
                    n = n + 1
                    pairs2d_2(1, n) = i
                    pairs2d_2(2, n) = j
                end if
            end do
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

        ! calculate pairs for finding pairs between two sets of particles
        n = 0
        do i = 1, 16
            do j = 1, 25
                if (sum((x2d(:, i) - x2d_other(:, j))**2) < (sqrt(0.6_fp)*dx)**2) then
                    n = n + 1
                    pairs2d_other(1, n) = i
                    pairs2d_other(2, n) = j
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

    end subroutine test_setup

    subroutine check_pairs(pairs, x_lhs, x_rhs, correct_pairs, ncorrect_pairs, case_string, kernel)

        type(particle_pairs), intent(in):: pairs
        real(fp), intent(in):: x_lhs(:, :), x_rhs(:, :)
        integer, intent(in):: ncorrect_pairs, correct_pairs(2, ncorrect_pairs)
        character(*), intent(in):: case_string
        type(grasph_cubic_bspline_kernel), intent(in):: kernel
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
                    call kernel%values(x_lhs(:, i) - x_rhs(:, j), w, dwdx)
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

        type(particle_pairs):: pairs
        type(grasph_cubic_bspline_kernel):: kernel

        call kernel%init(2, 1._fp)
        call pairs%init(16, 15, 2)

        ! 2d - cutoff selected for all particles to be paired with eachother
        call dsearch(x2d, sqrt(2._fp), kernel, pairs)
        call check_pairs(pairs, x2d, x2d, pairs2d_1, 120, "brute-force (2d - all pairs)", kernel)

        ! 2d -cutoff selected so pairs are adjacent (incl. diagonal) on grid
        call dsearch(x2d, dx*1.8_fp, kernel, pairs)
        call check_pairs(pairs, x2d, x2d, pairs2d_2, 42, "brute-force (2d - adj pairs)", kernel)

        call kernel%init(3, 1._fp)
        call pairs%init(27, 26, 3)

        ! 3d - cutoff selected for all particles to be paired with eachother
        call dsearch(x3d, sqrt(3._fp), kernel, pairs)
        call check_pairs(pairs, x3d, x3d, pairs3d_1, 351, "brute-force (3d - all pairs)", kernel)

        ! 3d - cutoff selected so pairs are adjacent (incl. diagonal) on grid
        call dsearch(x3d, dx*1.8_fp, kernel, pairs)
        call check_pairs(pairs, x3d, x3d, pairs3d_2, 158, "brute-force (3d - adj pairs)", kernel)

    end subroutine test_dsearch

    subroutine test_dsearch_other()

        type(particle_pairs):: pairs
        type(grasph_cubic_bspline_kernel):: kernel
        call kernel%init(2, 1._fp)
        call pairs%init(16, 4, 2)

        call dsearch(x2d, x2d_other, 25, dx, kernel, pairs)
        call check_pairs(pairs, x2d, x2d_other, pairs2d_other, 64, "brute-force (2d - 2sets)", kernel)

        call kernel%init(3, 1._fp)
        call pairs%init(27, 8, 3)

        call dsearch(x3d, x3d_other, 64, dx, kernel, pairs)
        call check_pairs(pairs, x3d, x3d_other, pairs3d_other, 216, "brute-force (3d - 2sets)", kernel)

    end subroutine test_dsearch_other

    subroutine test_cell_list()

        type(particle_pairs):: pairs
        type(grasph_cubic_bspline_kernel):: kernel

        call kernel%init(2, 1._fp)
        call pairs%init(16, 16, 2)

        ! 2d - cutoff selected for all particles to be paired with eachother
        call cell_list_search(x2d, sqrt(2._fp), kernel, pairs)
        call check_pairs(pairs, x2d, x2d, pairs2d_1, 120, "cell-lists (2d - all pairs)", kernel)

        ! 2d -cutoff selected so pairs are adjacent (incl. diagonal) on grid
        call cell_list_search(x2d, dx*1.8_fp, kernel, pairs)
        call check_pairs(pairs, x2d, x2d, pairs2d_2, 42, "cell-lists (3d - adj pairs)", kernel)

        call kernel%init(3, 1._fp)
        call pairs%init(27, 27, 3)

        ! 3d - cutoff selected for all particles to be paired with eachother
        call cell_list_search(x3d, sqrt(3._fp), kernel, pairs)
        call check_pairs(pairs, x3d, x3d, pairs3d_1, 351, "cell-lists (3d - all pairs)", kernel)

        ! 3d - cutoff selected so pairs are adjacent (incl. diagonal) on grid
        call cell_list_search(x3d, dx*1.8_fp, kernel, pairs)
        call check_pairs(pairs, x3d, x3d, pairs3d_2, 158, "cell-lists (3d - adj pairs)", kernel)

    end subroutine test_cell_list

    subroutine test_cell_list_other()

        type(particle_pairs):: pairs
        type(grasph_cubic_bspline_kernel):: kernel
        call kernel%init(2, 1._fp)
        call pairs%init(16, 4, 2)

        call cell_list_search(x2d, x2d_other, 25, sqrt(0.6_fp)*dx, kernel, pairs)
        call check_pairs(pairs, x2d, x2d_other, pairs2d_other, 64, "cell-lists (2d - 2sets)", kernel)

        call kernel%init(3, 1._fp)
        call pairs%init(27, 8, 3)

        call cell_list_search(x3d, x3d_other, 64, dx, kernel, pairs)
        call check_pairs(pairs, x3d, x3d_other, pairs3d_other, 216, "cell-lists (3d - 2sets)", kernel)

    end subroutine test_cell_list_other

end module test_pair_finding

program run_tests

    use test_pair_finding, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
