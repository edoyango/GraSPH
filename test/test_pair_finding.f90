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
    real(fp):: x2d(2, 16), x3d(3, 27)
    integer:: pairs2d_1(2, 120), pairs2d_2(2, 42), pairs3d_1(2, 351), pairs3d_2(2, 158)
    type(grasph_cubic_bspline_kernel):: kernel


contains

    type(test_list) function tests()

        call test_setup()

        tests = test_list([ &
            test("test_dsearch", test_dsearch), &
            test("test_cell_list", test_cell_list) &
        ])

    end function tests

    subroutine test_setup()

        integer:: i, j, k, ii, n

        ! particles on grid such that x, y ∈ (0, 1)
        do i = 0, 3
            do j = 0, 3
                ii = i*4+j+1
                x2d(1, ii) = (i+0.5_fp)*dx
                x2d(2, ii) = (j+0.5_fp)*dx
            enddo
        enddo

        do i = 0, 2
            do j = 0, 2
                do k = 0, 2
                    ii = i*9+j*3+k+1
                    x3d(1, ii) = (i+0.5_fp)*dx
                    x3d(2, ii) = (j+0.5_fp)*dx
                    x3d(3, ii) = (k+0.5_fp)*dx
                enddo
            enddo
        enddo

        call kernel%init(2, 1._fp)

        n = 0
        do i = 1, 15
            do j = i+1, 16
                n = n + 1
                pairs2d_1(1, n) = i
                pairs2d_1(2, n) = j
            enddo
        enddo

        n = 0
        do i = 1, 15
            do j = i+1, 16
                if (sum((x2d(:, i) - x2d(:, j))**2) < (dx*1.8_fp)**2) then
                    n = n + 1
                    pairs2d_2(1, n) = i
                    pairs2d_2(2, n) = j
                endif
            enddo
        enddo

        n = 0
        do i = 1, 26
            do j = i+1, 27
                n = n + 1
                pairs3d_1(1, n) = i
                pairs3d_1(2, n) = j
            enddo
        enddo

        n = 0
        do i = 1, 26
            do j = i+1, 27
                if (sum((x3d(:, i) - x3d(:, j))**2) < (dx*1.8_fp)**2) then
                    n = n + 1
                    pairs3d_2(1, n) = i
                    pairs3d_2(2, n) = j
                endif
            enddo
        enddo

    end subroutine test_setup

    subroutine check_pairs(pairs, correct_pairs, ncorrect_pairs, case_string)

        type(particle_pairs), intent(in):: pairs
        integer, intent(in):: ncorrect_pairs, correct_pairs(2, ncorrect_pairs)
        character(*), intent(in):: case_string
        integer:: n, i, j, jj, ii
        character(3):: c, ic, jc

        call check( &
            is_equal(pairs%npairs_total, ncorrect_pairs), &
            "Not all pairs found!" &
        )

        n = 0
        do i = 1, pairs%n
            i_rhs_sweep: do jj = pairs%offsets(i)+1, pairs%offsets(i+1)
                j = pairs%rhs(jj)
                n = n + 1
                write(c, "(I3)") n
                do ii = 1, ncorrect_pairs
                    if ((correct_pairs(1, ii) == i .and. correct_pairs(2, ii) == j) .or. &
                        (correct_pairs(1, ii) == j .and. correct_pairs(2, ii) == i)) then
                        call check(.true., "")
                        cycle i_rhs_sweep
                    endif
                enddo
                write(ic, "(I3)") i
                write(jc, "(I3)") j
                call check(.false., "Case: " // case_string // ". Couldn't find pair - i: " // ic // ", j: " // jc)
            enddo i_rhs_sweep
        enddo

    end subroutine check_pairs

    subroutine test_dsearch()

        type(particle_pairs):: pairs
        integer:: n, i, j, jj

        ! 2d - cutoff selected for all particles to be paired with eachother
        pairs = dsearch(x2d, 2, 16, sqrt(2._fp), kernel, 15)
        call check_pairs(pairs, pairs2d_1, 120, "brute-force (2d - all pairs)")

        ! 2d -cutoff selected so pairs are adjacent (incl. diagonal) on grid
        pairs = dsearch(x2d, 2, 16, dx*1.8_fp, kernel, 8)
        call check_pairs(pairs, pairs2d_2, 42, "brute-force (2d - adj pairs)")

        ! 3d - cutoff selected for all particles to be paired with eachother
        pairs = dsearch(x3d, 3, 27, sqrt(3._fp), kernel, 26)
        call check_pairs(pairs, pairs3d_1, 351, "brute-force (3d - all pairs)")

        ! 3d - cutoff selected so pairs are adjacent (incl. diagonal) on grid
        pairs = dsearch(x3d, 3, 27, dx*1.8_fp, kernel, 8)
        call check_pairs(pairs, pairs3d_2, 158, "brute-force (3d - adj pairs)")

    end subroutine test_dsearch
    
    subroutine test_cell_list()

        type(particle_pairs):: pairs
        integer:: n, i, j, jj
        character(2):: c

        ! 2d - cutoff selected for all particles to be paired with eachother
        pairs = cell_list_search(x2d, 2, 16, sqrt(2._fp), kernel, 16)
        call check_pairs(pairs, pairs2d_1, 120, "cell-lists (2d - all pairs)")

        ! 2d -cutoff selected so pairs are adjacent (incl. diagonal) on grid
        pairs = cell_list_search(x2d, 2, 16, dx*1.8_fp, kernel, 8)
        call check_pairs(pairs, pairs2d_2, 42, "cell-lists (3d - adj pairs)")

        ! 3d - cutoff selected for all particles to be paired with eachother
        pairs = cell_list_search(x3d, 3, 27, sqrt(3._fp), kernel, 27)
        call check_pairs(pairs, pairs3d_1, 351, "cell-lists (3d - all pairs)")

        ! 3d - cutoff selected so pairs are adjacent (incl. diagonal) on grid
        pairs = cell_list_search(x3d, 3, 27, dx*1.8_fp, kernel, 8)
        call check_pairs(pairs, pairs3d_2, 158, "cell-lists (3d - adj pairs)")

    end subroutine test_cell_list

end module test_pair_finding

program run_tests

    use test_pair_finding, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests