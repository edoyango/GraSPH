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

    subroutine check_pairs(pairs, x, correct_pairs, ncorrect_pairs, case_string, kernel)

        type(particle_pairs), intent(in):: pairs
        real(fp), intent(in):: x(:, :)
        integer, intent(in):: ncorrect_pairs, correct_pairs(2, ncorrect_pairs)
        character(*), intent(in):: case_string
        type(grasph_cubic_bspline_kernel), intent(in):: kernel
        integer:: n, i, j, jj, ii
        character(3):: c, ic, jc
        real(fp):: w, dwdx(kernel%d)

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
                write(ic, "(I3)") i
                write(jc, "(I3)") j
                do ii = 1, ncorrect_pairs
                    if ((correct_pairs(1, ii) == i .and. correct_pairs(2, ii) == j) .or. &
                        (correct_pairs(1, ii) == j .and. correct_pairs(2, ii) == i)) then
                        ! just to make sure the check index is updated properly
                        call check(.true., "")
                        ! checking the kernel values in pairs have been updated properly
                        call kernel%values(x(:, i) - x(:, j), w, dwdx)
                        call check( &
                            is_close(w, pairs%w(jj)), &
                            "Case: " // case_string // ". Incorrectly updated w value: " // ic // ", j :" // jc &
                        )
                        cycle i_rhs_sweep
                    endif
                enddo
                call check(.false., "Case: " // case_string // ". Couldn't find pair - i: " // ic // ", j: " // jc)
            enddo i_rhs_sweep
        enddo

    end subroutine check_pairs

    subroutine test_dsearch()

        type(particle_pairs):: pairs
        type(grasph_cubic_bspline_kernel):: kernel

        call kernel%init(2, 1._fp)

        ! 2d - cutoff selected for all particles to be paired with eachother
        pairs = dsearch(x2d, 2, 16, sqrt(2._fp), kernel, 15)
        call check_pairs(pairs, x2d, pairs2d_1, 120, "brute-force (2d - all pairs)", kernel)

        ! 2d -cutoff selected so pairs are adjacent (incl. diagonal) on grid
        pairs = dsearch(x2d, 2, 16, dx*1.8_fp, kernel, 8)
        call check_pairs(pairs, x2d, pairs2d_2, 42, "brute-force (2d - adj pairs)", kernel)

        call kernel%init(3, 1._fp)

        ! 3d - cutoff selected for all particles to be paired with eachother
        pairs = dsearch(x3d, 3, 27, sqrt(3._fp), kernel, 26)
        call check_pairs(pairs, x3d, pairs3d_1, 351, "brute-force (3d - all pairs)", kernel)

        ! 3d - cutoff selected so pairs are adjacent (incl. diagonal) on grid
        pairs = dsearch(x3d, 3, 27, dx*1.8_fp, kernel, 8)
        call check_pairs(pairs, x3d, pairs3d_2, 158, "brute-force (3d - adj pairs)", kernel)

    end subroutine test_dsearch
    
    subroutine test_cell_list()

        type(particle_pairs):: pairs
        type(grasph_cubic_bspline_kernel):: kernel

        call kernel%init(2, 1._fp)

        ! 2d - cutoff selected for all particles to be paired with eachother
        pairs = cell_list_search(x2d, 2, 16, sqrt(2._fp), kernel, 16)
        call check_pairs(pairs, x2d, pairs2d_1, 120, "cell-lists (2d - all pairs)", kernel)

        ! 2d -cutoff selected so pairs are adjacent (incl. diagonal) on grid
        pairs = cell_list_search(x2d, 2, 16, dx*1.8_fp, kernel, 8)
        call check_pairs(pairs, x2d, pairs2d_2, 42, "cell-lists (3d - adj pairs)", kernel)

        call kernel%init(3, 1._fp)

        ! 3d - cutoff selected for all particles to be paired with eachother
        pairs = cell_list_search(x3d, 3, 27, sqrt(3._fp), kernel, 27)
        call check_pairs(pairs, x3d, pairs3d_1, 351, "cell-lists (3d - all pairs)", kernel)

        ! 3d - cutoff selected so pairs are adjacent (incl. diagonal) on grid
        pairs = cell_list_search(x3d, 3, 27, dx*1.8_fp, kernel, 8)
        call check_pairs(pairs, x3d, pairs3d_2, 158, "cell-lists (3d - adj pairs)", kernel)

    end subroutine test_cell_list

end module test_pair_finding

program run_tests

    use test_pair_finding, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests