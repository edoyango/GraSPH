module test_kernel_values

    use grasph_constants, only: fp, pi
    use grasph_kernels_m, only: base_kernel_t, cubic_bspline_kernel_t
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

contains

    type(test_list) function tests()

        tests = test_list([ &
                          test("test_cubic_spline_values", test_cubic_spline_values) &
                          ])

    end function tests

    subroutine test_cubic_spline_values()
        implicit none
        type(cubic_bspline_kernel_t):: my_kernel
        real(fp):: w, dwdx(3)
        integer:: d
        real(fp), parameter:: h = 1.2_fp, alpha2d = 10._fp/(7._fp*pi*h*h), alpha3d = 1._fp/(pi*h*h*h)

        ! start off with 2D tests
        ! test correct initialization
        call my_kernel%init(2, 1.2_fp)

        call check( &
            is_equal(my_kernel%d, 2), &
            "Cubic spline kernel dimension not correctly initialized!" &
            )
        call check( &
            is_close(my_kernel%h, h), &
            "Cubic spline kernel smoothing length not correctly initialized!" &
            )
        call check( &
            is_close(my_kernel%alpha, alpha2d), &
            "Cubic spline kernel normalization factor not correctly initialized!" &
            )

        ! check values when dx = [0, 0]
        ! not checking dwdx at 0, as we assume that never happens
        call my_kernel%values([0._fp, 0._fp], w, dwdx)

        call check( &
            is_close(w, alpha2d*(0.25_fp*8._fp - 1._fp)), &
            "Cubic spline kernel incorrect at 0!" &
            )

        ! check values are 0 at kernel perimeter
        call my_kernel%values([2*h, 2*h], w, dwdx)

        call check( &
            is_close(w, 0.d0), &
            "Cubic spline kernel incorrect at 2h!" &
            )
        do d = 1, 2
            call check( &
                is_close(dwdx(d), 0.d0), &
                "Cubic spline kernel gradient incorrect at 2h!" &
                )
        end do

        ! check values are 0 beyond kernel perimeter
        call my_kernel%values([3*h, 3*h], w, dwdx)

        call check( &
            is_close(w, 0.d0), &
            "Cubic spline kernel incorrect at 3h!" &
            )
        do d = 1, 2
            call check( &
                is_close(dwdx(d), 0.d0), &
                "Cubic spline kernel gradient incorrect at 3h!" &
                )
        end do

        ! check values at 1h
        call my_kernel%values([h, 0._fp], w, dwdx)

        call check( &
            is_close(w, alpha2d*0.25_fp), &
            "Cubic spline kernel incorrect at 1h!" &
            )
        call check( &
            is_close(dwdx(1), -alpha2d*0.75_fp/h), &
            "Cubic spline kernel x-gradient incorrect at 1h!" &
            )
        call check( &
            is_close(dwdx(2), 0._fp), &
            "Cubic spline kernel y-gradient incorrect at 1h!" &
            )

        call my_kernel%values([0._fp, h], w, dwdx)

        call check( &
            is_close(w, alpha2d*0.25_fp), &
            "Cubic spline kernel incorrect at 1h!" &
            )
        call check( &
            is_close(dwdx(1), 0._fp), &
            "Cubic spline kernel x-gradient incorrect at 1h!" &
            )
        call check( &
            is_close(dwdx(2), -alpha2d*0.75_fp/h), &
            "Cubic spline kernel y-gradient incorrect at 1h!" &
            )

        call my_kernel%values([-h/sqrt(2._fp), -h/sqrt(2._fp)], w, dwdx)

        call check( &
            is_close(w, alpha2d*0.25_fp), &
            "Cubic spline kernel incorrect at 1h!" &
            )
        call check( &
            is_close(dwdx(1), alpha2d*0.75_fp/(h*sqrt(2._fp))), &
            "Cubic spline kernel x-gradient incorrect at 1h!" &
            )
        call check( &
            is_close(dwdx(2), alpha2d*0.75_fp/(h*sqrt(2._fp))), &
            "Cubic spline kernel y-gradient incorrect at 1h!" &
            )

        ! do 3d tests
        ! test correct initialization
        call my_kernel%init(3, 1.2_fp)

        call check( &
            is_equal(my_kernel%d, 3), &
            "Cubic spline kernel dimension not correctly initialized!" &
            )
        call check( &
            is_close(my_kernel%h, h), &
            "Cubic spline kernel smoothing length not correctly initialized!" &
            )
        call check( &
            is_close(my_kernel%alpha, alpha3d), &
            "Cubic spline kernel normalization factor not correctly initialized!" &
            )

        ! do single check to confirm correct alpha and update of 3d dwdx
        call my_kernel%values([-h/sqrt(3._fp), -h/sqrt(3._fp), -h/sqrt(3._fp)], w, dwdx)

        call check( &
            is_close(w, alpha3d*0.25_fp), &
            "Cubic spline kernel incorrect at 1h!" &
            )
        do d = 1, 3
            call check( &
                is_close(dwdx(d), alpha3d*0.75_fp/(h*sqrt(3._fp))), &
                "Cubic spline kernel gradient incorrect at 1h!" &
                )
        end do

    end subroutine test_cubic_spline_values

end module test_kernel_values

program run_tests

    use test_kernel_values, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
