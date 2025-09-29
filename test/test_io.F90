module test_io

    use grasph_constants, only: fp, ndims
    use grasph_particles, only: base_particles
    use weakly_compressible_particles, only: wcp => linear_eos_particles, linear_eos_particle
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

contains

    type(test_list) function tests()

        tests = test_list([ &
                          test("test_base_dump", test_base_dump), &
                          test("test_wcp_dump", test_wcp_dump) &
                          ])

    end function tests

    subroutine test_base_dump()

        type(base_particles):: ps, ps2
        integer:: i, d
        character(2):: ic
        character:: dc

        call ps%base_init(10, "test_base_particles")

        do i = 1, 10
            do d = 1, ndims
                ps%ps(i)%x(d) = (i - 1)*ndims + d
                ps%ps(i)%v(d) = ndims*10 + (i - 1)*ndims + d
                ps%ps(i)%dvxdt(d) = ndims*10*2 + 30 + (i - 1)*ndims + d
            end do
            ps%ps(i)%rho = ndims*10*2 + i
            ps%ps(i)%mass = ndims*10*2 + 10 + i
            ps%ps(i)%c = ndims*10*2 + 20 + i
            ps%ps(i)%drhodt = ndims*10*3 + 30 + i
            ps%ps(i)%id = ndims*10*4 + 50 + i
            ps%ps(i)%type = ndims*10*4 + 60 + i
        end do

        call ps%dump(1, "/tmp")

        call ps2%read("/tmp/grasph_particles_0000000001.h5", "test_base_particles")

        do i = 1, 10
            write (ic, "(I2)") i
            do d = 1, 3
                write (dc, "(I1)") d
                call check(is_close(ps%ps(i)%x(d), ps2%ps(i)%x(d)), "Incorrect x for dim "//dc//", particle "//ic)
                call check(is_close(ps%ps(i)%v(d), ps2%ps(i)%v(d)), "Incorrect v for dim "//dc//", particle "//ic)
                call check(is_close(ps%ps(i)%dvxdt(d), ps2%ps(i)%dvxdt(d)), "Incorrect dvxdt for dim "//dc//", particle "//ic)
            end do
            call check(is_equal(ps%ps(i)%id, ps2%ps(i)%id), "Incorrect id for particle "//ic)
            call check(is_equal(ps%ps(i)%type, ps2%ps(i)%type), "Incorrect type for particle "//ic)
            call check(is_close(ps%ps(i)%rho, ps2%ps(i)%rho), "Incorrect rho for particle "//ic)
            call check(is_close(ps%ps(i)%mass, ps2%ps(i)%mass), "Incorrect mass for particle "//ic)
            call check(is_close(ps%ps(i)%c, ps2%ps(i)%c), "Incorrect c for particle "//ic)
            call check(is_close(ps%ps(i)%drhodt, ps2%ps(i)%drhodt), "Incorrect drhodt for particle "//ic)
        end do

    end subroutine test_base_dump

    subroutine test_wcp_dump()

        type(wcp):: ps, ps2
        integer:: i, d
        character(2):: ic
        character:: dc
        type(linear_eos_particle), pointer:: ps_lhs(:), ps_rhs(:)

        call ps%init(10, "test_wcp_particles", 1000._fp)

        select type (psf => ps%ps)
        class is (linear_eos_particle)
            ps_lhs => psf
        class default
            error stop "Expected linear_eos_particle for ps%ps."
        end select

        do i = 1, 10
            do d = 1, ndims
                ps%ps(i)%x(d) = (i - 1)*ndims + d
                ps%ps(i)%v(d) = ndims*10 + (i - 1)*ndims + d
                ps%ps(i)%dvxdt(d) = ndims*10*2 + 30 + (i - 1)*ndims + d
            end do
            ps%ps(i)%rho = ndims*10*2 + i
            ps%ps(i)%mass = ndims*10*2 + 10 + i
            ps%ps(i)%c = ndims*10*2 + 20 + i
            ps%ps(i)%drhodt = ndims*10*3 + 30 + i
            ps%ps(i)%id = ndims*10*4 + 50 + i
            ps%ps(i)%type = ndims*10*4 + 60 + i
            ps_lhs(i)%p = ndims*10*4 + 70 + i
        end do

        call ps%dump(1, "/tmp", "test-grasph_particles", 4)

        call ps2%read("/tmp/test-grasph_particles_0000000001.h5", "test_wcp_particles")

        select type (psf => ps2%ps)
        class is (linear_eos_particle)
            ps_rhs => psf
        class default
            error stop "Expected linear_eos_particle for ps2%ps."
        end select

        do i = 1, 10
            write (ic, "(I2)") i
            do d = 1, 3
                write (dc, "(I1)") d
                call check(is_close(ps%ps(i)%x(d), ps2%ps(i)%x(d)), "Incorrect x for dim "//dc//", particle "//ic)
                call check(is_close(ps%ps(i)%v(d), ps2%ps(i)%v(d)), "Incorrect v for dim "//dc//", particle "//ic)
                call check(is_close(ps%ps(i)%dvxdt(d), ps2%ps(i)%dvxdt(d)), "Incorrect dvxdt for dim "//dc//", particle "//ic)
            end do
            call check(is_equal(ps%ps(i)%id, ps2%ps(i)%id), "Incorrect id for particle "//ic)
            call check(is_equal(ps%ps(i)%type, ps2%ps(i)%type), "Incorrect type for particle "//ic)
            call check(is_close(ps%ps(i)%rho, ps2%ps(i)%rho), "Incorrect rho for particle "//ic)
            call check(is_close(ps%ps(i)%mass, ps2%ps(i)%mass), "Incorrect mass for particle "//ic)
            call check(is_close(ps%ps(i)%c, ps2%ps(i)%c), "Incorrect c for particle "//ic)
            call check(is_close(ps%ps(i)%drhodt, ps2%ps(i)%drhodt), "Incorrect drhodt for particle "//ic)
            call check(is_close(ps_lhs(i)%p, ps_rhs(i)%p), "Incorrect p for particle "//ic)
        end do

    end subroutine test_wcp_dump

end module test_io

program run_tests

    use test_io, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
