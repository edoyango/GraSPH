module test_io

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles, wcp => weakly_compressible_particles
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

        call ps%base_init(10, 3, "test_base_particles")

        ps%x(:, :) = reshape([(real(i, kind=fp), i=1, 30)], shape(ps%x))
        ps%v(:, :) = reshape([(real(i, kind=fp), i=31, 60)], shape(ps%v))
        ps%rho(:) = [(real(i, kind=fp), i=61, 70)]
        ps%mass(:) = [(real(i, kind=fp), i=71, 80)]
        ps%c(:) = [(real(i, kind=fp), i=81, 90)]
        ps%dvxdt(:, :) = reshape([(real(i, kind=fp), i=91, 120)], shape(ps%dvxdt))
        ps%drhodt(:) = [(real(i, kind=fp), i=121, 130)]
        ps%v0(:, :) = reshape([(real(i, kind=fp), i=131, 160)], shape(ps%v0))
        ps%rho0(:) = [(real(i, kind=fp), i=161, 170)]
        ps%id(:) = [(i, i=171, 180)]
        ps%type(:) = [(i, i=181, 190)]

        call ps%dump(1, "/tmp", "test-", 4)

        call ps2%read("/tmp/test-grasph_particles_0000000001.h5", "test_base_particles")

        do i = 1, 10
            write(ic, "(I2)") i
            do d = 1, 3
                write(dc, "(I1)") d
                call check(is_close(ps%x(d, i), ps2%x(d, i)), "Incorrect x for dim " // dc // ", particle " // ic)
                call check(is_close(ps%v(d, i), ps2%v(d, i)), "Incorrect v for dim " // dc // ", particle " // ic)
                call check(is_close(ps%v0(d, i), ps2%v0(d, i)), "Incorrect v0 for dim " // dc // ", particle " // ic)
                call check(is_close(ps%dvxdt(d, i), ps2%dvxdt(d, i)), "Incorrect dvxdt for dim " // dc // ", particle " // ic)
            enddo
            call check(is_equal(ps%id(i), ps2%id(i)), "Incorrect id for particle " // ic)
            call check(is_equal(ps%type(i), ps2%type(i)), "Incorrect type for particle " // ic)
            call check(is_close(ps%rho(i), ps2%rho(i)), "Incorrect rho for particle " // ic)
            call check(is_close(ps%rho0(i), ps2%rho0(i)), "Incorrect rho0 for particle " // ic)
            call check(is_close(ps%mass(i), ps2%mass(i)), "Incorrect mass for particle " // ic)
            call check(is_close(ps%c(i), ps2%c(i)), "Incorrect c for particle " // ic)
            call check(is_close(ps%drhodt(i), ps2%drhodt(i)), "Incorrect drhodt for particle " // ic)
        enddo

    end subroutine test_base_dump

    subroutine test_wcp_dump()

        type(wcp):: ps, ps2
        integer:: i, d
        character(2):: ic
        character:: dc

        call ps%init(10, 3, "test_wcp_particles", 1000._fp)

        ps%x(:, :) = reshape([(real(i, kind=fp), i=1, 30)], shape(ps%x))
        ps%v(:, :) = reshape([(real(i, kind=fp), i=31, 60)], shape(ps%v))
        ps%rho(:) = [(real(i, kind=fp), i=61, 70)]
        ps%mass(:) = [(real(i, kind=fp), i=71, 80)]
        ps%c(:) = [(real(i, kind=fp), i=81, 90)]
        ps%dvxdt(:, :) = reshape([(real(i, kind=fp), i=91, 120)], shape(ps%dvxdt))
        ps%drhodt(:) = [(real(i, kind=fp), i=121, 130)]
        ps%v0(:, :) = reshape([(real(i, kind=fp), i=131, 160)], shape(ps%v0))
        ps%rho0(:) = [(real(i, kind=fp), i=161, 170)]
        ps%id(:) = [(i, i=171, 180)]
        ps%type(:) = [(i, i=181, 190)]
        ps%p(:) = [(i, i=191, 200)]

        call ps%dump(1, "/tmp", "test-", 4)

        call ps2%read("/tmp/test-grasph_particles_0000000001.h5", "test_wcp_particles")

        do i = 1, 10
            write(ic, "(I2)") i
            do d = 1, 3
                write(dc, "(I1)") d
                call check(is_close(ps%x(d, i), ps2%x(d, i)), "Incorrect x for dim " // dc // ", particle " // ic)
                call check(is_close(ps%v(d, i), ps2%v(d, i)), "Incorrect v for dim " // dc // ", particle " // ic)
                call check(is_close(ps%v0(d, i), ps2%v0(d, i)), "Incorrect v0 for dim " // dc // ", particle " // ic)
                call check(is_close(ps%dvxdt(d, i), ps2%dvxdt(d, i)), "Incorrect dvxdt for dim " // dc // ", particle " // ic)
            enddo
            call check(is_equal(ps%id(i), ps2%id(i)), "Incorrect id for particle " // ic)
            call check(is_equal(ps%type(i), ps2%type(i)), "Incorrect type for particle " // ic)
            call check(is_close(ps%rho(i), ps2%rho(i)), "Incorrect rho for particle " // ic)
            call check(is_close(ps%rho0(i), ps2%rho0(i)), "Incorrect rho0 for particle " // ic)
            call check(is_close(ps%mass(i), ps2%mass(i)), "Incorrect mass for particle " // ic)
            call check(is_close(ps%c(i), ps2%c(i)), "Incorrect c for particle " // ic)
            call check(is_close(ps%drhodt(i), ps2%drhodt(i)), "Incorrect drhodt for particle " // ic)
            call check(is_close(ps%p(i), ps2%p(i)), "Incorrect p for particle " // ic)
        enddo

    end subroutine test_wcp_dump

end module test_io

program run_tests

    use test_io, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests