module test_io

    use grasph_constants, only: fp, ndims
    use grasph_particles, only: base_particles
    use weakly_compressible_particles, only: eos_particle
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
        character(*), parameter:: name = "test_base_particles"

        call ps%base_init(10, name)

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

        call ps%register_io%register_variable(ps%ps(1), "x", ps%ps(1)%x)
        call ps%register_io%register_variable(ps%ps(1), "v", ps%ps(1)%v)
        call ps%register_io%register_variable(ps%ps(1), "rho", ps%ps(1)%rho)
        call ps%register_io%register_variable(ps%ps(1), "mass", ps%ps(1)%mass)
        call ps%register_io%register_variable(ps%ps(1), "c", ps%ps(1)%c)
        call ps%register_io%register_variable(ps%ps(1), "dvxdt", ps%ps(1)%dvxdt)
        call ps%register_io%register_variable(ps%ps(1), "drhodt", ps%ps(1)%drhodt)

        call ps%dump(1, "/tmp")

        call ps2%base_init(10, name)

        call ps2%register_io%register_variable(ps2%ps(1), "x", ps2%ps(1)%x)
        call ps2%register_io%register_variable(ps2%ps(1), "v", ps2%ps(1)%v)
        call ps2%register_io%register_variable(ps2%ps(1), "rho", ps2%ps(1)%rho)
        call ps2%register_io%register_variable(ps2%ps(1), "mass", ps2%ps(1)%mass)
        call ps2%register_io%register_variable(ps2%ps(1), "c", ps2%ps(1)%c)
        call ps2%register_io%register_variable(ps2%ps(1), "dvxdt", ps2%ps(1)%dvxdt)
        call ps2%register_io%register_variable(ps2%ps(1), "drhodt", ps2%ps(1)%drhodt)

        call ps2%read("/tmp/grasph_particles_0000000001.h5", name)

        do i = 1, 10
            write (ic, "(I2)") i
            do d = 1, ndims
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

        type(base_particles):: ps, ps2
        integer:: i, d
        character(2):: ic
        character:: dc
        type(eos_particle), pointer:: ps_lhs(:), ps_rhs(:)
        type(eos_particle):: ps_template
        character(*), parameter:: name = "test_wcp_particles"

        call ps%base_init(n=10, name=name, ps_template=ps_template)

        select type (psf => ps%ps)
        class is (eos_particle)
            ps_lhs => psf
            call ps%register_io%register_variable(psf(1), "x", psf(1)%x)
            call ps%register_io%register_variable(psf(1), "v", psf(1)%v)
            call ps%register_io%register_variable(psf(1), "rho", psf(1)%rho)
            call ps%register_io%register_variable(psf(1), "mass", psf(1)%mass)
            call ps%register_io%register_variable(psf(1), "c", psf(1)%c)
            call ps%register_io%register_variable(psf(1), "dvxdt", psf(1)%dvxdt)
            call ps%register_io%register_variable(psf(1), "drhodt", psf(1)%drhodt)
            call ps%register_io%register_variable(psf(1), "p", psf(1)%p)
        class default
            error stop "Expected eos_particle for ps%ps."
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

        call ps%dump(1, "/tmp")

        call ps2%base_init(n=10, name=name, ps_template=ps_template)

        select type (psf => ps2%ps)
        class is (eos_particle)
            ps_rhs => psf
            call ps2%register_io%register_variable(psf(1), "x", psf(1)%x)
            call ps2%register_io%register_variable(psf(1), "v", psf(1)%v)
            call ps2%register_io%register_variable(psf(1), "rho", psf(1)%rho)
            call ps2%register_io%register_variable(psf(1), "mass", psf(1)%mass)
            call ps2%register_io%register_variable(psf(1), "c", psf(1)%c)
            call ps2%register_io%register_variable(psf(1), "dvxdt", psf(1)%dvxdt)
            call ps2%register_io%register_variable(psf(1), "drhodt", psf(1)%drhodt)
            call ps2%register_io%register_variable(psf(1), "p", psf(1)%p)
        class default
            error stop "Expected eos_particle for ps2%ps."
        end select

        call ps2%read("/tmp/grasph_particles_0000000001.h5", name)

        do i = 1, 10
            write (ic, "(I2)") i
            do d = 1, ndims
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
