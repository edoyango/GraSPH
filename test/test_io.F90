module test_io

    use grasph_constants, only: fp, ndims
    use grasph_particles, only: particle_system_t
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

        type(particle_system_t):: psys, psys2
        integer:: i, d
        character(2):: ic
        character:: dc
        character(*), parameter:: name = "test_base_particles"

        call psys%base_init(10, name)

        do i = 1, 10
            do d = 1, ndims
                psys%particles(i)%x(d) = (i - 1)*ndims + d
                psys%particles(i)%v(d) = ndims*10 + (i - 1)*ndims + d
                psys%particles(i)%dvxdt(d) = ndims*10*2 + 30 + (i - 1)*ndims + d
            end do
            psys%particles(i)%rho = ndims*10*2 + i
            psys%particles(i)%mass = ndims*10*2 + 10 + i
            psys%particles(i)%c = ndims*10*2 + 20 + i
            psys%particles(i)%drhodt = ndims*10*3 + 30 + i
            psys%particles(i)%id = ndims*10*4 + 50 + i
            psys%particles(i)%type = ndims*10*4 + 60 + i
        end do

        call psys%register_io%register_variable(psys%particles(1), "x", psys%particles(1)%x)
        call psys%register_io%register_variable(psys%particles(1), "v", psys%particles(1)%v)
        call psys%register_io%register_variable(psys%particles(1), "rho", psys%particles(1)%rho)
        call psys%register_io%register_variable(psys%particles(1), "mass", psys%particles(1)%mass)
        call psys%register_io%register_variable(psys%particles(1), "c", psys%particles(1)%c)
        call psys%register_io%register_variable(psys%particles(1), "dvxdt", psys%particles(1)%dvxdt)
        call psys%register_io%register_variable(psys%particles(1), "drhodt", psys%particles(1)%drhodt)

        call psys%dump(1, "/tmp")

        call psys2%base_init(10, name)

        call psys2%register_io%register_variable(psys2%particles(1), "x", psys2%particles(1)%x)
        call psys2%register_io%register_variable(psys2%particles(1), "v", psys2%particles(1)%v)
        call psys2%register_io%register_variable(psys2%particles(1), "rho", psys2%particles(1)%rho)
        call psys2%register_io%register_variable(psys2%particles(1), "mass", psys2%particles(1)%mass)
        call psys2%register_io%register_variable(psys2%particles(1), "c", psys2%particles(1)%c)
        call psys2%register_io%register_variable(psys2%particles(1), "dvxdt", psys2%particles(1)%dvxdt)
        call psys2%register_io%register_variable(psys2%particles(1), "drhodt", psys2%particles(1)%drhodt)

        call psys2%read("/tmp/grasph_particles_0000000001.h5", name)

        do i = 1, 10
            write (ic, "(I2)") i
            do d = 1, ndims
                write (dc, "(I1)") d
                call check(is_close(psys%particles(i)%x(d), psys2%particles(i)%x(d)), "Incorrect x for dim "//dc//", particle "//ic)
                call check(is_close(psys%particles(i)%v(d), psys2%particles(i)%v(d)), "Incorrect v for dim "//dc//", particle "//ic)
                call check(is_close(psys%particles(i)%dvxdt(d), psys2%particles(i)%dvxdt(d)), "Incorrect dvxdt for dim "//dc// &
                           ", particle "//ic)
            end do
            call check(is_equal(psys%particles(i)%id, psys2%particles(i)%id), "Incorrect id for particle "//ic)
            call check(is_equal(psys%particles(i)%type, psys2%particles(i)%type), "Incorrect type for particle "//ic)
            call check(is_close(psys%particles(i)%rho, psys2%particles(i)%rho), "Incorrect rho for particle "//ic)
            call check(is_close(psys%particles(i)%mass, psys2%particles(i)%mass), "Incorrect mass for particle "//ic)
            call check(is_close(psys%particles(i)%c, psys2%particles(i)%c), "Incorrect c for particle "//ic)
            call check(is_close(psys%particles(i)%drhodt, psys2%particles(i)%drhodt), "Incorrect drhodt for particle "//ic)
        end do

    end subroutine test_base_dump

    subroutine test_wcp_dump()

        type(particle_system_t):: psys, psys2
        integer:: i, d
        character(2):: ic
        character:: dc
        type(eos_particle), pointer:: ps_lhs(:), ps_rhs(:)
        type(eos_particle):: ps_template
        character(*), parameter:: name = "test_wcp_particles"

        call psys%base_init(n=10, name=name, particle_template=ps_template)

        select type (psf => psys%particles)
        class is (eos_particle)
            ps_lhs => psf
            call psys%register_io%register_variable(psf(1), "x", psf(1)%x)
            call psys%register_io%register_variable(psf(1), "v", psf(1)%v)
            call psys%register_io%register_variable(psf(1), "rho", psf(1)%rho)
            call psys%register_io%register_variable(psf(1), "mass", psf(1)%mass)
            call psys%register_io%register_variable(psf(1), "c", psf(1)%c)
            call psys%register_io%register_variable(psf(1), "dvxdt", psf(1)%dvxdt)
            call psys%register_io%register_variable(psf(1), "drhodt", psf(1)%drhodt)
            call psys%register_io%register_variable(psf(1), "p", psf(1)%p)
        class default
            error stop "Expected eos_particle for psys%particles."
        end select

        do i = 1, 10
            do d = 1, ndims
                psys%particles(i)%x(d) = (i - 1)*ndims + d
                psys%particles(i)%v(d) = ndims*10 + (i - 1)*ndims + d
                psys%particles(i)%dvxdt(d) = ndims*10*2 + 30 + (i - 1)*ndims + d
            end do
            psys%particles(i)%rho = ndims*10*2 + i
            psys%particles(i)%mass = ndims*10*2 + 10 + i
            psys%particles(i)%c = ndims*10*2 + 20 + i
            psys%particles(i)%drhodt = ndims*10*3 + 30 + i
            psys%particles(i)%id = ndims*10*4 + 50 + i
            psys%particles(i)%type = ndims*10*4 + 60 + i
            ps_lhs(i)%p = ndims*10*4 + 70 + i
        end do

        call psys%dump(1, "/tmp")

        call psys2%base_init(n=10, name=name, particle_template=ps_template)

        select type (psf => psys2%particles)
        class is (eos_particle)
            ps_rhs => psf
            call psys2%register_io%register_variable(psf(1), "x", psf(1)%x)
            call psys2%register_io%register_variable(psf(1), "v", psf(1)%v)
            call psys2%register_io%register_variable(psf(1), "rho", psf(1)%rho)
            call psys2%register_io%register_variable(psf(1), "mass", psf(1)%mass)
            call psys2%register_io%register_variable(psf(1), "c", psf(1)%c)
            call psys2%register_io%register_variable(psf(1), "dvxdt", psf(1)%dvxdt)
            call psys2%register_io%register_variable(psf(1), "drhodt", psf(1)%drhodt)
            call psys2%register_io%register_variable(psf(1), "p", psf(1)%p)
        class default
            error stop "Expected eos_particle for psys2%particles."
        end select

        call psys2%read("/tmp/grasph_particles_0000000001.h5", name)

        do i = 1, 10
            write (ic, "(I2)") i
            do d = 1, ndims
                write (dc, "(I1)") d
                call check(is_close(psys%particles(i)%x(d), psys2%particles(i)%x(d)), "Incorrect x for dim "//dc//", particle "//ic)
                call check(is_close(psys%particles(i)%v(d), psys2%particles(i)%v(d)), "Incorrect v for dim "//dc//", particle "//ic)
                call check(is_close(psys%particles(i)%dvxdt(d), psys2%particles(i)%dvxdt(d)), "Incorrect dvxdt for dim "//dc// &
                           ", particle "//ic)
            end do
            call check(is_equal(psys%particles(i)%id, psys2%particles(i)%id), "Incorrect id for particle "//ic)
            call check(is_equal(psys%particles(i)%type, psys2%particles(i)%type), "Incorrect type for particle "//ic)
            call check(is_close(psys%particles(i)%rho, psys2%particles(i)%rho), "Incorrect rho for particle "//ic)
            call check(is_close(psys%particles(i)%mass, psys2%particles(i)%mass), "Incorrect mass for particle "//ic)
            call check(is_close(psys%particles(i)%c, psys2%particles(i)%c), "Incorrect c for particle "//ic)
            call check(is_close(psys%particles(i)%drhodt, psys2%particles(i)%drhodt), "Incorrect drhodt for particle "//ic)
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
