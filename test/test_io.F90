module test_io

    use grasph_constants_m, only: fp, ndims
    use grasph_particle_system_m, only: particle_system_t
    use weakly_compressible_particles_m, only: eos_particles_t
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

        call psys%init(10, name)

        do i = 1, 10
            do d = 1, ndims
                psys%particles%x(d, i) = (i - 1)*ndims + d
                psys%particles%v(d, i) = ndims*10 + (i - 1)*ndims + d
                psys%particles%dvxdt(d, i) = ndims*10*2 + 30 + (i - 1)*ndims + d
            end do
            psys%particles%rho(i) = ndims*10*2 + i
            psys%particles%mass(i) = ndims*10*2 + 10 + i
            psys%particles%c(i) = ndims*10*2 + 20 + i
            psys%particles%drhodt(i) = ndims*10*3 + 30 + i
            psys%particles%id(i) = ndims*10*4 + 50 + i
            psys%particles%type(i) = ndims*10*4 + 60 + i
        end do

        call psys%register_io%register_variable("x", psys%particles%x)
        call psys%register_io%register_variable("v", psys%particles%v)
        call psys%register_io%register_variable("rho", psys%particles%rho)
        call psys%register_io%register_variable("mass", psys%particles%mass)
        call psys%register_io%register_variable("c", psys%particles%c)
        call psys%register_io%register_variable("dvxdt", psys%particles%dvxdt)
        call psys%register_io%register_variable("drhodt", psys%particles%drhodt)

        call psys%dump(1, "/tmp")

        call psys2%init(10, name)

        call psys2%register_io%register_variable("x", psys2%particles%x)
        call psys2%register_io%register_variable("v", psys2%particles%v)
        call psys2%register_io%register_variable("rho", psys2%particles%rho)
        call psys2%register_io%register_variable("mass", psys2%particles%mass)
        call psys2%register_io%register_variable("c", psys2%particles%c)
        call psys2%register_io%register_variable("dvxdt", psys2%particles%dvxdt)
        call psys2%register_io%register_variable("drhodt", psys2%particles%drhodt)

        call psys2%read("/tmp/grasph_particles_0000000001.h5", name)

        do i = 1, 10
            write (ic, "(I2)") i
            do d = 1, ndims
                write (dc, "(I1)") d
                call check(is_close(psys%particles%x(d, i), psys2%particles%x(d, i)), "Incorrect x for dim "//dc//", particle "//ic)
                call check(is_close(psys%particles%v(d, i), psys2%particles%v(d, i)), "Incorrect v for dim "//dc//", particle "//ic)
                call check(is_close(psys%particles%dvxdt(d, i), psys2%particles%dvxdt(d, i)), "Incorrect dvxdt for dim "//dc// &
                           ", particle "//ic)
            end do
            call check(is_equal(psys%particles%id(i), psys2%particles%id(i)), "Incorrect id for particle "//ic)
            call check(is_equal(psys%particles%type(i), psys2%particles%type(i)), "Incorrect type for particle "//ic)
            call check(is_close(psys%particles%rho(i), psys2%particles%rho(i)), "Incorrect rho for particle "//ic)
            call check(is_close(psys%particles%mass(i), psys2%particles%mass(i)), "Incorrect mass for particle "//ic)
            call check(is_close(psys%particles%c(i), psys2%particles%c(i)), "Incorrect c for particle "//ic)
            call check(is_close(psys%particles%drhodt(i), psys2%particles%drhodt(i)), "Incorrect drhodt for particle "//ic)
        end do

    end subroutine test_base_dump

    subroutine test_wcp_dump()

        type(particle_system_t), target:: psys, psys2
        integer:: i, d
        character(2):: ic
        character:: dc
        type(eos_particles_t), pointer:: ps_lhs, ps_rhs
        type(eos_particles_t):: ps_template
        character(*), parameter:: name = "test_wcp_particles"

        call psys%init(n=10, name=name, particle_template=ps_template)

        select type (psf => psys%particles)
        class is (eos_particles_t)
            ps_lhs => psf
            call psys%register_io%register_variable("x", psf%x)
            call psys%register_io%register_variable("v", psf%v)
            call psys%register_io%register_variable("rho", psf%rho)
            call psys%register_io%register_variable("mass", psf%mass)
            call psys%register_io%register_variable("c", psf%c)
            call psys%register_io%register_variable("dvxdt", psf%dvxdt)
            call psys%register_io%register_variable("drhodt", psf%drhodt)
            call psys%register_io%register_variable("p", psf%p)
        class default
            error stop "Expected eos_particles_t for psys%particles."
        end select

        do i = 1, 10
            do d = 1, ndims
                psys%particles%x(d, i) = (i - 1)*ndims + d
                psys%particles%v(d, i) = ndims*10 + (i - 1)*ndims + d
                psys%particles%dvxdt(d, i) = ndims*10*2 + 30 + (i - 1)*ndims + d
            end do
            psys%particles%rho(i) = ndims*10*2 + i
            psys%particles%mass(i) = ndims*10*2 + 10 + i
            psys%particles%c(i) = ndims*10*2 + 20 + i
            psys%particles%drhodt(i) = ndims*10*3 + 30 + i
            psys%particles%id(i) = ndims*10*4 + 50 + i
            psys%particles%type(i) = ndims*10*4 + 60 + i
            ps_lhs%p(i) = ndims*10*4 + 70 + i
        end do

        call psys%dump(1, "/tmp")

        call psys2%init(n=10, name=name, particle_template=ps_template)

        select type (psf => psys2%particles)
        class is (eos_particles_t)
            ps_rhs => psf
            call psys2%register_io%register_variable("x", psf%x)
            call psys2%register_io%register_variable("v", psf%v)
            call psys2%register_io%register_variable("rho", psf%rho)
            call psys2%register_io%register_variable("mass", psf%mass)
            call psys2%register_io%register_variable("c", psf%c)
            call psys2%register_io%register_variable("dvxdt", psf%dvxdt)
            call psys2%register_io%register_variable("drhodt", psf%drhodt)
            call psys2%register_io%register_variable("p", psf%p)
        class default
            error stop "Expected eos_particles_t for psys2%particles."
        end select

        call psys2%read("/tmp/grasph_particles_0000000001.h5", name)

        do i = 1, 10
            write (ic, "(I2)") i
            do d = 1, ndims
                write (dc, "(I1)") d
                call check(is_close(psys%particles%x(d, i), psys2%particles%x(d, i)), "Incorrect x for dim "//dc//", particle "//ic)
                call check(is_close(psys%particles%v(d, i), psys2%particles%v(d, i)), "Incorrect v for dim "//dc//", particle "//ic)
                call check(is_close(psys%particles%dvxdt(d, i), psys2%particles%dvxdt(d, i)), "Incorrect dvxdt for dim "//dc// &
                           ", particle "//ic)
            end do
            call check(is_equal(psys%particles%id(i), psys2%particles%id(i)), "Incorrect id for particle "//ic)
            call check(is_equal(psys%particles%type(i), psys2%particles%type(i)), "Incorrect type for particle "//ic)
            call check(is_close(psys%particles%rho(i), psys2%particles%rho(i)), "Incorrect rho for particle "//ic)
            call check(is_close(psys%particles%mass(i), psys2%particles%mass(i)), "Incorrect mass for particle "//ic)
            call check(is_close(psys%particles%c(i), psys2%particles%c(i)), "Incorrect c for particle "//ic)
            call check(is_close(psys%particles%drhodt(i), psys2%particles%drhodt(i)), "Incorrect drhodt for particle "//ic)
            call check(is_close(ps_lhs%p(i), ps_rhs%p(i)), "Incorrect p for particle "//ic)
        end do

    end subroutine test_wcp_dump

end module test_io

program run_tests

    use test_io, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
