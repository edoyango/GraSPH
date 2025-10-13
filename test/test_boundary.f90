module test_boundary

    use grasph_constants_m, only: fp
    use grasph_particle_system_m, only: particle_system_t
    use grasph_system_interactions_m, only: system_interaction_t
    use weakly_compressible_particles_m, only: eos_particle_t, eos_ghost_particle_t
    use weakly_compressible_interactions_m, only: ghost_timestep_setuper_t
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

contains

    type(test_list) function tests()

        tests = test_list([ &
                          test("test_vertical_walls", test_vertical_walls), &
                          test("test_diagonal_walls", test_diagonal_walls) &
                          ])

    end function tests

    subroutine test_vertical_walls()

        type(eos_particle_t):: ps_template
        type(eos_ghost_particle_t):: ps_ghost_template
        type(particle_system_t):: psys, psys_ghost
        type(ghost_timestep_setuper_t):: ghost_setuper
        type(system_interaction_t):: interaction
        integer:: i

        ! Expected    real       Expected
        ! LHS ghost  particles   RHS ghost
        ! particles              particles
        ! y=5      |         o | x
        !   4      |       o   |   x
        !   3      |     o     |
        !   2  x   |   o       |
        !   1    x | o         |
        !   x=-2-1 0 1 2 3 4 5 6 7 8

        call psys%init(5, "test_ghost_real", particle_template=ps_template)
        call psys_ghost%init(5, "test_ghost_ghost", particle_template=ps_ghost_template)

        ! generate real particles
        do i = 1, 5
            psys%particles(i)%id = i
            psys%particles(i)%x(:) = real(i, kind=fp)
        end do

        ! first test left wall
        ghost_setuper%cutoff = 2.5_fp
        ghost_setuper%surface_normal = [1._fp, 0._fp] ! normal vector pointing right-wards
        ghost_setuper%point = [0._fp, 0._fp]          ! plane/line travels through (0, 0)

        ! initialize system interaction with ghost setuper
        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 2), "Number of ghost particles generated for left wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), -1._fp), "1st left-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 1._fp), "1st left-wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), -2._fp), "2nd left-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 2._fp), "2nd left-wall ghost particle y-coordinate incorrect.")

        ! test right wall
        ! don't need to set cutoff as that's already set
        ghost_setuper%surface_normal = [-1._fp, 0._fp]
        ghost_setuper%point = [6._fp, 0._fp]

        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 2), "Number of ghost particles generated for right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), 8._fp), "1st right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 4._fp), "1st right-wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), 7._fp), "2nd right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 5._fp), "2nd right-wall ghost particle y-coordinate incorrect.")

        ! test particle penetration - real particles inside ghost region shouldn't be mirrored
        ! left wall
        ghost_setuper%cutoff = 2._fp ! tighten cutoff so only two particles are mirrored
        ghost_setuper%surface_normal = [1._fp, 0._fp]
        ghost_setuper%point = [1.5_fp, 0._fp] ! particle at (1,1) should be in ghost region

        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 2), "Number of ghost particles generated for right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), 1._fp), "1st right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 2._fp), "1st right-wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), 0._fp), "2nd right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 3._fp), "2nd right-wall ghost particle y-coordinate incorrect.")

        ! right wall
        ghost_setuper%surface_normal = [-1._fp, 0._fp]
        ghost_setuper%point = [4.5_fp, 0._fp] ! particle at (5,5) should be in ghost region

        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 2), "Number of ghost particles generated for right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), 6._fp), "1st right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 3._fp), "1st right-wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), 5._fp), "2nd right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 4._fp), "2nd right-wall ghost particle y-coordinate incorrect.")

    end subroutine test_vertical_walls

    subroutine test_diagonal_walls()

        type(eos_particle_t):: ps_template
        type(eos_ghost_particle_t):: ps_ghost_template
        type(particle_system_t):: psys, psys_ghost
        type(ghost_timestep_setuper_t):: ghost_setuper
        type(system_interaction_t):: interaction
        integer:: i
        character(1):: ic
        ! real particles generated along x = 1 y = [1, 5]
        ! config 1:         config 2:
        ! y = x - 1         y = 1 - x
        ! y = x + 5         y = 7 - x
        !       /              \
        !   x /                  \ x
        ! x / o                  o \ x
        ! /   o                  o   \
        !     o     /            o
        !     o   /          \   o
        !     o / x          x \ o
        !     / x              x \
        !                          \

        call psys%init(5, "test_ghost_real", particle_template=ps_template)
        call psys_ghost%init(5, "test_ghost_ghost", particle_template=ps_ghost_template)

        ! generate real particles along x = 1
        do i = 1, 5
            psys%particles(i)%id = i
            psys%particles(i)%x(1) = 1._fp
            psys%particles(i)%x(2) = real(i, kind=fp)
        end do

        ! test bottom-right wall (y = x - 1)
        ghost_setuper%cutoff = 2._fp
        ghost_setuper%surface_normal = [-sqrt(0.5_fp), sqrt(0.5_fp)] ! normal vector pointing top-left
        ghost_setuper%point = [0._fp, -1._fp]          ! plane/line travels through (0, 0)

        ! initialize system interaction with ghost setuper
        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 2), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), 2._fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 0._fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), 3._fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 0._fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle y-coordinate incorrect.")

        ! check particle in ghost region isn't mirrored
        ! move boundary to y = x + 0.5 to make particle at (1, 1) in ghost region
        ghost_setuper%point = [0._fp, 0.5_fp]
        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 3), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), 1.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 1.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), 2.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 1.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(3)%x(1), 3.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(3)%x(2), 1.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle y-coordinate incorrect.")

        ! test top-left wall (y = x + 5)
        ghost_setuper%surface_normal = [sqrt(0.5_fp), -sqrt(0.5_fp)]
        ghost_setuper%point = [0._fp, 5._fp]

        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 2), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), -1._fp, atol=1.d-14), &
                   "1st top-left wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 6._fp, atol=1.d-14), &
                   "1st top-left wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), 0._fp, atol=1.d-14), &
                   "2nd top-left wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 6._fp, atol=1.d-14), &
                   "2nd top-left wall ghost particle y-coordinate incorrect.")

        ! check particle in ghost region isn't mirrored
        ! move boundary to y = x + 3.5 to make particle at (5, 5) in ghost region
        ghost_setuper%point = [0._fp, 3.5_fp]
        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 3), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), -1.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 4.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), -0.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 4.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(3)%x(1), 0.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(3)%x(2), 4.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle y-coordinate incorrect.")

        ! test bottom-left wall (y = 1 - x)
        ghost_setuper%surface_normal = [sqrt(0.5_fp), sqrt(0.5_fp)] ! normal vector pointing top-right
        ghost_setuper%point = [0._fp, 1._fp]                        ! plane/line travels through (0, 5)

        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 2), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), 0._fp, atol=1.d-14), &
                   "1st bottom-left wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 0._fp, atol=1.d-14), &
                   "1st bottom-left wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), -1._fp, atol=1.d-14), &
                   "2nd bottom-left wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 0._fp, atol=1.d-14), &
                   "2nd bottom-left wall ghost particle y-coordinate incorrect.")

        ! check particle in ghost region isn't mirrored
        ! move boundary to y = 2.5 - x to make particle at (1, 1) in ghost region
        ghost_setuper%point = [0._fp, 2.5_fp]
        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 3), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), 0.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 1.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), -0.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 1.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(3)%x(1), -1.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(3)%x(2), 1.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle y-coordinate incorrect.")

        ! test top-right wall (y = 7 - x)
        ! ghost_setuper%cutoff = 1._fp
        ghost_setuper%surface_normal = [-sqrt(0.5_fp), -sqrt(0.5_fp)] ! normal vector pointing top-right
        ghost_setuper%point = [0._fp, 7._fp]                          ! plane/line travels through (0, 8)

        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 2), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), 3._fp, atol=1.d-14), &
                   "1st top-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 6._fp, atol=1.d-14), &
                   "1st top-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), 2._fp, atol=1.d-14), &
                   "2nd top-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 6._fp, atol=1.d-14), &
                   "2nd top-right wall ghost particle y-coordinate incorrect.")

        ! check particle in ghost region isn't mirrored
        ! move boundary to y = 5.5 - x to make particle at (5, 5) in ghost region
        ghost_setuper%point = [0._fp, 5.5_fp]
        call interaction%init(1, psys, psys_ghost, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size, 3), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(1), 3.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(1)%x(2), 4.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(1), 2.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(2)%x(2), 4.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(3)%x(1), 1.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles(3)%x(2), 4.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle y-coordinate incorrect.")

    end subroutine test_diagonal_walls

end module test_boundary

program run_tests

    use test_boundary, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
