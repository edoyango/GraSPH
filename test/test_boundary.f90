module test_boundary

    use grasph_constants_m, only: fp, ndims
    use grasph_particle_m, only: base_particles_t
    use grasph_particle_system_m, only: particle_system_t
    use grasph_system_interactions_m, only: system_interaction_t, sweeper_container_t
    use weakly_compressible_particles_m, only: eos_particles_t, eos_ghost_particles_t
    use grasph_kernels_m, only: cubic_bspline_kernel_t
    use weakly_compressible_interactions_m, only: ghost_timestep_setuper_t, morris_boundary_sweeper_t
    use fortuno_serial, only: is_equal, is_close, test => serial_case_item, check => serial_check, test_list

    implicit none

    private
    public:: tests

contains

    type(test_list) function tests()

        tests = test_list([ &
                          test("test_coaxial_ghost_walls", test_coaxial_ghost_walls), &
                          test("test_diagonal_ghost_walls", test_diagonal_ghost_walls), &
                          test("test_morris_walls", test_morris_walls) &
                          ])

    end function tests

    subroutine test_coaxial_ghost_walls()

        type(eos_particles_t):: ps_template
        type(eos_ghost_particles_t):: ps_ghost_template
        type(particle_system_t):: psys, psys_ghost
        type(ghost_timestep_setuper_t):: ghost_setuper
        type(system_interaction_t):: interaction
        type(sweeper_container_t):: sweepers(0)
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
            psys%particles%id(i) = i
            psys%particles%x(:, i) = real(i, kind=fp)
        end do

        ! first test left wall
        ghost_setuper%cutoff = 2.5_fp
        ghost_setuper%surface_normal = [1._fp, 0._fp] ! normal vector pointing right-wards
        ghost_setuper%point = [0._fp, 0._fp]          ! plane/line travels through (0, 0)

        ! initialize system interaction with ghost setuper
        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 2), "Number of ghost particles generated for left wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), -1._fp), "1st left-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 1._fp), "1st left-wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), -2._fp), "2nd left-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 2._fp), "2nd left-wall ghost particle y-coordinate incorrect.")

        ! test right wall
        ! don't need to set cutoff as that's already set
        ghost_setuper%surface_normal = [-1._fp, 0._fp]
        ghost_setuper%point = [6._fp, 0._fp]

        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 2), "Number of ghost particles generated for right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), 8._fp), "1st right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 4._fp), "1st right-wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), 7._fp), "2nd right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 5._fp), "2nd right-wall ghost particle y-coordinate incorrect.")

        ! test particle penetration - real particles inside ghost region shouldn't be mirrored
        ! left wall
        ghost_setuper%cutoff = 2._fp ! tighten cutoff so only two particles are mirrored
        ghost_setuper%surface_normal = [1._fp, 0._fp]
        ghost_setuper%point = [1.5_fp, 0._fp] ! particle at (1,1) should be in ghost region

        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 2), "Number of ghost particles generated for right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), 1._fp), "1st right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 2._fp), "1st right-wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), 0._fp), "2nd right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 3._fp), "2nd right-wall ghost particle y-coordinate incorrect.")

        ! right wall
        ghost_setuper%surface_normal = [-1._fp, 0._fp]
        ghost_setuper%point = [4.5_fp, 0._fp] ! particle at (5,5) should be in ghost region

        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 2), "Number of ghost particles generated for right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), 6._fp), "1st right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 3._fp), "1st right-wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), 5._fp), "2nd right-wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 4._fp), "2nd right-wall ghost particle y-coordinate incorrect.")

    end subroutine test_coaxial_ghost_walls

    subroutine test_diagonal_ghost_walls()

        type(eos_particles_t):: ps_template
        type(eos_ghost_particles_t):: ps_ghost_template
        type(particle_system_t):: psys, psys_ghost
        type(ghost_timestep_setuper_t):: ghost_setuper
        type(system_interaction_t):: interaction
        type(sweeper_container_t):: sweepers(0)
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
            psys%particles%id(i) = i
            psys%particles%x(1, i) = 1._fp
            psys%particles%x(2, i) = real(i, kind=fp)
        end do

        ! test bottom-right wall (y = x - 1)
        ghost_setuper%cutoff = 2._fp
        ghost_setuper%surface_normal = [-sqrt(0.5_fp), sqrt(0.5_fp)] ! normal vector pointing top-left
        ghost_setuper%point = [0._fp, -1._fp]          ! plane/line travels through (0, 0)

        ! initialize system interaction with ghost setuper
        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 2), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), 2._fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 0._fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), 3._fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 0._fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle y-coordinate incorrect.")

        ! check particle in ghost region isn't mirrored
        ! move boundary to y = x + 0.5 to make particle at (1, 1) in ghost region
        ghost_setuper%point = [0._fp, 0.5_fp]
        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 3), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), 1.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 1.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), 2.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 1.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 3), 3.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 3), 1.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle y-coordinate incorrect.")

        ! test top-left wall (y = x + 5)
        ghost_setuper%surface_normal = [sqrt(0.5_fp), -sqrt(0.5_fp)]
        ghost_setuper%point = [0._fp, 5._fp]

        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 2), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), -1._fp, atol=1.d-14), &
                   "1st top-left wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 6._fp, atol=1.d-14), &
                   "1st top-left wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), 0._fp, atol=1.d-14), &
                   "2nd top-left wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 6._fp, atol=1.d-14), &
                   "2nd top-left wall ghost particle y-coordinate incorrect.")

        ! check particle in ghost region isn't mirrored
        ! move boundary to y = x + 3.5 to make particle at (5, 5) in ghost region
        ghost_setuper%point = [0._fp, 3.5_fp]
        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 3), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), -1.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 4.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), -0.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 4.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 3), 0.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 3), 4.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle y-coordinate incorrect.")

        ! test bottom-left wall (y = 1 - x)
        ghost_setuper%surface_normal = [sqrt(0.5_fp), sqrt(0.5_fp)] ! normal vector pointing top-right
        ghost_setuper%point = [0._fp, 1._fp]                        ! plane/line travels through (0, 5)

        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 2), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), 0._fp, atol=1.d-14), &
                   "1st bottom-left wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 0._fp, atol=1.d-14), &
                   "1st bottom-left wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), -1._fp, atol=1.d-14), &
                   "2nd bottom-left wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 0._fp, atol=1.d-14), &
                   "2nd bottom-left wall ghost particle y-coordinate incorrect.")

        ! check particle in ghost region isn't mirrored
        ! move boundary to y = 2.5 - x to make particle at (1, 1) in ghost region
        ghost_setuper%point = [0._fp, 2.5_fp]
        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 3), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), 0.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 1.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), -0.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 1.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 3), -1.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 3), 1.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle y-coordinate incorrect.")

        ! test top-right wall (y = 7 - x)
        ! ghost_setuper%cutoff = 1._fp
        ghost_setuper%surface_normal = [-sqrt(0.5_fp), -sqrt(0.5_fp)] ! normal vector pointing top-right
        ghost_setuper%point = [0._fp, 7._fp]                          ! plane/line travels through (0, 8)

        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 2), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), 3._fp, atol=1.d-14), &
                   "1st top-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 6._fp, atol=1.d-14), &
                   "1st top-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), 2._fp, atol=1.d-14), &
                   "2nd top-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 6._fp, atol=1.d-14), &
                   "2nd top-right wall ghost particle y-coordinate incorrect.")

        ! check particle in ghost region isn't mirrored
        ! move boundary to y = 5.5 - x to make particle at (5, 5) in ghost region
        ghost_setuper%point = [0._fp, 5.5_fp]
        call interaction%init(1, psys, psys_ghost, sweepers=sweepers, timestep_setuper=ghost_setuper)

        call interaction%do_timestep_setup()

        call check(is_equal(psys_ghost%size(), 3), "Number of ghost particles generated for bottom-right wall incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 1), 3.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 1), 4.5_fp, atol=1.d-14), &
                   "1st bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 2), 2.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 2), 4.5_fp, atol=1.d-14), &
                   "2nd bottom-right wall ghost particle y-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(1, 3), 1.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle x-coordinate incorrect.")
        call check(is_close(psys_ghost%particles%x(2, 3), 4.5_fp, atol=1.d-14), &
                   "3rd bottom-right wall ghost particle y-coordinate incorrect.")

    end subroutine test_diagonal_ghost_walls

    subroutine test_morris_walls()

        type(particle_system_t):: fluid_psys, wall_psys
        type(eos_particles_t):: fluid_template
        type(system_interaction_t):: interaction
        type(sweeper_container_t):: sweepers(1)
        type(cubic_bspline_kernel_t):: kernel
        type(morris_boundary_sweeper_t), pointer:: sweeper
        integer:: i

        ! position fluid particle
        call fluid_psys%init(1, "fluid", fluid_template)
        fluid_psys%particles%x(:, 1) = [(real(i, kind=fp)/10._fp, i=1, ndims)] ! [0.1, 0.2, 0.3]
        fluid_psys%particles%v(:, 1) = [(2._fp*real(i, kind=fp), i=1, ndims)] ! [2, 4, 6]
        fluid_psys%particles%rho(1) = 1000._fp
        fluid_psys%particles%mass(1) = 10._fp
        fluid_psys%particles%c(1) = 200._fp
        select type (ps => fluid_psys%particles)
        class is (eos_particles_t)
            ps%p(1) = 5._fp
        end select

        ! position wall particles
        call wall_psys%init(1, "wall")
        wall_psys%particles%x(:, 1) = -fluid_psys%particles%x(:, 1) ! [-0.1, -0.2, -0.3]

        call kernel%init(ndims, 1.2_fp)

        allocate (morris_boundary_sweeper_t::sweepers(1)%sweeper)
        select type (s => sweepers(1)%sweeper)
        type is (morris_boundary_sweeper_t)
            sweeper => s
        end select

        ! test horizontal wall
        sweeper%h = kernel%h
        sweeper%artvisc_alpha = 0.1_fp
        sweeper%artvisc_beta = 0.1_fp
        sweeper%normal(:) = 0._fp
        sweeper%normal(ndims) = 1._fp
        sweeper%point(:) = 0._fp
        call interaction%init(1, fluid_psys, wall_psys, sweepers=sweepers)

        call run_check(interaction, kernel, 1._fp, sweepers, "")

        ! move wall particle further away
        wall_psys%particles%x(:, 1) = -2._fp*fluid_psys%particles%x(:, 1) ! [-0.2, -0.4, -0.6]

        call run_check(interaction, kernel, expected_dbda=2._fp, sweepers=sweepers, &
                       msg_suffix="after changing wall particle position.")

        ! horizontal upper wall
        fluid_psys%particles%x(ndims, 1) = 40._fp - fluid_psys%particles%x(ndims, 1)
        wall_psys%particles%x(:, 1) = [(-real(i, kind=fp)/10._fp, i=1, ndims)]
        wall_psys%particles%x(ndims, 1) = 40._fp - wall_psys%particles%x(ndims, 1)

        sweeper%normal(ndims) = -sweeper%normal(ndims)
        sweeper%point(ndims) = 40._fp

        call interaction%init(1, fluid_psys, wall_psys, sweepers=sweepers)

        call run_check(interaction, kernel, 1._fp, sweepers, "for horizontal upper wall.")

        ! move wall particle further away
        wall_psys%particles%x(:, 1) = [(-real(i, kind=fp)/5._fp, i=1, ndims)]
        wall_psys%particles%x(ndims, :) = 40._fp - wall_psys%particles%x(ndims, :)

        call run_check(interaction, kernel, 2._fp, sweepers, "for horizontal upper wall.")

        ! vertical lower wall
        fluid_psys%particles%x(:, 1) = [(real(i, kind=fp)/10._fp, i=1, ndims)]
        wall_psys%particles%x(:, 1) = -fluid_psys%particles%x(:, 1)
        sweeper%normal(:) = 0._fp
        sweeper%normal(1) = 1._fp
        sweeper%point(:) = 0._fp

        call interaction%init(1, fluid_psys, wall_psys, sweepers=sweepers)

        call run_check(interaction, kernel, 1._fp, sweepers, "for vertical lower wall.")

        ! vertical upper wall
        fluid_psys%particles%x(:, 1) = [(real(i, kind=fp)/10._fp, i=1, ndims)]
        fluid_psys%particles%x(ndims, :) = 40._fp - fluid_psys%particles%x(ndims, :)
        wall_psys%particles%x(:, 1) = [(real(i, kind=fp)/10._fp, i=1, ndims)]
        wall_psys%particles%x(ndims, :) = 40._fp + wall_psys%particles%x(ndims, :)

        sweeper%normal(:) = 0._fp
        sweeper%normal(1) = 1._fp
        sweeper%point(:) = 0._fp
        sweeper%point(1) = 40._fp

        call interaction%init(1, fluid_psys, wall_psys, sweepers=sweepers)

        call run_check(interaction, kernel, 1._fp, sweepers, "for vertical upper wall.")

        ! test angular wall
        fluid_psys%particles%x(:, 1) = [(real(i, kind=fp)/10._fp, i=1, ndims)]
        wall_psys%particles%x(:, 1) = -fluid_psys%particles%x(:, 1)
        sweeper%normal(:) = 0._fp
        sweeper%normal(1) = 2._fp/sqrt(5._fp)
        sweeper%normal(2) = 1._fp/sqrt(5._fp)
        sweeper%point(:) = 0._fp

        call interaction%init(1, fluid_psys, wall_psys, sweepers=sweepers)

        call run_check(interaction, kernel, 1._fp, sweepers, "for vertical upper wall.")

        ! move wall particle a little
        wall_psys%particles%x(1, 1) = 2._fp*wall_psys%particles%x(1, 1)

        call run_check(interaction, kernel, 1.5_fp, sweepers, "for vertical upper wall, after moving wall.")

    end subroutine test_morris_walls

    subroutine run_check(interaction, kernel, expected_dbda, sweepers, msg_suffix)

        use grasph_pair_interactions_m, only: continuity_density, artificial_viscosity_monaghan1994, isotropic_pressure_force

        type(system_interaction_t), intent(inout):: interaction
        type(cubic_bspline_kernel_t), intent(in):: kernel
        real(fp), intent(in):: expected_dbda
        class(sweeper_container_t), intent(in):: sweepers(1)
        character(*), intent(in):: msg_suffix
        class(eos_particles_t), pointer:: ps_fluid
        class(base_particles_t), pointer:: ps_wall
        real(fp):: w, dwdx(ndims), expected_wall_v(ndims), expected_wall_rho, expected_wall_mass, expected_wall_p, &
                   expected_wall_c, expected_drhodt, expected_dvxdt(ndims), dummy_dvxdt(ndims), dummy_drhodt
        integer:: d
        character:: dc

        if (.not. interaction%is_pair_set) error stop "Not pair set"
        if (interaction%psys_lhs%size() /= 1) error stop "Too many particles in psys_lhs"
        if (interaction%psys_rhs%size() /= 1) error stop "Too many particles in psys_rhs"

        select type (ps => interaction%psys_lhs%particles)
        class is (eos_particles_t)
            ps_fluid => ps
        class default
            error stop "interaction%psys_lhs%particles is not class eos_particles_t"
        end select
        ps_wall => interaction%psys_rhs%particles

        call kernel%values( &
            ps_fluid%x(:, 1) - ps_wall%x(:, 1), &
            w, dwdx &
            )

        expected_wall_v(:) = -expected_dbda*ps_fluid%v(:, 1)
        expected_wall_rho = ps_fluid%rho(1)
        expected_wall_mass = ps_fluid%mass(1)
        expected_wall_p = ps_fluid%p(1)
        expected_wall_c = ps_fluid%c(1)

        ps_fluid%drhodt(1) = 0._fp
        ps_fluid%dvxdt(:, 1) = 0._fp
        expected_drhodt = 0._fp
        expected_dvxdt = 0._fp

        call interaction%find_pairs(kernel%cutoff, kernel)

        call check(is_equal(interaction%pairs%npairs_total, 1))
        call interaction%do_sweep(1)

        call continuity_density( &
            ps_fluid%v(:, 1), &
            expected_wall_v(:), &
            ps_fluid%mass(1), &
            expected_wall_mass, &
            expected_drhodt, &
            dummy_drhodt, &
            dwdx(:) &
            )
        select type (sweeper => sweepers(1)%sweeper)
        type is (morris_boundary_sweeper_t)
            call artificial_viscosity_monaghan1994( &
                ps_fluid%x(:, 1), &
                ps_wall%x(:, 1), &
                ps_fluid%v(:, 1), &
                expected_wall_v, &
                ps_fluid%rho(1), &
                expected_wall_rho, &
                kernel%h, &
                kernel%h, &
                ps_fluid%c(1), &
                expected_wall_c, &
                ps_fluid%mass(1), &
                expected_wall_mass, &
                expected_dvxdt(:), &
                dummy_dvxdt(:), &
                dwdx, &
                sweeper%artvisc_alpha, &
                sweeper%artvisc_beta &
                )
        end select
        call isotropic_pressure_force( &
            ps_fluid%p(1), &
            expected_wall_p, &
            ps_fluid%rho(1), &
            expected_wall_rho, &
            ps_fluid%mass(1), &
            expected_wall_mass, &
            expected_dvxdt, &
            dummy_dvxdt, &
            dwdx &
            )

        call check(is_close(expected_drhodt, ps_fluid%drhodt(1)), "Incorrect drhodt "//msg_suffix)

        do d = 1, ndims
            write (dc, "(I1)") d
            call check(is_close(expected_dvxdt(d), ps_fluid%dvxdt(d, 1)), "Incorrect dvxdt for d="//dc//" "//msg_suffix)
        end do

    end subroutine run_check

end module test_boundary

program run_tests

    use test_boundary, only: tests
    use fortuno_serial, only: execute => execute_serial_cmd_app
    implicit none

    call execute(tests())

end program run_tests
