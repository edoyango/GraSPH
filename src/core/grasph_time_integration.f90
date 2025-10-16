!> @file grasph_time_integration.f90
!> @brief Module containing useful generic time-integration schemes to be used in an SPH simulation.
!> @author Edward Yang
!> @date 2025-06-09
module grasph_time_integration_m

    use grasph_constants_m, only: fp
    use grasph_particle_system_m, only: particle_system_t
    use grasph_system_interactions_m, only: system_interaction_t
    use grasph_kernels_m, only: base_kernel_t
    use grasph_misc_m, only: print_summary, system_timer_t
    use grasph_common_m, only: array_pointer_container_t

    implicit none
    private

    public:: leap_frog_time_integration

contains

    !> @brief Leap-Frog time-integration (e.g. https://en.wikipedia.org/wiki/Leapfrog_integration).
    !> @param maxtimestep The maximum number of time-steps to run the time-integration for.
    !> @param print_step The interval number of time-steps to update the terminal with run information.
    !> @param save_step The interval number of time-steps to save data to disk using particles' inbuilt dump method.
    !> @param particles The list of particles who's state is being updated over time.
    !> @param interactions The list of particle_pair_set objects which describe particles' relationship with one another.
    !> @param CFL The Courant-Freidrichs-Lewy coefficient for time-stepping.
    !> @param kernel The kernel to use.
    !> @param output_path The directory to store saved data.
    !> @param output_prefix The filename prefix to use in the output files.
    !> @param output_comp_level The level of GZIP compression to use when writing the output HDF5 files.
    subroutine leap_frog_time_integration(maxtimestep, print_step, save_step, psystems, interactions, CFL, kernel, &
                                          output_path, output_prefix, output_comp_level, damping_coef)

        integer, intent(in):: maxtimestep, print_step, save_step
        class(particle_system_t):: psystems(:)
        class(system_interaction_t):: interactions(:)
        real(fp), intent(in):: CFL
        class(base_kernel_t), intent(in):: kernel
        character(*), intent(in):: output_path
        character(*), optional, intent(in):: output_prefix
        integer, optional, intent(in):: output_comp_level
        real(fp), optional, intent(in):: damping_coef
        integer:: nparticle_sets, nparticle_interactions, itimestep, i, j, k, max_registrations
        real(fp):: dt, time, maxc
        type(system_timer_t):: timer
        type(array_pointer_container_t), allocatable:: vars0(:, :)
        real(fp), pointer:: var_ptr(:), deriv_ptr(:)

        nparticle_sets = size(psystems)
        nparticle_interactions = size(interactions)

        max_registrations = 0
        do i = 1, nparticle_sets
            max_registrations = max(psystems(i)%register_v%nregistrations, max_registrations)
        end do
        allocate (vars0(max_registrations, nparticle_sets))

        do i = 1, nparticle_sets
            do j = 1, psystems(i)%register_v%nregistrations
                allocate (vars0(j, i)%p(psystems(i)%register_v%dims(j), psystems(i)%size))
            end do
        end do

        time = 0._fp

        call timer%start()

        do itimestep = 1, maxtimestep

            ! perform any setup needed for each particle interaction e.g. create ghost particles
            do i = 1, nparticle_interactions
                call interactions(i)%do_timestep_setup
            end do

            ! calculate timestep to use
            maxc = psystems(1)%particles(1)%c
            do i = 1, nparticle_sets
                do j = 1, psystems(i)%size
                    maxc = max(maxc, psystems(i)%particles(j)%c)
                end do
            end do
            dt = CFL*kernel%h/maxc

            ! save data at start of timestep for each particles
            do i = 1, nparticle_sets
                do j = 1, psystems(i)%register_v%nregistrations
                    do k = 1, psystems(i)%size
                        call psystems(i)%register_v%get(psystems(i)%particles(k), j, var_ptr, deriv_ptr)
                        vars0(j, i)%p(:, k) = var_ptr(:)
                    end do
                end do
            end do

            ! find pairs between provided particle interaction sets
            do i = 1, nparticle_interactions
                call interactions(i)%find_pairs(kernel%cutoff, kernel)
            end do

            ! update particles to mid-timestep
            do i = 1, nparticle_sets
                do j = 1, psystems(i)%register_v%nregistrations
                    if (present(damping_coef) .and. trim(psystems(i)%register_v%names(j)) == "v") then
                        do k = 1, psystems(i)%size
                            psystems(i)%particles(k)%v(:) = (1._fp - 0.5_fp*damping_coef*dt)*psystems(i)%particles(k)%v(:)
                        end do
                    end if
                    do k = 1, psystems(i)%size
                        call psystems(i)%register_v%get(psystems(i)%particles(k), j, var_ptr, deriv_ptr)
                        var_ptr(:) = var_ptr(:) + 0.5_fp*dt*deriv_ptr(:)
                    end do
                end do
            end do

            ! perform pre-sweep prologue e.g. to update boundary particles' state
            do i = 1, nparticle_interactions
                call interactions(i)%do_sweep_prologue
            end do

            ! Update particle state e.g. pressure/stress
            do i = 1, nparticle_sets
                call psystems(i)%do_state_update(0.5_fp*dt)
            end do

            ! perform actual sweep i.e., calculate acceleration, density change etc.
            do i = 1, nparticle_interactions
                call interactions(i)%do_sweep
            end do

            ! update states to full-timestep
            do i = 1, nparticle_sets
                do j = 1, psystems(i)%register_v%nregistrations
                    do k = 1, psystems(i)%size
                        call psystems(i)%register_v%get(psystems(i)%particles(k), j, var_ptr, deriv_ptr)
                        var_ptr(:) = vars0(j, i)%p(:, k) + dt*deriv_ptr(:)
                    end do
                    if (present(damping_coef) .and. trim(psystems(i)%register_v%names(j)) == "v") then
                        do k = 1, psystems(i)%size
                            psystems(i)%particles(k)%v(:) = (1._fp - damping_coef*dt)*psystems(i)%particles(k)%v(:)
                        end do
                    end if
                end do
                do j = 1, psystems(i)%register_x%nregistrations
                    do k = 1, psystems(i)%size
                        call psystems(i)%register_x%get(psystems(i)%particles(k), j, var_ptr, deriv_ptr)
                        var_ptr(:) = var_ptr(:) + dt*deriv_ptr(:)
                    end do
                end do
            end do

            ! perform shifting
            do i = 1, nparticle_interactions
                call interactions(i)%do_shift(dt)
            end do

            ! write data
            if (mod(itimestep, save_step) == 0) then
                do i = 1, nparticle_sets
                    call psystems(i)%dump(itimestep, output_path, output_prefix, output_comp_level)
                end do
            end if

            time = time + dt

            ! update number of interactions over loop lifetime
            do i = 1, nparticle_interactions
                call timer%update_interactions(interactions(i)%pairs%npairs_total)
            end do

            ! print data to screen
            if (mod(itimestep, print_step) == 0) then
                call print_summary(itimestep, "Leap-Frog", psystems, timer, time)
            end if

        end do

        deallocate (vars0)

    end subroutine leap_frog_time_integration

end module grasph_time_integration_m
