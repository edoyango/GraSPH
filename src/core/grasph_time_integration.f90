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
    !> @param psystems The list of particle systems who's state is being updated over time.
    !> @param interactions The list of particle_pair_set objects which describe particles' relationship with one another.
    !> @param CFL The Courant-Freidrichs-Lewy coefficient for time-stepping.
    !> @param kernel The kernel to use.
    !> @param output_path The directory to store saved data.
    !> @param output_prefix The filename prefix to use in the output files.
    !> @param output_comp_level The level of GZIP compression to use when writing the output HDF5 files.
    !> @param damping_coef The strength of damping to apply during time stepping.
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
        integer:: nparticle_sets, nparticle_interactions, itimestep, i, j, k, max_registrations, n_sweep_updates, istage
        real(fp):: dt, time, maxc
        type(system_timer_t):: timer
        type(array_pointer_container_t), allocatable:: vars0(:, :)
        real(fp), pointer:: var_ptr(:), deriv_ptr(:)

        n_sweep_updates = count_sweeps_and_updates(psystems, interactions)

        nparticle_sets = size(psystems)
        nparticle_interactions = size(interactions)

        max_registrations = 0
        do i = 1, nparticle_sets
            max_registrations = max(psystems(i)%register_v%nregistrations, max_registrations)
        end do
        allocate (vars0(max_registrations, nparticle_sets))

        do i = 1, nparticle_sets
            do j = 1, psystems(i)%register_v%nregistrations
                if (psystems(i)%register_v%dims(j) > 1) then
                    allocate (vars0(j, i)%p(psystems(i)%register_v%dims(j), psystems(i)%size()))
                else
                    allocate (vars0(j, i)%p(psystems(i)%size(), 1))
                end if
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
            maxc = psystems(1)%particles%c(1)
            do i = 1, nparticle_sets
                do j = 1, psystems(i)%size()
                    maxc = max(maxc, psystems(i)%particles%c(j))
                end do
            end do
            dt = CFL*kernel%h/maxc

            ! save data at start of timestep for each particles
            do i = 1, nparticle_sets
                do j = 1, psystems(i)%register_v%nregistrations
                    vars0(j, i)%p(:, :) = psystems(i)%register_v%variables(j)%p(:, :)
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
                        psystems(i)%particles%v(:, :) = (1._fp - 0.5_fp*damping_coef*dt)*psystems(i)%particles%v(:, :)
                    end if
                    psystems(i)%register_v%variables(j)%p(:, :) = psystems(i)%register_v%variables(j)%p(:, :) + &
                                                                  0.5_fp*dt*psystems(i)%register_v%derivatives(j)%p(:, :)
                end do
            end do

            do istage = 1, n_sweep_updates
                do i = 1, nparticle_sets
                    if (allocated(psystems(i)%state_updaters) .and. size(psystems(i)%state_updaters) > 0) &
                        call psystems(i)%do_state_update(istage, 0.5_fp*dt)
                end do
                do i = 1, nparticle_interactions
                    call interactions(i)%do_sweep(istage)
                end do
            end do

            ! update states to full-timestep
            do i = 1, nparticle_sets
                do j = 1, psystems(i)%register_v%nregistrations

                    psystems(i)%register_v%variables(j)%p(:, :) = vars0(j, i)%p(:, :) + &
                                                                  dt*psystems(i)%register_v%derivatives(j)%p(:, :)

                    if (present(damping_coef) .and. trim(psystems(i)%register_v%names(j)) == "v") then
                        psystems(i)%particles%v(:, :) = (1._fp - damping_coef*dt)*psystems(i)%particles%v(:, :)
                    end if
                end do
                do j = 1, psystems(i)%register_x%nregistrations
                    psystems(i)%register_x%variables(j)%p(:, :) = psystems(i)%register_x%variables(j)%p(:, :) + &
                                                                  dt*psystems(i)%register_x%derivatives(j)%p(:, :)
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

    !> @brief Ensures the count of 'system interaction sweepers' matches the count of 'particle system state updaters'.
    !>        A particle system may have 0 zero state updaters (i.e. unallocated or allocated and size 0).
    integer function count_sweeps_and_updates(psystems, sinteractions)
        type(particle_system_t), intent(in):: psystems(:)
        type(system_interaction_t), intent(in):: sinteractions(:)
        integer:: i, nsweeps, nupdates

        nsweeps = size(sinteractions(1)%sweepers)

        do i = 2, size(sinteractions)
            if (nsweeps /= size(sinteractions(i)%sweepers)) error stop "Not all interactions have same number of sweepers."
        end do

        do i = 1, size(psystems)
            if (allocated(psystems(i)%state_updaters)) then
                nupdates = size(psystems(i)%state_updaters)
                if (nupdates > 0 .and. nupdates /= nsweeps) &
                    error stop "Number of state updaters in particle systems don't match number of sweeps."
            end if
        end do

        count_sweeps_and_updates = nsweeps

    end function count_sweeps_and_updates

end module grasph_time_integration_m
