!> @file grasph_time_integration.f90
!> @brief Module containing useful generic time-integration schemes to be used in an SPH simulation.
!> @author Edward Yang
!> @date 2025-06-09
module grasph_time_integration

    use grasph_constants, only: fp
    use grasph_particles, only: particles_container
    use grasph_pair_sets, only: particle_interactions
    use grasph_kernels, only: grasph_base_kernel
    use grasph_misc, only: print_summary, system_timer

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
    subroutine leap_frog_time_integration(maxtimestep, print_step, save_step, particles, interactions, CFL, kernel, &
                                          output_path, output_prefix, output_comp_level)

        integer, intent(in):: maxtimestep, print_step, save_step
        class(particles_container):: particles(:)
        class(particle_interactions):: interactions(:)
        real(fp), intent(in):: CFL
        class(grasph_base_kernel), intent(in):: kernel
        character(*), intent(in):: output_path
        character(*), optional, intent(in):: output_prefix
        integer, optional, intent(in):: output_comp_level
        integer:: nparticle_sets, nparticle_interactions, itimestep, i
        real(fp):: dt
        type(system_timer):: timer

        nparticle_sets = size(particles)
        nparticle_interactions = size(interactions)

        call timer%start()

        do itimestep = 1, maxtimestep

            ! calculate timestep to use
            dt = huge(1._fp)
            do i = 1, nparticle_sets
                dt = min(dt, CFL*kernel%h/maxval(particles(i)%p%c(:)))
            end do

            ! save data at start of timestep for each particles
            do i = 1, nparticle_sets
                call particles(i)%p%start_timestep()
            end do

            ! find pairs between provided particle interaction sets
            do i = 1, nparticle_interactions
                call interactions(i)%find_pairs(kernel%cutoff, kernel)
            end do

            ! update particles to mid-timestep
            do i = 1, nparticle_sets
                if (particles(i)%p%evolve) call particles(i)%p%mid_timestep_update(0.5_fp*dt)
            end do

            ! perform pre-sweep prologue e.g. to update boundary particles' state
            do i = 1, nparticle_interactions
                call interactions(i)%do_sweep_prologue
            end do

            ! Update particle state e.g. pressure/stress
            do i = 1, nparticle_sets
                call particles(i)%p%state_update(0.5_fp*dt)
            end do

            ! perform actual sweep i.e., calculate acceleration, density change etc.
            do i = 1, nparticle_interactions
                call interactions(i)%do_sweep
            end do

            ! update states to full-timestep
            do i = 1, nparticle_sets
                if (particles(i)%p%evolve) call particles(i)%p%full_timestep_update(dt, update_position=.true.)
            end do

            ! perform shifting
            do i = 1, nparticle_interactions
                call interactions(i)%do_shift(dt)
            end do

            ! write data
            if (mod(itimestep, save_step) == 0) then
                do i = 1, nparticle_sets
                    call particles(i)%p%dump(itimestep, output_path, output_prefix, output_comp_level)
                end do
            end if

            ! print data to screen
            if (mod(itimestep, print_step) == 0) then
                call print_summary(itimestep, "Leap-Frog", particles, timer)
            end if

        end do

    end subroutine leap_frog_time_integration

end module grasph_time_integration
