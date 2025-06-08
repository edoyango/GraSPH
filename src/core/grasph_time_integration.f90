module grasph_time_integration

    use grasph_constants, only: fp
    use grasph_particles, only: particles_container
    use grasph_pair_sets, only: particle_interactions_container
    use grasph_kernels, only: grasph_base_kernel
    
    implicit none
    private

    public:: leap_frog_time_integration

contains

    subroutine leap_frog_time_integration(maxtimestep, print_step, save_step, particles, particle_interactions, CFL, kernel, &
                                          output_path, output_prefix, output_comp_level)

        integer, intent(in):: maxtimestep, print_step, save_step
        class(particles_container):: particles(:)
        class(particle_interactions_container):: particle_interactions(:)
        real(fp), intent(in):: CFL
        class(grasph_base_kernel), intent(in):: kernel
        character(*), intent(in):: output_path, output_prefix
        integer, intent(in):: output_comp_level
        integer:: nparticle_sets, nparticle_interactions, itimestep, i
        real(fp):: dt

        nparticle_sets = size(particles)
        nparticle_interactions = size(particle_interactions)

        do itimestep = 1, maxtimestep

            ! calculate timestep to use
            dt = huge(1._fp)
            do i = 1, nparticle_sets
                dt = min(dt, CFL*kernel%h/maxval(particles(i)%p%c(:)))
            enddo

            ! find pairs between provided particle interaction sets
            do i = 1, nparticle_interactions
                call particle_interactions(i)%pi%find_pairs(kernel%cutoff, kernel)
            enddo

            ! save data at start of timestep for each particles
            do i = 1, nparticle_sets
                call particles(i)%p%start_timestep()
            enddo

            ! update particles to mid-timestep
            do i = 1, nparticle_sets
                call particles(i)%p%mid_timestep_update(0.5_fp*dt)
                call particles(i)%p%state_update(0.5_fp*dt)
            enddo

            ! perform pre-sweep prologue e.g. to update boundary particles' state
            do i = 1, nparticle_interactions
                call particle_interactions(i)%pi%sweep_prologue()
            enddo
            ! perform actual sweep i.e., calculate acceleration, density change etc.
            do i = 1, nparticle_interactions
                call particle_interactions(i)%pi%sweep
            enddo

            ! update states to full-timestep
            do i = 1, nparticle_sets
                call particles(i)%p%full_timestep_update(dt, update_position=.true.)
            enddo

            ! write data
            if (mod(itimestep, save_step) == 0) then
                do i = 1, nparticle_sets
                    call particles(i)%p%dump(itimestep, output_path, output_prefix, output_comp_level)
                enddo
            endif

            ! print data to screen
            if (mod(itimestep, print_step) == 0) then
                write(*, "(A)")       "------------------------- GraSPH Output -------------------------"
                write(*, "(A,I13,A)") "                         ", itimestep, "                         "
            endif

        enddo

    end subroutine leap_frog_time_integration

end module grasph_time_integration