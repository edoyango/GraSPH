module grasph_misc

    use grasph_particles, only: particles_container

    implicit none
    private
    public:: print_summary

contains

    subroutine print_summary(itimestep, time_integration_scheme, particles)

        integer, intent(in):: itimestep
        character(*), intent(in):: time_integration_scheme
        class(particles_container):: particles(:)
        character(:), allocatable:: psummary
        integer:: i

        write(*, "(A)") "========================= GraSPH Output ========================="
        write(*, "(A, I13)") time_integration_scheme // " time-intregration, time-step: ", itimestep
        do i = 1, size(particles)
            write(*, "(A)") "  Summary data for: " // trim(particles(i)%p%name)
            call particles(i)%p%generate_summary(psummary)
            write(*, "(A)") psummary
        enddo
        write(*, "(A)") "================================================================="

    end subroutine print_summary

end module grasph_misc