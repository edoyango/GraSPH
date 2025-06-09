module grasph_misc

    use grasph_particles, only: particles_container
    use iso_fortran_env, only: int64, real64

    implicit none
    private

    type system_timer
        integer(int64):: start_count, rate
    contains
        procedure:: start => start_timer, stop => stop_timer
    end type system_timer

    public:: print_summary, system_timer

contains

    subroutine print_summary(itimestep, time_integration_scheme, particles, timer)

        integer, intent(in):: itimestep
        character(*), intent(in):: time_integration_scheme
        type(system_timer), optional, intent(in):: timer
        class(particles_container):: particles(:)
        character(:), allocatable:: psummary
        integer:: i

        write(*, "(A)") "========================= GraSPH Output ========================="
        write(*, "(A, I13)") time_integration_scheme // " time-intregration, time-step: ", itimestep
        do i = 1, size(particles)
            write(*, "(A)") "  Summary data for: " // trim(particles(i)%p%name)
            if (particles(i)%p%to_print_summary) then
                call particles(i)%p%generate_summary(psummary)
                write(*, "(A)") psummary
            else
                write(*, "(A)") "    skipped"
            endif
        enddo
        if (present(timer)) then
            write(*, "(A, f12.5)") "Elapsed wall time: ", timer%stop()
        endif
        write(*, "(A)") "================================================================="

    end subroutine print_summary

    subroutine start_timer(self)
        class(system_timer), intent(out):: self

        call system_clock(count=self%start_count, count_rate=self%rate)
    end subroutine start_timer

    real(real64) function stop_timer(self) result(elapsed)
        class(system_timer), intent(in):: self
        integer(int64):: end_count
        call system_clock(count=end_count)
        elapsed = real(end_count - self%start_count, kind=real64)/real(self%rate, kind=real64)
    end function stop_timer

end module grasph_misc