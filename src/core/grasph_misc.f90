!> @file grasph_misc.f90
!> @brief Module miscelaneous helper routines and types.
!> @author Edward Yang
!> @date 2025-06-09
module grasph_misc

    use grasph_constants, only: fp
    use grasph_particles, only: base_particles
    use grasph_pairs, only: particle_pairs
    use iso_fortran_env, only: int64, real64

    implicit none
    private

    !> @brief A simple timer derived type that uses the more accurate system_clock intrinsic to
    !>        measure elapsed time.
    type system_timer
        !> @brief System clock counter when the timer was started.
        integer(int64):: start_count
        !> @brief System clock rate set when timer is started.
        integer(int64):: rate
        !> @brief A counter for number of particle interactions. Used to calculate interactions per wall-second.
        integer(int64):: ninteractions = 0
    contains
        !> @brief Starts the system clock timer.
        procedure:: start => start_timer
        !> @brief Stop the system clock timer and returns the elapsed time.
        procedure:: stop => stop_timer
        !> @brief Updates the interaction counter.
        procedure:: update_interactions
    end type system_timer

    public:: print_summary, system_timer

contains

    !> @brief A generic helper function to print summary data about the running simulation.
    !>        It leverages particles' "generate_summary" procedures to print summary data about the
    !>        the particles to the terminal.
    !> @param itimestep The current time-step.
    !> @param time_integration_scheme The name of the time-integration scheme used in the simulation.
    !> @param particles The list of particles whose generate_summary methods to use.
    !> @param timer The system_timer object used to track time.
    subroutine print_summary(itimestep, time_integration_scheme, particles, timer, time)

        integer, intent(in):: itimestep
        character(*), intent(in):: time_integration_scheme
        type(system_timer), optional, intent(in):: timer
        class(base_particles), intent(in):: particles(:)
        real(fp), intent(in):: time
        character(:), allocatable:: psummary
        integer:: i
        real(real64):: t
        character(20):: tc, nips

        write (*, "(A)") "========================= GraSPH Output ========================="
        write (*, "(A, I13)") time_integration_scheme//" time-intregration, time-step: ", itimestep
        write (*, "(A, f14.7)") "  In-simulation time: ", time
        do i = 1, size(particles)
            write (*, "(A)") "  Summary data for: "//trim(particles(i)%name)
            if (particles(i)%to_print_summary) then
                call particles(i)%generate_summary(psummary)
                write (*, "(A)") psummary
            else
                write (*, "(A)") "    skipped"
            end if
        end do
        if (present(timer)) then
            t = timer%stop()
            write (tc, '(f12.5)') t
            write (*, "(A)") "Elapsed wall time (s)       : "//trim(adjustl(tc))
            write (nips, "(es15.3)") real(timer%ninteractions, kind=real64)/t
            write (*, "(A)") "Interactions per wall-second: "//trim(adjustl(nips))
        end if
        write (*, "(A)") "================================================================="

    end subroutine print_summary

    !> @brief Starts the system_timer by recording the current system_clock time.
    !> @param self The timer to start.
    subroutine start_timer(self)
        class(system_timer), intent(out):: self

        call system_clock(count=self%start_count, count_rate=self%rate)
    end subroutine start_timer

    !> @brief Returns the elapsed time, relative to the currently recorded start time.
    !> @param self The timer.
    !> @returns elapsed The time elapsed since the given timer was started.
    real(real64) function stop_timer(self) result(elapsed)
        class(system_timer), intent(in):: self
        integer(int64):: end_count
        call system_clock(count=end_count)
        elapsed = real(end_count - self%start_count, kind=real64)/real(self%rate, kind=real64)
    end function stop_timer

    !> @brief Updates the interaction counter.
    !> @param The timer object to update.
    !> @param The number of interactions to increment by.
    subroutine update_interactions(self, ninteractions)
        class(system_timer), intent(inout):: self
        integer, intent(in):: ninteractions

        self%ninteractions = self%ninteractions + int(ninteractions, kind=int64)
    end subroutine update_interactions

end module grasph_misc
