!> @file grasph_register.f90
!> @brief Module containing classes for registering particle variables e.g. for io or time-evolution.
!> @author Edward Yang
!> @date 2025-10-02
module grasph_register

    use grasph_constants, only: fp, ndims
    use grasph_particle, only: base_particle
    use iso_c_binding, only: c_intptr_t, c_f_pointer, c_ptr, c_loc

    implicit none

    private

    !> @brief Max allowable registrations. Will be removed for a more dynamic approach.
    integer, parameter:: max_registrations = 20

    !> @brief Registers particle's variables for access through pointer e.g. for dynamically selecting variables for IO.
    type variable_register
        !> @brief Number of variables currently registered.
        integer:: nregistrations = 0
        !> @brief Dimension of each variable registered.
        integer:: dims(max_registrations)
        !> @brief Name of variables registere. Used in IO.
        character(20):: names(max_registrations)
        !> @brief Offset in memory of registered variables.
        integer(c_intptr_t):: offsets(max_registrations)
    contains
        !> @brief Registers a scalar member of a particle e.g. density.
        procedure, private:: register_variable_vector
        !> @brief Registers a vector member of a particle e.g. position.
        procedure, private:: register_variable_scalar
        !> @brief Registers a member of a particle.
        generic, public:: register_variable => register_variable_vector, register_variable_scalar
        !> @brief Associates pointer to a registered variable.
        procedure, public:: get_variable
    end type variable_register

    !> @brief Registers a particle's varaible, along with its derivative.
    type, extends(variable_register):: variable_deriv_register
        !> @brief Offset in memory of registered derivative variables.
        integer(c_intptr_t):: deriv_offsets(max_registrations)
    contains
        !> @brief Registers a scalar member of a particle and its derivative e.g. rho and drhodt.
        procedure, private:: register_variable_deriv_vector
        !> @brief Registers a vector member of a particle and its derivative e.g. v and dvxdt.
        procedure, private:: register_variable_deriv_scalar
        !> @brief Registers a member of particle and its derivative.
        generic, public:: register => register_variable_deriv_vector, register_variable_deriv_scalar
        !> @brief Associates pointers to a registered variable and its derivative.
        procedure, public:: get
    end type variable_deriv_register

    public:: max_registrations, variable_register, variable_deriv_register

contains

    !> @brief Register vector member of base and its derivative.
    !> @param self The register.
    !> @param base The particle who's member is being registered.
    !> @param name The name of the variable being registered.
    !> @param member The member variable of "base" being registered.
    !> @param member_deriv The derivative of the member being registered. Should also be a member of "base".
    subroutine register_variable_deriv_vector(self, base, name, member, member_deriv)
        class(variable_deriv_register), intent(inout):: self
        class(base_particle), target, intent(in):: base
        character(*), intent(in):: name
        real(fp), target, intent(in):: member(:), member_deriv(:)

        if (size(member) == 0) error stop "Cannot register 0-size member variable."
        if (size(member_deriv) == 0) error stop "Cannot register 0-size member derivative variable."
        if (size(member) /= size(member_deriv)) error stop "member and member_deriv are not same size."

        call register_variable_deriv_scalar(self, base, name, member(1), member_deriv(1))

        self%dims(self%nregistrations) = size(member)

    end subroutine register_variable_deriv_vector

    !> @brief Register scalar member of base and its derivative.
    !> @param self The register.
    !> @param base The particle who's member is being registered.
    !> @param name The name of the variable being registered.
    !> @param member The member variable of "base" being registered.
    !> @param member_deriv The derivative of the member being registered. Should also be a member of "base".
    subroutine register_variable_deriv_scalar(self, base, name, member, member_deriv)
        class(variable_deriv_register), intent(inout):: self
        class(base_particle), target, intent(in):: base
        character(*), intent(in):: name
        real(fp), target, intent(in):: member, member_deriv
        integer(c_intptr_t):: base_addr

        if (self%nregistrations == max_registrations) error stop "Exceeded maximum variable registrations."

        base_addr = transfer(c_loc(base%id), base_addr)

        self%nregistrations = self%nregistrations + 1
        self%names(self%nregistrations) = name
        self%offsets(self%nregistrations) = get_offset_(base_addr, member)
        self%deriv_offsets(self%nregistrations) = get_offset_(base_addr, member_deriv)
        self%dims(self%nregistrations) = 1

        if (self%offsets(self%nregistrations) >= sizeof(base)) error stop "Member is not a subset of base."
        if (self%deriv_offsets(self%nregistrations) >= sizeof(base)) error stop "member_deriv is not a subset of base."

    end subroutine register_variable_deriv_scalar

    !> @brief Calculates memory address offset of member, relative to base_addr.
    !> @param base_addr The address to calculate the offset from.
    !> @param member The member to calculate the address offset of.
    elemental integer(c_intptr_t) function get_offset_(base_addr, member)
        integer(c_intptr_t), intent(in):: base_addr
        real(fp), target, intent(in):: member

        get_offset_ = transfer(c_loc(member), base_addr) - base_addr

    end function get_offset_

    !> @brief Associates ptr and ptr_deriv to a registered variable and its derivative, respectively.
    !> @param self The register.
    !> @param base The particle who's registered variables to assign the pointers to.
    !> @param idx The index of the registered variable of interest.
    !> @param ptr The pointer that will be associated to base's member.
    !> @param ptr_deriv The pointer that will be associated to base's member's derivative.
    subroutine get(self, base, idx, ptr, ptr_deriv)
        class(variable_deriv_register), intent(in):: self
        class(base_particle), target, intent(in):: base
        integer, intent(in):: idx
        real(fp), pointer, intent(out):: ptr(:), ptr_deriv(:)
        integer(c_intptr_t):: base_addr

        base_addr = transfer(c_loc(base%id), base_addr)

        call ptr_from_offset_(base_addr, self%offsets(idx), self%dims(idx), ptr)
        call ptr_from_offset_(base_addr, self%deriv_offsets(idx), self%dims(idx), ptr_deriv)

    end subroutine get

    !> @brief Associates "ptr" to the variable located at base_addr + offset.
    !> @brief base_addr The base address that offset is calculated from.
    !> @brief offset The memory address offset that the variable is located at.
    !> @brief dims The dimension of the variable (1 for scalar, >1 for vector).
    !> @brief ptr The pointer to be assigned to the variable located at base_addr + offset.
    subroutine ptr_from_offset_(base_addr, offset, dims, ptr)
        integer(c_intptr_t), intent(in):: base_addr, offset
        integer, intent(in):: dims
        real(fp), pointer, intent(out):: ptr(:)
        integer(c_intptr_t):: member_offset
        type(c_ptr):: member_c_ptr

        ! calculate address of member
        member_offset = base_addr + offset
        ! convert address to c_ptr
        member_c_ptr = transfer(member_offset, member_c_ptr)
        ! convert c_ptr to Fortran pointer
        call c_f_pointer(member_c_ptr, ptr, [dims])

    end subroutine ptr_from_offset_

    !> @brief Register vector member of base.
    !> @param self The register.
    !> @param base The particle who's member is being registered.
    !> @param name The name of the variable being registered.
    !> @param member The member variable of "base" being registered.
    subroutine register_variable_vector(self, base, name, member)
        class(variable_register), intent(inout):: self
        class(base_particle), target, intent(in):: base
        character(*), intent(in):: name
        real(fp), target, intent(in):: member(:)

        if (size(member) == 0) error stop "Cannot register 0-size member variable."

        call register_variable_scalar(self, base, name, member(1))

        self%dims(self%nregistrations) = size(member)

    end subroutine register_variable_vector

    !> @brief Register scalar member of base.
    !> @param self The register.
    !> @param base The particle who's member is being registered.
    !> @param name The name of the variable being registered.
    !> @param member The member variable of "base" being registered.
    subroutine register_variable_scalar(self, base, name, member)
        class(variable_register), intent(inout):: self
        class(base_particle), target, intent(in):: base
        character(*), intent(in):: name
        real(fp), target, intent(in):: member
        integer(c_intptr_t):: base_addr

        if (self%nregistrations == max_registrations) error stop "Exceeded maximum variable registrations."

        base_addr = transfer(c_loc(base%id), base_addr)

        self%nregistrations = self%nregistrations + 1
        self%names(self%nregistrations) = name
        self%offsets(self%nregistrations) = get_offset_(base_addr, member)
        self%dims(self%nregistrations) = 1

        if (self%offsets(self%nregistrations) >= sizeof(base)) error stop "Member is not a subset of base."

    end subroutine register_variable_scalar

    !> @brief Associates ptr to a registered variable.
    !> @param self The register.
    !> @param base The particle who's registered variables to assign the pointers to.
    !> @param idx The index of the registered variable of interest.
    !> @param ptr The pointer that will be associated to base's member.
    subroutine get_variable(self, base, idx, ptr)
        class(variable_register), intent(in):: self
        class(base_particle), target, intent(in):: base
        integer, intent(in):: idx
        real(fp), pointer, intent(out):: ptr(:)
        integer(c_intptr_t):: base_addr

        base_addr = transfer(c_loc(base%id), base_addr)

        call ptr_from_offset_(base_addr, self%offsets(idx), self%dims(idx), ptr)

    end subroutine get_variable

end module grasph_register
