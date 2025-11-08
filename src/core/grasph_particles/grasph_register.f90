!> @file grasph_register.f90
!> @brief Module containing classes for registering particle variables e.g. for io or time-evolution.
!> @author Edward Yang
!> @date 2025-10-02
module grasph_register_m

    use grasph_constants_m, only: fp, ndims
    use grasph_particle_m, only: base_particle_t
    use iso_c_binding, only: c_intptr_t, c_f_pointer, c_ptr, c_loc

    implicit none

    private

    !> @brief Max allowable registrations. Will be removed for a more dynamic approach.
    integer, parameter:: max_variable_name = 20

    !> @brief Registers particle's variables for access through pointer e.g. for dynamically selecting variables for IO.
    type variable_register_t
        !> @brief Number of variables currently registered.
        integer:: nregistrations = 0
        !> @brief Dimension of each variable registered.
        integer, allocatable:: dims(:)
        !> @brief Name of variables registere. Used in IO.
        character(max_variable_name), allocatable:: names(:)
        !> @brief Offset in memory of registered variables.
        integer(c_intptr_t), allocatable:: offsets(:)
    contains
        !> @brief Registers a scalar member of a particle e.g. density.
        procedure, private:: register_variable_vector
        !> @brief Registers a vector member of a particle e.g. position.
        procedure, private:: register_variable_scalar
        !> @brief Registers a member of a particle.
        generic, public:: register_variable => register_variable_vector, register_variable_scalar
        !> @brief Associates pointer to a registered variable.
        procedure, public:: get_variable
        !> @brief Utility to automatically resize internal arrays and return a "safe" last index.
        procedure, private:: safe_size_plus_1 => variable_register_safe_size_plus_1
        !> @brief Given a variable name, deregisters that variable from the register.
        procedure, public:: deregister => deregister_variable
    end type variable_register_t

    !> @brief Registers a particle's varaible, along with its derivative.
    type, extends(variable_register_t):: variable_deriv_register_t
        !> @brief Offset in memory of registered derivative variables.
        integer(c_intptr_t), allocatable:: deriv_offsets(:)
    contains
        !> @brief Registers a scalar member of a particle and its derivative e.g. rho and drhodt.
        procedure, private:: register_variable_deriv_vector
        !> @brief Registers a vector member of a particle and its derivative e.g. v and dvxdt.
        procedure, private:: register_variable_deriv_scalar
        !> @brief Registers a member of particle and its derivative.
        generic, public:: register => register_variable_deriv_vector, register_variable_deriv_scalar
        !> @brief Associates pointers to a registered variable and its derivative.
        procedure, public:: get
        !> @brief Utility to automatically resize internal arrays and return a "safe" last index.
        procedure, private:: safe_size_plus_1 => deriv_register_safe_size_plus_1
        !> @brief Given a variable name, deregisters that variable and its derivative from the register.
        procedure, public:: deregister => deregister_variable_deriv
    end type variable_deriv_register_t

    public:: variable_register_t, variable_deriv_register_t

contains

    !> @brief Register vector member of base and its derivative.
    !> @param self The register.
    !> @param base The particle who's member is being registered.
    !> @param name The name of the variable being registered.
    !> @param member The member variable of "base" being registered.
    !> @param member_deriv The derivative of the member being registered. Should also be a member of "base".
    subroutine register_variable_deriv_vector(self, base, name, member, member_deriv)
        class(variable_deriv_register_t), intent(inout):: self
        class(base_particle_t), target, intent(in):: base
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
        class(variable_deriv_register_t), intent(inout):: self
        class(base_particle_t), target, intent(in):: base
        character(*), intent(in):: name
        real(fp), target, intent(in):: member, member_deriv
        integer(c_intptr_t):: base_addr
        integer:: nregs

        base_addr = transfer(c_loc(base%id), base_addr)

        nregs = self%safe_size_plus_1()
        self%names(nregs) = name
        self%offsets(nregs) = get_offset_(base_addr, member)
        self%deriv_offsets(nregs) = get_offset_(base_addr, member_deriv)
        self%dims(nregs) = 1

        if (self%offsets(nregs) >= sizeof(base)) error stop "Member is not a subset of base."
        if (self%deriv_offsets(nregs) >= sizeof(base)) error stop "member_deriv is not a subset of base."

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
        class(variable_deriv_register_t), intent(in):: self
        class(base_particle_t), target, intent(in):: base
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
        class(variable_register_t), intent(inout):: self
        class(base_particle_t), target, intent(in):: base
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
        class(variable_register_t), intent(inout):: self
        class(base_particle_t), target, intent(in):: base
        character(*), intent(in):: name
        real(fp), target, intent(in):: member
        integer(c_intptr_t):: base_addr
        integer:: nregs

        base_addr = transfer(c_loc(base%id), base_addr)

        nregs = self%safe_size_plus_1()
        self%names(nregs) = name
        self%offsets(nregs) = get_offset_(base_addr, member)
        self%dims(nregs) = 1

        if (self%offsets(nregs) >= sizeof(base)) error stop "Member is not a subset of base."

    end subroutine register_variable_scalar

    !> @brief Given a variable name, deregister that variable from the register.
    !> @param self The register to deregister the variable from.
    !> @param name The name of the variable to deregister.
    subroutine deregister_variable(self, name)
        class(variable_register_t), intent(inout):: self
        character(*), intent(in):: name
        logical:: match(self%nregistrations)
        integer:: n_old_regs, n_new_regs

        n_old_regs = self%nregistrations
        match(:) = self%names(1:n_old_regs) /= name

        n_new_regs = count(match)

        self%dims(1:n_new_regs) = pack(self%dims(1:n_old_regs), match(:))
        self%offsets(1:n_new_regs) = pack(self%offsets(1:n_old_regs), match(:))
        self%names(1:n_new_regs) = pack(self%names(1:n_old_regs), match(:))
        self%nregistrations = n_new_regs
    end subroutine deregister_variable

    !> @brief Given a variable name, deregister that variable and its derivative from the register.
    !> @param self The register to deregister the variable and derivative from.
    !> @param name The name of the variable to deregister.
    subroutine deregister_variable_deriv(self, name)
        class(variable_deriv_register_t), intent(inout):: self
        character(*), intent(in):: name
        logical:: match(self%nregistrations)
        integer:: n_old_regs, n_new_regs

        n_old_regs = self%nregistrations
        match(:) = self%names(1:n_old_regs) /= name

        n_new_regs = count(match)

        self%dims(1:n_new_regs) = pack(self%dims(1:n_old_regs), match(:))
        self%offsets(1:n_new_regs) = pack(self%offsets(1:n_old_regs), match(:))
        self%names(1:n_new_regs) = pack(self%names(1:n_old_regs), match(:))
        self%deriv_offsets(1:n_new_regs) = pack(self%deriv_offsets(1:n_old_regs), match(:))
        self%nregistrations = n_new_regs
    end subroutine deregister_variable_deriv

    !> @brief Associates ptr to a registered variable.
    !> @param self The register.
    !> @param base The particle who's registered variables to assign the pointers to.
    !> @param idx The index of the registered variable of interest.
    !> @param ptr The pointer that will be associated to base's member.
    subroutine get_variable(self, base, idx, ptr)
        class(variable_register_t), intent(in):: self
        class(base_particle_t), target, intent(in):: base
        integer, intent(in):: idx
        real(fp), pointer, intent(out):: ptr(:)
        integer(c_intptr_t):: base_addr

        base_addr = transfer(c_loc(base%id), base_addr)

        call ptr_from_offset_(base_addr, self%offsets(idx), self%dims(idx), ptr)

    end subroutine get_variable

    subroutine realloc_integer_r1_(arr, n)
        integer, allocatable, intent(inout):: arr(:)
        integer, intent(in):: n
        integer, allocatable:: tmp_arr(:)

        if (.not. allocated(arr)) then
            allocate (arr(1))
        else if (size(arr) == 0) then
            deallocate (arr)
            allocate (arr(1))
        else if (size(arr) == n) then
            call move_alloc(arr, tmp_arr)
            allocate (arr(2*size(tmp_arr)))
            arr(1:n) = tmp_arr(1:n)
            deallocate (tmp_arr)
        end if
    end subroutine realloc_integer_r1_

    subroutine realloc_character_r1_(arr, n)
        character(max_variable_name), allocatable, intent(inout):: arr(:)
        integer, intent(in):: n
        character(max_variable_name), allocatable:: tmp_arr(:)

        if (.not. allocated(arr)) then
            allocate (arr(1))
        else if (size(arr) == 0) then
            deallocate (arr)
            allocate (arr(1))
        else if (size(arr) == n) then
            call move_alloc(arr, tmp_arr)
            allocate (arr(2*size(tmp_arr)))
            arr(1:n) = tmp_arr(1:n)
            deallocate (tmp_arr)
        end if
    end subroutine realloc_character_r1_

    subroutine realloc_c_intptr_t_r1_(arr, n)
        integer(c_intptr_t), allocatable, intent(inout):: arr(:)
        integer, intent(in):: n
        integer(c_intptr_t), allocatable:: tmp_arr(:)

        if (.not. allocated(arr)) then
            allocate (arr(1))
        else if (size(arr) == 0) then
            deallocate (arr)
            allocate (arr(1))
        else if (size(arr) == n) then
            call move_alloc(arr, tmp_arr)
            allocate (arr(2*size(tmp_arr)))
            arr(1:n) = tmp_arr(1:n)
            deallocate (tmp_arr)
        end if
    end subroutine realloc_c_intptr_t_r1_

    function variable_register_safe_size_plus_1(self) result(new_size)
        class(variable_register_t), intent(inout):: self
        integer:: new_size

        call realloc_integer_r1_(self%dims, self%nregistrations)
        call realloc_character_r1_(self%names, self%nregistrations)
        call realloc_c_intptr_t_r1_(self%offsets, self%nregistrations)

        self%nregistrations = self%nregistrations + 1
        new_size = self%nregistrations

    end function variable_register_safe_size_plus_1

    function deriv_register_safe_size_plus_1(self) result(new_size)
        class(variable_deriv_register_t), intent(inout):: self
        integer:: new_size

        call realloc_integer_r1_(self%dims, self%nregistrations)
        call realloc_character_r1_(self%names, self%nregistrations)
        call realloc_c_intptr_t_r1_(self%offsets, self%nregistrations)
        call realloc_c_intptr_t_r1_(self%deriv_offsets, self%nregistrations)

        self%nregistrations = self%nregistrations + 1
        new_size = self%nregistrations

    end function deriv_register_safe_size_plus_1

end module grasph_register_m
