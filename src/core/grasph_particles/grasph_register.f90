!> @file grasph_register.f90
!> @brief Module containing classes for registering particle variables e.g. for io or time-evolution.
!> @author Edward Yang
!> @date 2025-10-02
module grasph_register_m

    use grasph_constants_m, only: fp, ndims
    use grasph_particle_m, only: base_particles_t
    use grasph_common_m, only: array_pointer_container_t

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
        !> @brief Name of variables registered. Used in IO.
        character(max_variable_name), allocatable:: names(:)
        !> @brief The variables that have been registered.
        type(array_pointer_container_t), allocatable:: variables(:)
    contains
        !> @brief Registers a scalar particle variable array e.g. density.
        procedure, private:: register_variable_vector
        !> @brief Registers a vector particle variable array e.g. position.
        procedure, private:: register_variable_scalar
        !> @brief Registers a particle variable array.
        generic, public:: register_variable => register_variable_vector, register_variable_scalar
        !> @brief Utility to automatically resize internal arrays and return a "safe" last index.
        procedure, private:: safe_size_plus_1 => variable_register_safe_size_plus_1
        !> @brief Given a variable name, deregisters that variable from the register.
        procedure, public:: deregister => deregister_variable
    end type variable_register_t

    !> @brief Registers a particle's varaible, along with its derivative.
    type, extends(variable_register_t):: variable_deriv_register_t
        !> @brief Offset in memory of registered derivative variables.
        type(array_pointer_container_t), allocatable:: derivatives(:)
    contains
        !> @brief Registers a particle scalar variable array and its derivative e.g. rho and drhodt.
        procedure, private:: register_variable_deriv_vector
        !> @brief Registers a particle vector variable array and its derivative e.g. v and dvxdt.
        procedure, private:: register_variable_deriv_scalar
        !> @brief Registers a particle variable array and its derivative.
        generic, public:: register => register_variable_deriv_vector, register_variable_deriv_scalar
        !> @brief Utility to automatically resize internal arrays and return a "safe" last index.
        procedure, private:: safe_size_plus_1 => deriv_register_safe_size_plus_1
        !> @brief Given a variable name, deregisters that variable and its derivative from the register.
        procedure, public:: deregister => deregister_variable_deriv
    end type variable_deriv_register_t

    public:: variable_register_t, variable_deriv_register_t

contains

    !> @brief Register vector variable and its derivative.
    !> @param self The register.
    !> @param name The name of the variable being registered.
    !> @param variable The variable being registered.
    !> @param derivative The derivative of the variable being registered.
    subroutine register_variable_deriv_vector(self, name, variable, derivative)
        class(variable_deriv_register_t), intent(inout):: self
        character(*), intent(in):: name
        real(fp), target, intent(in):: variable(:, :), derivative(:, :)
        integer:: nregs

        if (size(variable) == 0) error stop "Cannot register 0-size variable."
        if (size(derivative) == 0) error stop "Cannot register 0-size derivative."
        if (size(variable) /= size(derivative)) error stop "variable and derivative are not same size."

        nregs = self%safe_size_plus_1()
        self%names(nregs) = name
        self%dims(nregs) = size(variable, dim=1)
        self%variables(nregs)%p => variable
        self%derivatives(nregs)%p => derivative

    end subroutine register_variable_deriv_vector

    !> @brief Register scalar variable and its derivative.
    !> @param self The register.
    !> @param name The name of the variable being registered.
    !> @param variable The variable being registered.
    !> @param derivative The derivative of the variable being registered.
    subroutine register_variable_deriv_scalar(self, name, variable, derivative)
        class(variable_deriv_register_t), intent(inout):: self
        character(*), intent(in):: name
        real(fp), target, intent(in):: variable(:), derivative(:)
        integer:: nregs

        if (size(variable) == 0) error stop "Cannot register 0-size variable."
        if (size(derivative) == 0) error stop "Cannot register 0-size derivative."
        if (size(variable) /= size(derivative)) error stop "variable and derivative are not same size."

        nregs = self%safe_size_plus_1()
        self%names(nregs) = name
        self%variables(nregs)%p(1:size(variable), 1:1) => variable
        self%derivatives(nregs)%p(1:size(derivative), 1:1) => derivative
        self%dims(nregs) = 1

    end subroutine register_variable_deriv_scalar

    !> @brief Register vector variable.
    !> @param self The register.
    !> @param name The name of the variable being registered.
    !> @param variable The variable being registered.
    subroutine register_variable_vector(self, name, variable)
        class(variable_register_t), intent(inout):: self
        character(*), intent(in):: name
        real(fp), target, intent(in):: variable(:, :)
        integer:: nregs

        if (size(variable) == 0) error stop "Cannot register 0-size variable."

        nregs = self%safe_size_plus_1()
        self%names(nregs) = name
        self%dims(nregs) = size(variable, dim=1)
        self%variables(nregs)%p => variable
        self%nregistrations = nregs

    end subroutine register_variable_vector

    !> @brief Register scalar variable.
    !> @param self The register.
    !> @param name The name of the variable being registered.
    !> @param variable The variable being registered.
    subroutine register_variable_scalar(self, name, variable)
        class(variable_register_t), intent(inout):: self
        character(*), intent(in):: name
        real(fp), target, intent(in):: variable(:)
        integer:: nregs

        if (size(variable) == 0) error stop "Cannot register 0-size variable."

        nregs = self%safe_size_plus_1()
        self%names(nregs) = name
        self%variables(nregs)%p(1:size(variable), 1:1) => variable
        self%dims(nregs) = 1
        self%nregistrations = nregs

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
        self%variables(1:n_new_regs) = pack(self%variables(1:n_old_regs), match(:))
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
        self%variables(1:n_new_regs) = pack(self%variables(1:n_old_regs), match(:))
        self%names(1:n_new_regs) = pack(self%names(1:n_old_regs), match(:))
        self%derivatives(1:n_new_regs) = pack(self%derivatives(1:n_old_regs), match(:))
        self%nregistrations = n_new_regs
    end subroutine deregister_variable_deriv

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

    subroutine realloc_ptr_container_t_r1_(arr, n)
        type(array_pointer_container_t), allocatable, intent(inout):: arr(:)
        integer, intent(in):: n
        type(array_pointer_container_t), allocatable:: tmp_arr(:)

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
    end subroutine realloc_ptr_container_t_r1_

    function variable_register_safe_size_plus_1(self) result(new_size)
        class(variable_register_t), intent(inout):: self
        integer:: new_size

        call realloc_integer_r1_(self%dims, self%nregistrations)
        call realloc_character_r1_(self%names, self%nregistrations)
        call realloc_ptr_container_t_r1_(self%variables, self%nregistrations)

        self%nregistrations = self%nregistrations + 1
        new_size = self%nregistrations

    end function variable_register_safe_size_plus_1

    function deriv_register_safe_size_plus_1(self) result(new_size)
        class(variable_deriv_register_t), intent(inout):: self
        integer:: new_size

        call realloc_integer_r1_(self%dims, self%nregistrations)
        call realloc_character_r1_(self%names, self%nregistrations)
        call realloc_ptr_container_t_r1_(self%variables, self%nregistrations)
        call realloc_ptr_container_t_r1_(self%derivatives, self%nregistrations)

        self%nregistrations = self%nregistrations + 1
        new_size = self%nregistrations

    end function deriv_register_safe_size_plus_1

end module grasph_register_m
