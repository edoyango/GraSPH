!> @file grasph_common.f90
!> @brief Contains miscelaneous utilities used across the package.
!> @author Edward Yang
!> @date 2025-09-26
module grasph_common_m

    use grasph_constants_m, only: fp

    implicit none

    public

    !> @brief A type that contina a pointer to a 2d array.
    type array_pointer_container_t
        !> @brief The 2d array pointer.
        real(fp), pointer:: p(:, :)
    end type array_pointer_container_t

end module grasph_common_m
