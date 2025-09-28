!> @file grasph_common.f90
!> @brief Contains miscelaneous utilities used across the package.
!> @author Edward Yang
!> @date 2025-09-26
module grasph_common

    use grasph_constants, only: fp

    implicit none

    public

    type array_pointer_container
        real(fp), pointer:: p(:, :)
    end type array_pointer_container

end module grasph_common
