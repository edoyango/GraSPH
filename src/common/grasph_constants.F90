!> @file grasph_constants.F90
!> @brief Holds constants used throughout GraSPH.
!> @author Edward Yang
!> @date 2025-06-01
module grasph_constants_m

    use iso_fortran_env, only: real64, real32

    implicit none

    public

#ifdef SINGLE_PRECISION
    !> @brief The kind type parameter to use for floats (4-byte)
    integer, parameter:: fp = real32
#else
    !> @brief The kind type parameter to use for floats (8-byte)
    integer, parameter:: fp = real64
#endif

    !> @brief pi
    real(fp), parameter:: pi = 4._fp*atan(1._fp)

#ifdef THREED
    !> @brief The number of dimensions of the problem (3)
    integer, parameter:: ndims = 3
#else
    !> @brief The number of dimensions of the problem (2)
    integer, parameter:: ndims = 2
#endif

    !> @brief Maximum allowable characters for a variable name.
    integer, parameter:: max_name_len = 30

end module grasph_constants_m
