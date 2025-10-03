!> @file grasph_constants.F90
!> @brief Holds constants used throughout GraSPH.
!> @author Edward Yang
!> @date 2025-06-01
module grasph_constants_m

    use iso_fortran_env, only: real64, real32

    implicit none

    public

#ifdef SINGLE_PRECISION
    integer, parameter:: fp = real32
#else
    integer, parameter:: fp = real64
#endif

    real(fp), parameter:: pi = 4._fp*atan(1._fp)

#ifdef THREED
    integer, parameter:: ndims = 3
#else
    integer, parameter:: ndims = 2
#endif

end module grasph_constants_m
