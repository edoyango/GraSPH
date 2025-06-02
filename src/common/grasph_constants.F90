module grasph_constants

    use iso_fortran_env, only: real64, real32

    implicit none

    public

#ifdef SINGLE_PRECISION
    integer, parameter:: fp = real32
#else
    integer, parameter:: fp = real64
#endif

    real(fp), parameter:: pi = 4._fp*atan(1._fp)

    integer, parameter:: ndims = 3

end module grasph_constants