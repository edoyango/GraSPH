!> @file grasph_kernels.f90
!> @brief Contains truncated kernels that can be used in SPH.
!> @author Edward Yang
!> @date 2025-06-01
module grasph_kernels_m

    use grasph_constants_m, only: fp, pi

    implicit none

    private

    !> @brief base kernel which describes required members and methods of a kernel
    type, abstract:: base_kernel_t
        !> @brief Number of spatial dimensions
        integer:: d = 0
        !> @brief Normalization factor base on kernel and spatial dimension.
        real(fp):: alpha = 0._fp
        !> @brief Smoothing length to be used in calculations (will need to be changed in the future when particles have variable
        !>        smoothing lengths...)
        real(fp):: h = 0._fp
        !> @brief The cutoff distance of the kernel i.e., the value of at and beyond which the kernel value is 0.
        real(fp):: cutoff = 0._fp
        !> @brief Whether the kernel has been initialized.
        logical:: initialized = .false.
    contains
        !> @brief The initalizer intended to set d, alpha, h, cutoff, initialized.
        procedure(kernel_init_interface), deferred:: init
        !> @brief A subroutine which calculates kernel value and gradients in one call.
        procedure:: values
        !> @brief The kernel function.
        procedure(w_interface), deferred:: w
        !> @brief The spatial gradient of the kernel gradient.
        procedure(w_interface), deferred:: gradw
    end type base_kernel_t

    interface
        !> @brief The kernel function interface.
        !> @param self The kernel class with required constants.
        !> @param q The distance value, normalized by smoothing length.
        real(fp) pure function w_interface(self, q)
            import:: fp, base_kernel_t
            class(base_kernel_t), intent(in):: self
            real(fp), intent(in):: q
        end function w_interface
        !> @brief The kernel initializer interface.
        !> @param self The kernel to be initialized.
        !> @param d The spatial dimension.
        !> @param h The constant smoothing length.
        subroutine kernel_init_interface(self, d, h)
            import:: fp, base_kernel_t
            class(base_kernel_t), intent(inout):: self
            integer, intent(in):: d
            real(fp), intent(in):: h
        end subroutine kernel_init_interface
    end interface

    !> @brief Implementation of the cubic B-spline SPH kernel.
    type, extends(base_kernel_t):: cubic_bspline_kernel_t
    contains
        !> @brief Overriden initializer for the cubic B-spline kernel.
        procedure:: init => init_cubic_bspline
        !> @brief Concrete definition of the cubic B-spline kernel function.
        procedure:: w => w_cubic_bspline
        !> @brief Concret definition of the cubic B-spline kernel gradient.
        procedure:: gradw => gradw_cubic_bspline
    end type cubic_bspline_kernel_t

    public:: base_kernel_t, cubic_bspline_kernel_t
contains

    !> @brief A call to calculate both kernel values and gradients.
    !> @param self The kernel with the kernel function and gradients implementation to use.
    !> @param dx The relative position vector between particles.
    !> @param w The return kernel value.
    !> @param dwdx The return kernel gradient values.
    pure subroutine values(self, dx, w, dwdx)
        class(base_kernel_t), intent(in):: self
        real(fp), intent(in):: dx(self%d)
        real(fp), intent(out):: w, dwdx(self%d)
        real(fp):: q, r
        r = sqrt(sum(dx**2))
        q = r/self%h
        w = self%w(q)
        dwdx(:) = self%gradw(q)*dx(:)/(r*self%h)
    end subroutine values

    !> @brief The implementation of the kernel initializer for cubic B-spline.
    !> @param self The cubic B-spline kernel to initialize.
    !> @param d The spatial dimension.
    !> @param h The smoothing length to use.
    subroutine init_cubic_bspline(self, d, h)
        class(cubic_bspline_kernel_t), intent(inout):: self
        integer, intent(in):: d
        real(fp), intent(in):: h
        self%d = d
        self%h = h
        self%cutoff = 2._fp*h
        if (d == 2) then
            self%alpha = 10._fp/(7._fp*pi*h*h)
        elseif (d == 3) then
            self%alpha = 1._fp/(pi*h*h*h)
        end if
        self%initialized = .true.
    end subroutine init_cubic_bspline

    !> @brief The implementation of the kernel function for cubic B-spline.
    !> @param self The cubic B-spline kernel with needed constants.
    !> @param q The normalized distance between particles.
    real(fp) pure function w_cubic_bspline(self, q)
        class(cubic_bspline_kernel_t), intent(in):: self
        real(fp), intent(in):: q
        w_cubic_bspline = self%alpha*(0.25_fp*dim(2._fp, q)**3 - dim(1._fp, q)**3)
    end function w_cubic_bspline

    !> @brief The implementation of the kernel gradient function for cubic B-spline.
    !> @param self The cubic B-spline kernel with needed constants.
    !> @param q The normalized distance between particles.
    real(fp) pure function gradw_cubic_bspline(self, q)
        class(cubic_bspline_kernel_t), intent(in):: self
        real(fp), intent(in):: q
        gradw_cubic_bspline = -self%alpha*3._fp*(0.25_fp*dim(2._fp, q)**2 - dim(1._fp, q)**2)
    end function gradw_cubic_bspline

end module grasph_kernels_m
