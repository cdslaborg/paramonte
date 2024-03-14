program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_optimization, only: isFailedMinPowell
    use pm_arrayFill, only: getFilled

    implicit none

    integer(IK), parameter :: ndim = 4
    real(RKC) :: xmin(ndim), fmin
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getSq(x) = sum(x - [real(RKC) :: 1, 2, 3, 4])**2)")
    call disp%show("xmin = getFilled(0., ndim)")
                    xmin = getFilled(0., ndim)
    call disp%show("if (isFailedMinPowell(getSq, xmin)) error stop 'minimization failed.'")
                    if (isFailedMinPowell(getSq, xmin)) error stop 'minimization failed.'
    call disp%show("xmin")
    call disp%show( xmin )
    call disp%skip()

    call disp%skip()
    call disp%show("getSq(x) = (x - [real(RKC) :: 1, 2, 3, 4])**2")
    call disp%show("xmin = getFilled(0., ndim); fmin = getSq(xmin)")
                    xmin = getFilled(0., ndim); fmin = getSq(xmin)
    call disp%show("if (isFailedMinPowell(getSq, xmin, fmin, tol = epsilon(xmin)**.8, niter = 100)) error stop 'minimization failed.'")
                    if (isFailedMinPowell(getSq, xmin, fmin, tol = epsilon(xmin)**.8, niter = 100)) error stop 'minimization failed.'
    call disp%show("[xmin, fmin]")
    call disp%show( [xmin, fmin] )
    call disp%skip()

contains

    function getSq(x) result(func)
        real(RKC)   , intent(in) :: x(ndim)
        real(RKC)   :: func
        integer(IK) :: idim
        func = sum((x - [real(RKC) :: (idim, idim = 1, ndim)])**2)
    end function

end program example