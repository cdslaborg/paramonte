program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_optimization, only: setMinBrent
    use pm_optimization, only: setBracketMin

    implicit none

    integer(IK) :: niter, retin
    real(RKC) :: xlow, xmin, xupp, fmin, tol
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getSq(x) = (x - 1)**2")
    call disp%show("xlow = -3; xupp = -1; tol = epsilon(xmin)**.8; retin = 100; niter = 100")
                    xlow = -3; xupp = -1; tol = epsilon(xmin)**.8; retin = 100; niter = 100
    call disp%show("call setBracketMin(getSq, retin, xmin, xlow, xupp, fmin) ! find a good bracket, though here the choice is obvious.")
                    call setBracketMin(getSq, retin, xmin, xlow, xupp, fmin) ! find a good bracket, though here the choice is obvious.
    call disp%show("if (100 < retin) error stop 'Bracketing failed.'")
                    if (100 < retin) error stop 'Bracketing failed.'
    call disp%show("call setMinBrent(getSq, xmin, xlow, xupp, fmin, tol, niter)")
                    call setMinBrent(getSq, xmin, xlow, xupp, fmin, tol, niter)
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("if (niter > 100) error stop 'minimization failed.'")
                    if (niter > 100) error stop 'minimization failed.'
    call disp%show("[xmin, fmin]")
    call disp%show( [xmin, fmin] )
    call disp%skip()

    call disp%skip()
    call disp%show("getSq(x) = (x - 1)**2")
    call disp%show("xlow = 1; xmin = 1; xupp = 3; tol = epsilon(xmin)**.8; niter = 100")
                    xlow = 1; xmin = 1; xupp = 3; tol = epsilon(xmin)**.8; niter = 100
    call disp%show("call setMinBrent(getSq, xmin, xlow, xupp, fmin, tol, niter)")
                    call setMinBrent(getSq, xmin, xlow, xupp, fmin, tol, niter)
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("if (niter > 100) error stop 'minimization failed.'")
                    if (niter > 100) error stop 'minimization failed.'
    call disp%show("[xmin, fmin]")
    call disp%show( [xmin, fmin] )
    call disp%skip()

    call disp%skip()
    call disp%show("getSq(x) = (x - 1)**2")
    call disp%show("xlow = 1; xmin = 1; xupp = 1; tol = epsilon(xmin)**.8; niter = 100")
                    xlow = 1; xmin = 1; xupp = 1; tol = epsilon(xmin)**.8; niter = 100
    call disp%show("call setMinBrent(getSq, xmin, xlow, xupp, fmin, tol, niter)")
                    call setMinBrent(getSq, xmin, xlow, xupp, fmin, tol, niter)
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("if (niter > 100) error stop 'minimization failed.'")
                    if (niter > 100) error stop 'minimization failed.'
    call disp%show("[xmin, fmin]")
    call disp%show( [xmin, fmin] )
    call disp%skip()

    call disp%skip()
    call disp%show("getSq(x) = (x - 1)**2")
    call disp%show("xmin = 0; xlow = -1; xupp = 3; fmin = getSq(xmin); tol = epsilon(xmin)**.8; niter = 100")
                    xmin = 0; xlow = -1; xupp = 3; fmin = getSq(xmin); tol = epsilon(xmin)**.8; niter = 100
    call disp%show("call setMinBrent(getSq, xmin, xlow, xupp, fmin, tol, niter)")
                    call setMinBrent(getSq, xmin, xlow, xupp, fmin, tol, niter)
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("if (niter > 100) error stop 'minimization failed.'")
                    if (niter > 100) error stop 'minimization failed.'
    call disp%show("[xmin, fmin]")
    call disp%show( [xmin, fmin] )
    call disp%skip()

contains

    function getSq(x) result(func)
        real(RKC), intent(in) :: x
        real(RKC) :: func
        func = (x - 1)**2
    end function

end program example