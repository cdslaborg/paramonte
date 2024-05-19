program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_optimization, only: getMinBrent

    implicit none

    integer(IK) :: niter
    real(RKG) :: xmin, fmin
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getSq(x) = (x - 1)**2")
    call disp%show("xmin = getMinBrent(getSq)")
                    xmin = getMinBrent(getSq)
    call disp%show("xmin")
    call disp%show( xmin )
    call disp%skip()

    call disp%skip()
    call disp%show("getSq(x) = (x - 1)**2")
    call disp%show("niter = 10")
                    niter = 10
    call disp%show("xmin = getMinBrent(getSq, fmin = fmin, tol = epsilon(xmin)**.8, niter = niter)")
                    xmin = getMinBrent(getSq, fmin = fmin, tol = epsilon(xmin)**.8, niter = niter)
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("if (niter > 10) error stop 'minimization failed.'")
                    if (niter > 10) error stop 'minimization failed.'
    call disp%show("[xmin, fmin]")
    call disp%show( [xmin, fmin] )
    call disp%skip()

contains

    function getSq(x) result(func)
        real(RKG), intent(in) :: x
        real(RKG) :: func
        func = (x - 1)**2
    end function

end program example