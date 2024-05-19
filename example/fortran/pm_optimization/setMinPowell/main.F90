program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_optimization, only: setMinPowell
    use pm_arrayFill, only: getFilled
    use pm_matrixInit, only: getMatInit, uppLowDia

    implicit none

    integer(IK) :: niter
    integer(IK), parameter :: ndim = 4
    real(RKG) :: xmin(ndim), dirset(ndim, ndim), fmin
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getSq(x) = (x - [real(RKG) :: 1, 2, 3, 4])**2")
    call disp%show("niter = 100; xmin = getFilled(0., ndim); fmin = getSq(xmin); dirset = getMatInit([ndim, ndim], uppLowDia, 0._RKG, 0._RKG, 1._RKG)")
                    niter = 100; xmin = getFilled(0., ndim); fmin = getSq(xmin); dirset = getMatInit([ndim, ndim], uppLowDia, 0._RKG, 0._RKG, 1._RKG)
    call disp%show("call setMinPowell(getSq, xmin, fmin, dirset, epsilon(xmin)**.5, niter)")
                    call setMinPowell(getSq, xmin, fmin, dirset, epsilon(xmin)**.5, niter)
    call disp%show("if (100 < niter) error stop 'minimization failed.'")
                    if (100 < niter) error stop 'minimization failed.'
    call disp%show("[xmin, fmin]")
    call disp%show( [xmin, fmin] )
    call disp%skip()

contains

    function getSq(x) result(func)
        real(RKG)   , intent(in) :: x(ndim)
        real(RKG)   :: func
        integer(IK) :: idim
        func = sum((x - [real(RKG) :: (idim, idim = 1, ndim)])**2)
    end function

end program example