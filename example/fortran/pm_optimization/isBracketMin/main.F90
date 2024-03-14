program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_optimization, only: setBracketMin, isBracketMin
    use pm_err, only: setAsserted
    use pm_val2str, only: getStr

    implicit none

    integer(IK) :: niter, minniter = 1000
    real(RKC) :: xlow, xmin, xupp, flow, fmin, fupp
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getSqNeg(x) = (x - 1)**2")
    call disp%show("xlow = -3; xupp = -1; niter = minniter")
                    xlow = -3; xupp = -1; niter = minniter
    call disp%show("call setBracketMin(getSqNeg, niter, xmin, xlow, xupp, fmin)")
                    call setBracketMin(getSqNeg, niter, xmin, xlow, xupp, fmin)
    call disp%show("if (niter > minniter) error stop 'Bracketing failed.'")
                    if (niter > minniter) error stop 'Bracketing failed.'
    call disp%show("[xlow, xmin, xupp, getSqNeg(xlow), fmin, getSqNeg(xupp)]")
    call disp%show( [xlow, xmin, xupp, getSqNeg(xlow), fmin, getSqNeg(xupp)] )
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("call setAsserted(isBracketMin(xmin, xlow, xupp, fmin, getSqNeg(xlow), getSqNeg(xupp)))")
                    call setAsserted(isBracketMin(xmin, xlow, xupp, fmin, getSqNeg(xlow), getSqNeg(xupp)))
    call disp%skip()

    call disp%skip()
    call disp%show("getSqNeg(x) = (x - 1)**2")
    call disp%show("xlow = -3; xupp = -1; niter = minniter")
                    xlow = -3; xupp = -1; niter = minniter
    call disp%show("call setBracketMin(getSqNeg, niter, xmin, xlow, xupp, fmin, flow, fupp)")
                    call setBracketMin(getSqNeg, niter, xmin, xlow, xupp, fmin, flow, fupp)
    call disp%show("if (niter > minniter) error stop 'Bracketing failed.'")
                    if (niter > minniter) error stop 'Bracketing failed.'
    call disp%show("[xlow, xmin, xupp, flow, fmin, fupp]")
    call disp%show( [xlow, xmin, xupp, flow, fmin, fupp] )
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("call setAsserted(isBracketMin(xmin, xlow, xupp, flow, fmin, fupp))")
                    call setAsserted(isBracketMin(xmin, xlow, xupp, fmin, flow, fupp))
    call disp%skip()

    call disp%skip()
    call disp%show("getSqNeg(x) = (x - 1)**2")
    call disp%show("xlow = 3; xupp = 9; niter = minniter")
                    xlow = 3; xupp = 9; niter = minniter
    call disp%show("call setBracketMin(getSqNeg, niter, xmin, xlow, xupp, fmin, flow, fupp)")
                    call setBracketMin(getSqNeg, niter, xmin, xlow, xupp, fmin, flow, fupp)
    call disp%show("if (niter > minniter) error stop 'Bracketing failed.'")
                    if (niter > minniter) error stop 'Bracketing failed.'
    call disp%show("[xlow, xmin, xupp, flow, fmin, fupp]")
    call disp%show( [xlow, xmin, xupp, flow, fmin, fupp] )
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("call setAsserted(isBracketMin(xmin, xlow, xupp, fmin, flow, fupp))")
                    call setAsserted(isBracketMin(xmin, xlow, xupp, fmin, flow, fupp))
    call disp%skip()

    call disp%skip()
    call disp%show("getSqNeg(x) = (x - 1)**2")
    call disp%show("xlow = -5; xupp = 4; niter = minniter")
                    xlow = -5; xupp = 4; niter = minniter
    call disp%show("call setBracketMin(getSqNeg, niter, xmin, xlow, xupp, fmin, flow, fupp)")
                    call setBracketMin(getSqNeg, niter, xmin, xlow, xupp, fmin, flow, fupp)
    call disp%show("if (niter > minniter) error stop 'Bracketing failed.'")
                    if (niter > minniter) error stop 'Bracketing failed.'
    call disp%show("[xlow, xmin, xupp, flow, fmin, fupp]")
    call disp%show( [xlow, xmin, xupp, flow, fmin, fupp] )
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("call setAsserted(isBracketMin(xmin, xlow, xupp, fmin, flow, fupp))")
                    call setAsserted(isBracketMin(xmin, xlow, xupp, fmin, flow, fupp))
    call disp%skip()

contains

    function getSqNeg(x) result(func)
        real(RKC), intent(in) :: x
        real(RKC) :: func
        func = (x - 1)**2
    end function

end program example