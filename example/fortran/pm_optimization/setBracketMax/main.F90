program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_optimization, only: setBracketMax, isBracketMax
    use pm_err, only: setAsserted
    use pm_val2str, only: getStr

    implicit none

    integer(IK) :: niter, minniter = 1000
    real(RKG) :: xlow, xmax, xupp, flow, fmax, fupp
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getSqNeg(x) = (x - 1)**2")
    call disp%show("xlow = -3; xupp = -1; niter = minniter")
                    xlow = -3; xupp = -1; niter = minniter
    call disp%show("call setBracketMax(getSqNeg, niter, xmax, xlow, xupp, fmax)")
                    call setBracketMax(getSqNeg, niter, xmax, xlow, xupp, fmax)
    call disp%show("if (niter > minniter) error stop 'Bracketing failed.'")
                    if (niter > minniter) error stop 'Bracketing failed.'
    call disp%show("[xlow, xmax, xupp, getSqNeg(xlow), fmax, getSqNeg(xupp)]")
    call disp%show( [xlow, xmax, xupp, getSqNeg(xlow), fmax, getSqNeg(xupp)] )
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("call setAsserted(isBracketMax(xmax, xlow, xupp, fmax, getSqNeg(xlow), getSqNeg(xupp)))")
                    call setAsserted(isBracketMax(xmax, xlow, xupp, fmax, getSqNeg(xlow), getSqNeg(xupp)))
    call disp%skip()

    call disp%skip()
    call disp%show("getSqNeg(x) = (x - 1)**2")
    call disp%show("xlow = -3; xupp = -1; niter = minniter")
                    xlow = -3; xupp = -1; niter = minniter
    call disp%show("call setBracketMax(getSqNeg, niter, xmax, xlow, xupp, fmax, flow, fupp)")
                    call setBracketMax(getSqNeg, niter, xmax, xlow, xupp, fmax, flow, fupp)
    call disp%show("if (niter > minniter) error stop 'Bracketing failed.'")
                    if (niter > minniter) error stop 'Bracketing failed.'
    call disp%show("[xlow, xmax, xupp, flow, fmax, fupp]")
    call disp%show( [xlow, xmax, xupp, flow, fmax, fupp] )
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("call setAsserted(isBracketMax(xmax, xlow, xupp, flow, fmax, fupp))")
                    call setAsserted(isBracketMax(xmax, xlow, xupp, fmax, flow, fupp))
    call disp%skip()

    call disp%skip()
    call disp%show("getSqNeg(x) = (x - 1)**2")
    call disp%show("xlow = 3; xupp = 9; niter = minniter")
                    xlow = 3; xupp = 9; niter = minniter
    call disp%show("call setBracketMax(getSqNeg, niter, xmax, xlow, xupp, fmax, flow, fupp)")
                    call setBracketMax(getSqNeg, niter, xmax, xlow, xupp, fmax, flow, fupp)
    call disp%show("if (niter > minniter) error stop 'Bracketing failed.'")
                    if (niter > minniter) error stop 'Bracketing failed.'
    call disp%show("[xlow, xmax, xupp, flow, fmax, fupp]")
    call disp%show( [xlow, xmax, xupp, flow, fmax, fupp] )
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("call setAsserted(isBracketMax(xmax, xlow, xupp, fmax, flow, fupp))")
                    call setAsserted(isBracketMax(xmax, xlow, xupp, fmax, flow, fupp))
    call disp%skip()

    call disp%skip()
    call disp%show("getSqNeg(x) = (x - 1)**2")
    call disp%show("xlow = -5; xupp = 4; niter = minniter")
                    xlow = -5; xupp = 4; niter = minniter
    call disp%show("call setBracketMax(getSqNeg, niter, xmax, xlow, xupp, fmax, flow, fupp)")
                    call setBracketMax(getSqNeg, niter, xmax, xlow, xupp, fmax, flow, fupp)
    call disp%show("if (niter > minniter) error stop 'Bracketing failed.'")
                    if (niter > minniter) error stop 'Bracketing failed.'
    call disp%show("[xlow, xmax, xupp, flow, fmax, fupp]")
    call disp%show( [xlow, xmax, xupp, flow, fmax, fupp] )
    call disp%show("niter")
    call disp%show( niter )
    call disp%show("call setAsserted(isBracketMax(xmax, xlow, xupp, fmax, flow, fupp))")
                    call setAsserted(isBracketMax(xmax, xlow, xupp, fmax, flow, fupp))
    call disp%skip()

contains

    function getSqNeg(x) result(func)
        real(RKG), intent(in) :: x
        real(RKG) :: func
        func = -(x - 1)**2
    end function

end program example