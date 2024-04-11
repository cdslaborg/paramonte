program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RKH ! all processor kinds are supported.
    use pm_mathRoot, only: setRoot, false, bisection, secant, brent, ridders, toms748, newton, halley, schroder
    use pm_io, only: display_type

    implicit none

    real(RKC)   :: root
    real(RKC)   :: lb, ub
    real(RKC)   :: lf, uf
    logical(LK) :: failed
    integer(IK) :: neval
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("lb = 2._RKC; lf = getSin(lb)")
                    lb = 2._RKC; lf = getSin(lb)
    call disp%show("ub = 4._RKC; uf = getSin(ub)")
                    ub = 4._RKC; uf = getSin(ub)
    call disp%show("[lf, uf]")
    call disp%show( [lf, uf] )
    call disp%show("call setRoot(false, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(false, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getSin(root)]")
    call disp%show( [root, getSin(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(bisection, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(bisection, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getSin(root)]")
    call disp%show( [root, getSin(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(secant, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(secant, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getSin(root)]")
    call disp%show( [root, getSin(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(brent, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(brent, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getSin(root)]")
    call disp%show( [root, getSin(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(ridders, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(ridders, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getSin(root)]")
    call disp%show( [root, getSin(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(toms748, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(toms748, getSin, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getSin(root)]")
    call disp%show( [root, getSin(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(newton, getSinDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(newton, getSinDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getSin(root)]")
    call disp%show( [root, getSin(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(halley, getSinDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(halley, getSinDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getSin(root)]")
    call disp%show( [root, getSin(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(schroder, getSinDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(schroder, getSinDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getSin(root)]")
    call disp%show( [root, getSin(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()



    call disp%skip()
    call disp%show("lb = 1._RKC; lf = getCos(lb)")
                    lb = 1._RKC; lf = getCos(lb)
    call disp%show("ub = 3._RKC; uf = getCos(ub)")
                    ub = 3._RKC; uf = getCos(ub)
    call disp%show("[lf, uf]")
    call disp%show( [lf, uf] )
    call disp%show("call setRoot(false, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(false, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getCos(root)]")
    call disp%show( [root, getCos(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(bisection, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(bisection, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getCos(root)]")
    call disp%show( [root, getCos(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(secant, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(secant, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getCos(root)]")
    call disp%show( [root, getCos(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(brent, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(brent, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getCos(root)]")
    call disp%show( [root, getCos(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(ridders, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(ridders, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getCos(root)]")
    call disp%show( [root, getCos(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(toms748, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(toms748, getCos, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getCos(root)]")
    call disp%show( [root, getCos(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(newton, getCosDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(newton, getCosDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getCos(root)]")
    call disp%show( [root, getCos(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(halley, getCosDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(halley, getCosDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getCos(root)]")
    call disp%show( [root, getCos(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(schroder, getCosDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(schroder, getCosDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getCos(root)]")
    call disp%show( [root, getCos(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()



    call disp%skip()
    call disp%show("lb = -1._RKC; lf = getQuad(lb)")
                    lb = -1._RKC; lf = getQuad(lb)
    call disp%show("ub = +4._RKC; uf = getQuad(ub)")
                    ub = +4._RKC; uf = getQuad(ub)
    call disp%show("[lf, uf]")
    call disp%show( [lf, uf] )
    call disp%show("call setRoot(false, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(false, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getQuad(root)]")
    call disp%show( [root, getQuad(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(bisection, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(bisection, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getQuad(root)]")
    call disp%show( [root, getQuad(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(secant, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(secant, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getQuad(root)]")
    call disp%show( [root, getQuad(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(brent, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(brent, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getQuad(root)]")
    call disp%show( [root, getQuad(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(ridders, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(ridders, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getQuad(root)]")
    call disp%show( [root, getQuad(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(toms748, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(toms748, getQuad, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getQuad(root)]")
    call disp%show( [root, getQuad(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(newton, getQuadDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(newton, getQuadDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getQuad(root)]")
    call disp%show( [root, getQuad(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(halley, getQuadDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(halley, getQuadDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getQuad(root)]")
    call disp%show( [root, getQuad(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()
    call disp%show("call setRoot(schroder, getQuadDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)")
                    call setRoot(schroder, getQuadDiff, root, lb, ub, lf, uf, abstol = epsilon(0._RKC)**.7, neval = neval)
    call disp%show("if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'")
                    if (neval < 0_IK) error stop 'The root-finding algorithm failed to converge.'
    call disp%show("[root, getQuad(root)]")
    call disp%show( [root, getQuad(root)] )
    call disp%show("neval")
    call disp%show( neval )
    call disp%skip()

contains

    pure function getSin(x) result(func)
        real(RKC), intent(in) :: x
        real(RKC) :: func
        func = sin(x)
    end function

    pure function getSinDiff(x, order) result(func)
        integer(IK), intent(in) :: order
        real(RKC), intent(in) :: x
        real(RKC) :: func
        if (order == 0) func = +getSin(x)
        if (order == 1) func = +cos(x)
        if (order == 2) func = -sin(x)
    end function

    pure function getCos(x) result(func)
        real(RKC), intent(in) :: x
        real(RKC) :: func
        func = cos(x)
    end function

    pure function getCosDiff(x, order) result(func)
        integer(IK), intent(in) :: order
        real(RKC), intent(in) :: x
        real(RKC) :: func
        if (order == 0) func = +getCos(x)
        if (order == 1) func = -sin(x)
        if (order == 2) func = -cos(x)
    end function

    pure function getQuad(x) result(func)
        real(RKC), intent(in) :: x
        real(RKC) :: func
        func = x * (x - 1._RKC) * (x - 2._RKC)
    end function

    pure function getQuadDiff(x, order) result(func)
        integer(IK), intent(in) :: order
        real(RKC), intent(in) :: x
        real(RKC) :: func
        if (order == 0) func = getQuad(x)
        if (order == 1) func = 3._RKC * x**2 - 6._RKC * x + 2._RKC
        if (order == 2) func = 6._RKC * x - 6._RKC
    end function

end program example