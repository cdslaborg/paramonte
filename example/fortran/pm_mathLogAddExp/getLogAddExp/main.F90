program example

    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_kind, only: IK, CK32, CK64, CK128, RK32, RK64, RK128
    use pm_mathLogAddExp, only: getLogAddExp
    use pm_distUnif, only: getUnifRand
    use pm_mathMinMax, only: getMinMax

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!32-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log(4._RK32), larger = log(10._RK32)) )")
    call disp%show( exp( getLogAddExp(smaller = log(4._RK32), larger = log(10._RK32)) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!64-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log(4._RK64), larger = log(10._RK64)) )")
    call disp%show( exp( getLogAddExp(smaller = log(4._RK64), larger = log(10._RK64)) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%")
    call disp%show("!128-bit real")
    call disp%show("!%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log(4._RK128), larger = log(10._RK128)) )")
    call disp%show( exp( getLogAddExp(smaller = log(4._RK128), larger = log(10._RK128)) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!32-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log((4._CK32, -4._CK32)), larger = log((10._CK32, -10._CK32))) )")
    call disp%show( exp( getLogAddExp(smaller = log((4._CK32, -4._CK32)), larger = log((10._CK32, -10._CK32))) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log((4._CK64, -4._CK64)), larger = log((10._CK64, -10._CK64))) )")
    call disp%show( exp( getLogAddExp(smaller = log((4._CK64, -4._CK64)), larger = log((10._CK64, -10._CK64))) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log((4._CK128, -4._CK128)), larger = log((10._CK128, -10._CK128))) )")
    call disp%show( exp( getLogAddExp(smaller = log((4._CK128, -4._CK128)), larger = log((10._CK128, -10._CK128))) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Input arguments as a vector of [minimum, maximum] by calling getMinMax()")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp(getLogAddExp(getMinMax(getUnifRand(-1., 1., 2_IK))))")
    call disp%show( exp(getLogAddExp(getMinMax(getUnifRand(-1., 1., 2_IK)))) )

end program example