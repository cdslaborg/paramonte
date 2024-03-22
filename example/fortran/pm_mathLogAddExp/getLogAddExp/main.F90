program example

    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_kind, only: IK, CKS, CKD, CKH, RKS, RKD, RKH
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

    call disp%show("exp( getLogAddExp(smaller = log(4._RKS), larger = log(10._RKS)) )")
    call disp%show( exp( getLogAddExp(smaller = log(4._RKS), larger = log(10._RKS)) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!64-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log(4._RKD), larger = log(10._RKD)) )")
    call disp%show( exp( getLogAddExp(smaller = log(4._RKD), larger = log(10._RKD)) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%")
    call disp%show("!128-bit real")
    call disp%show("!%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log(4._RKH), larger = log(10._RKH)) )")
    call disp%show( exp( getLogAddExp(smaller = log(4._RKH), larger = log(10._RKH)) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!32-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log((4._CKS, -4._CKS)), larger = log((10._CKS, -10._CKS))) )")
    call disp%show( exp( getLogAddExp(smaller = log((4._CKS, -4._CKS)), larger = log((10._CKS, -10._CKS))) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log((4._CKD, -4._CKD)), larger = log((10._CKD, -10._CKD))) )")
    call disp%show( exp( getLogAddExp(smaller = log((4._CKD, -4._CKD)), larger = log((10._CKD, -10._CKD))) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp( getLogAddExp(smaller = log((4._CKH, -4._CKH)), larger = log((10._CKH, -10._CKH))) )")
    call disp%show( exp( getLogAddExp(smaller = log((4._CKH, -4._CKH)), larger = log((10._CKH, -10._CKH))) ) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Input arguments as a vector of [minimum, maximum] by calling getMinMax()")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("exp(getLogAddExp(getMinMax(getUnifRand(-1., 1., 2_IK))))")
    call disp%show( exp(getLogAddExp(getMinMax(getUnifRand(-1., 1., 2_IK)))) )

end program example