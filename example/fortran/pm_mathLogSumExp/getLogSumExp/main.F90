program example

    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_kind, only: IK, CKS, CKD, CKH, RKS, RKD, RKH
    use pm_mathLogSumExp, only: getLogSumExp 

    implicit none

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%")
    call disp%show("32-bit real")
    call disp%show("!%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( array = log([1._RKS, 2._RKS, 3._RKS, 4._RKS]), maxArray = log(4._RKS) ) )")
    call disp%show( exp( getLogSumExp( array = log([1._RKS, 2._RKS, 3._RKS, 4._RKS]), maxArray = log(4._RKS) ) ) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%")
    call disp%show("64-bit real")
    call disp%show("!%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( array = log([1._RKD, 2._RKD, 3._RKD, 4._RKD]), maxArray = log(4._RKD) ) )")
    call disp%show( exp( getLogSumExp( array = log([1._RKD, 2._RKD, 3._RKD, 4._RKD]), maxArray = log(4._RKD) ) ) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%")
    call disp%show("128-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( array = log([1._RKH, 2._RKH, 3._RKH, 4._RKH]), maxArray = log(4._RKH) ) )")
    call disp%show( exp( getLogSumExp( array = log([1._RKH, 2._RKH, 3._RKH, 4._RKH]), maxArray = log(4._RKH) ) ) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("32-bit complex")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( array = cmplx(log([(1._CKS,-1._CKS), (2._CKS,-2._CKS), (3._CKS,-3._CKS)]), kind = CKS), maxArray = log((3._CKS,-3._CKS)) ) )")
    call disp%show( exp( getLogSumExp( array = cmplx(log([(1._CKS,-1._CKS), (2._CKS,-2._CKS), (3._CKS,-3._CKS)]), kind = CKS), maxArray = log((3._CKS,-3._CKS)) ) ) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("64-bit complex")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( array = cmplx(log([(1._CKD,-1._CKD), (2._CKD,-2._CKD), (3._CKD,-3._CKD)]), kind = CKD), maxArray = log((3._CKD,-3._CKD)) ) )")
    call disp%show( exp( getLogSumExp( array = cmplx(log([(1._CKD,-1._CKD), (2._CKD,-2._CKD), (3._CKD,-3._CKD)]), kind = CKD), maxArray = log((3._CKD,-3._CKD)) ) ) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( array = cmplx(log([(1._CKH,-1._CKH), (2._CKH,-2._CKH), (3._CKH,-3._CKH)]), kind = CKH), maxArray = log((3._CKH,-3._CKH)) ) )")
    call disp%show( exp( getLogSumExp( array = cmplx(log([(1._CKH,-1._CKH), (2._CKH,-2._CKH), (3._CKH,-3._CKH)]), kind = CKH), maxArray = log((3._CKH,-3._CKH)) ) ) )
    call disp%skip()

end program example