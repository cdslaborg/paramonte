program example

    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_kind, only: IK, CK32, CK64, CK128, RK32, RK64, RK128
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
    call disp%show("exp( getLogSumExp( Array = log([1._RK32, 2._RK32, 3._RK32, 4._RK32]), maxArray = log(4._RK32) ) )")
    call disp%show( exp( getLogSumExp( Array = log([1._RK32, 2._RK32, 3._RK32, 4._RK32]), maxArray = log(4._RK32) ) ) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%")
    call disp%show("64-bit real")
    call disp%show("!%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( Array = log([1._RK64, 2._RK64, 3._RK64, 4._RK64]), maxArray = log(4._RK64) ) )")
    call disp%show( exp( getLogSumExp( Array = log([1._RK64, 2._RK64, 3._RK64, 4._RK64]), maxArray = log(4._RK64) ) ) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%")
    call disp%show("128-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( Array = log([1._RK128, 2._RK128, 3._RK128, 4._RK128]), maxArray = log(4._RK128) ) )")
    call disp%show( exp( getLogSumExp( Array = log([1._RK128, 2._RK128, 3._RK128, 4._RK128]), maxArray = log(4._RK128) ) ) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("32-bit complex")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( Array = cmplx(log([(1._CK32,-1._CK32), (2._CK32,-2._CK32), (3._CK32,-3._CK32)]), kind = CK32), maxArray = log((3._CK32,-3._CK32)) ) )")
    call disp%show( exp( getLogSumExp( Array = cmplx(log([(1._CK32,-1._CK32), (2._CK32,-2._CK32), (3._CK32,-3._CK32)]), kind = CK32), maxArray = log((3._CK32,-3._CK32)) ) ) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("64-bit complex")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( Array = cmplx(log([(1._CK64,-1._CK64), (2._CK64,-2._CK64), (3._CK64,-3._CK64)]), kind = CK64), maxArray = log((3._CK64,-3._CK64)) ) )")
    call disp%show( exp( getLogSumExp( Array = cmplx(log([(1._CK64,-1._CK64), (2._CK64,-2._CK64), (3._CK64,-3._CK64)]), kind = CK64), maxArray = log((3._CK64,-3._CK64)) ) ) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("exp( getLogSumExp( Array = cmplx(log([(1._CK128,-1._CK128), (2._CK128,-2._CK128), (3._CK128,-3._CK128)]), kind = CK128), maxArray = log((3._CK128,-3._CK128)) ) )")
    call disp%show( exp( getLogSumExp( Array = cmplx(log([(1._CK128,-1._CK128), (2._CK128,-2._CK128), (3._CK128,-3._CK128)]), kind = CK128), maxArray = log((3._CK128,-3._CK128)) ) ) )
    call disp%skip()

end program example