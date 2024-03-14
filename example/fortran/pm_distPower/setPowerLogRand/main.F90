program example

    use pm_kind, only: SK, IK, LK
    use pm_distPower, only: getPowerLogCDFNF
    use pm_distPower, only: setPowerLogRand
    use pm_distNegExp, only: getNegExpRand
    use pm_io, only: display_type

    implicit none

    real                    :: LogX(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogRand(LogX(1), negExpRand = getNegExpRand(1.), alpha = +2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.")
                    call setPowerLogRand(LogX(1), negExpRand = getNegExpRand(1.), alpha = +2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogRand(LogX(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = [+2., +3., +4.], logCDFNF = getPowerLogCDFNF(alpha = [+2., +3., +4.], logMaxX = -2.)) ! Power distribution.")
                    call setPowerLogRand(LogX(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = [+2., +3., +4.], logCDFNF = getPowerLogCDFNF(alpha = [+2., +3., +4.], logMaxX = -2.)) ! Power distribution.
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Truncated Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogRand(LogX(1), negExpRand = getNegExpRand(1.), alpha = +2., logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.")
                    call setPowerLogRand(LogX(1), negExpRand = getNegExpRand(1.), alpha = +2., logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogRand(LogX(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = +[+2., +3., +4.], logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.")
                    call setPowerLogRand(LogX(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = +[+2., +3., +4.], logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example LogX array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: Alpha(3), LogMinX(3), LogMaxX(3), LogX(3)
        integer(IK) :: fileUnit, i
        Alpha = [.1, 2., .5]
        LogMinX = log([3., tiny(0.), tiny(0.)])
        LogMaxX = log([10., 5., 8.])
        open(newunit = fileUnit, file = "setPowerLogRand.RK.txt")
        do i = 1, 2000_IK
            call setPowerLogRand(LogX(1), getNegExpRand(1.), Alpha(1), LogMinX(1), getPowerLogCDFNF(Alpha(1), LogMinX(1), LogMaxX(1)))
            call setPowerLogRand(LogX(2:3), getNegExpRand([1., 1.]), Alpha(2:3), getPowerLogCDFNF(Alpha(2:3), LogMaxX(2:3)))
            write(fileUnit, "(*(g0,:,', '))") exp(LogX)
        end do
        close(fileUnit)
    end block

end program example