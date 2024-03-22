program example

    use pm_kind, only: SK, IK, LK
    use pm_distPower, only: getPowerLogCDFNF
    use pm_distPower, only: setPowerLogRand
    use pm_distNegExp, only: getNegExpRand
    use pm_io, only: display_type

    implicit none

    real                    :: logx(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = +2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.")
                    call setPowerLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = +2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = [+2., +3., +4.], logCDFNF = getPowerLogCDFNF(alpha = [+2., +3., +4.], logMaxX = -2.)) ! Power distribution.")
                    call setPowerLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = [+2., +3., +4.], logCDFNF = getPowerLogCDFNF(alpha = [+2., +3., +4.], logMaxX = -2.)) ! Power distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Truncated Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = +2., logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.")
                    call setPowerLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = +2., logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = +[+2., +3., +4.], logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.")
                    call setPowerLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = +[+2., +3., +4.], logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: alpha(3), logMinX(3), logMaxX(3), logx(3)
        integer(IK) :: fileUnit, i
        alpha = [.1, 2., .5]
        logMinX = log([3., tiny(0.), tiny(0.)])
        logMaxX = log([10., 5., 8.])
        open(newunit = fileUnit, file = "setPowerLogRand.RK.txt")
        do i = 1, 2000_IK
            call setPowerLogRand(logx(1), getNegExpRand(1.), alpha(1), logMinX(1), getPowerLogCDFNF(alpha(1), logMinX(1), logMaxX(1)))
            call setPowerLogRand(logx(2:3), getNegExpRand([1., 1.]), alpha(2:3), getPowerLogCDFNF(alpha(2:3), logMaxX(2:3)))
            write(fileUnit, "(*(g0,:,', '))") exp(logx)
        end do
        close(fileUnit)
    end block

end program example