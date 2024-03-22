program example

    use pm_kind, only: SK, IK, LK
    use pm_distPareto, only: getParetoLogCDFNF
    use pm_distPareto, only: setParetoLogRand
    use pm_distNegExp, only: getNegExpRand
    use pm_io, only: display_type

    implicit none

    real                    :: logx(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setParetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = -2., logMinX = -2.) ! Pareto distribution.")
                    call setParetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = -2., logMinX = -2.) ! Pareto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setParetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.")
                    call setParetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Truncated Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setParetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = -2., logMinX = -2., logCDFNF = getParetoLogCDFNF(alpha = -2., logMinX = -2., logMaxX = 5.)) ! Truncated Pareto distribution.")
                    call setParetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = -2., logMinX = -2., logCDFNF = getParetoLogCDFNF(alpha = -2., logMinX = -2., logMaxX = 5.)) ! Truncated Pareto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setParetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = -[+2., +3., +4.], logMinX = -2., logCDFNF = getParetoLogCDFNF(alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Pareto distribution.")
                    call setParetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = -[+2., +3., +4.], logMinX = -2., logCDFNF = getParetoLogCDFNF(alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Pareto distribution.
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
        alpha = -[.5, 1., 10.]
        logMinX = log([3., 1., 2.])
        logMaxX = log([5., 4., huge(0.)])
        open(newunit = fileUnit, file = "setParetoLogRand.RK.txt")
        do i = 1, 2000_IK
            call setParetoLogRand(logx(1:2), getNegExpRand([1., 1.]), alpha(1:2), logMinX(1:2), getParetoLogCDFNF(alpha(1:2), logMinX(1:2), logMaxX(1:2)))
            call setParetoLogRand(logx(3), getNegExpRand(1.), alpha(3), logMinX(3), 0.)
            write(fileUnit, "(*(g0,:,', '))") exp(logx)
        end do
        close(fileUnit)
    end block

end program example