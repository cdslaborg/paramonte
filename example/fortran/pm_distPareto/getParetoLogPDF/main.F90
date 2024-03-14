program example

    use pm_kind, only: SK, IK, LK
    use pm_distPareto, only: getParetoLogPDF
    use pm_io, only: display_type

    implicit none

    real                    :: LogPDF(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1) = getParetoLogPDF(logx = 3., alpha = -2., logMinX = -2.) ! Pareto distribution.")
                    LogPDF(1) = getParetoLogPDF(logx = 3., alpha = -2., logMinX = -2.) ! Pareto distribution.
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1:3) = getParetoLogPDF(logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.")
                    LogPDF(1:3) = getParetoLogPDF(logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.
    call disp%show("LogPDF(1:3)")
    call disp%show( LogPDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Truncated Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1) = getParetoLogPDF(logx = 3., alpha = -1., logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.")
                    LogPDF(1) = getParetoLogPDF(logx = 3., alpha = -1., logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1:3) = getParetoLogPDF(logx = [3., 0., -2.], alpha = -[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.")
                    LogPDF(1:3) = getParetoLogPDF(logx = [3., 0., -2.], alpha = -[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.
    call disp%show("LogPDF(1:3)")
    call disp%show( LogPDF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example LogPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: Alpha(2), logMinX, logMaxX, LogPDF(4), LogX(2000)
        integer(IK) :: fileUnit, i
        call setLinSpace(LogX, x1 = log(0.1), x2 = log(10.))
        Alpha = -[+0.5, +2.0]
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "getParetoLogPDF.RK.txt")
        do i = 1, size(LogX, 1, IK)
            LogPDF = -huge(0.)
            if (logMinX <= LogX(i)) then
                LogPDF(1:2) = getParetoLogPDF(LogX(i), alpha, logMinX)
                if (LogX(i) <= logMaxX) LogPDF(3:4) = getParetoLogPDF(LogX(i), alpha, logMinX, logMaxX)
            end if
            write(fileUnit, "(5(g0,:,', '))") exp(LogX(i)), exp(LogPDF)
        end do
        close(fileUnit)
    end block

end program example