program example

    use pm_kind, only: SK, IK, LK
    use pm_distPower, only: getPowerLogPDF
    use pm_io, only: display_type

    implicit none

    real                    :: LogPDF(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1) = getPowerLogPDF(logx = 3., alpha = +2., logMaxX = 5.) ! Power distribution.")
                    LogPDF(1) = getPowerLogPDF(logx = 3., alpha = +2., logMaxX = 5.) ! Power distribution.
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1:3) = getPowerLogPDF(logx = [3., 4., 5.], alpha = [+2., +3., +4.], logMaxX = 5.) ! Power distribution.")
                    LogPDF(1:3) = getPowerLogPDF(logx = [3., 4., 5.], alpha = [+2., +3., +4.], logMaxX = 5.) ! Power distribution.
    call disp%show("LogPDF(1:3)")
    call disp%show( LogPDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Truncated Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1) = getPowerLogPDF(logx = 3., alpha = 1., logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.")
                    LogPDF(1) = getPowerLogPDF(logx = 3., alpha = 1., logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1:3) = getPowerLogPDF(logx = [3., 0., -2.], alpha = [+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.")
                    LogPDF(1:3) = getPowerLogPDF(logx = [3., 0., -2.], alpha = [+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.
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
        alpha = [+0.5, +2.0]
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "getPowerLogPDF.RK.txt")
        do i = 1, size(LogX, 1, IK)
            LogPDF = -huge(0.)
            if (LogX(i) <= logMaxX) then
                LogPDF(1:2) = getPowerLogPDF(LogX(i), alpha, logMaxX)
                if (logMinX <= LogX(i)) LogPDF(3:4) = getPowerLogPDF(LogX(i), alpha, logMinX, logMaxX)
            end if
            write(fileUnit, "(5(g0,:,', '))") exp(LogX(i)), exp(LogPDF)
        end do
        close(fileUnit)
    end block

end program example