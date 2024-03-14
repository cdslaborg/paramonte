program example

    use pm_kind, only: SK, IK, LK
    use pm_distPareto, only: getParetoLogQuan
    use pm_io, only: display_type

    implicit none

    real                    :: LogX(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1) = getParetoLogQuan(logCDF = -3., Alpha = -2., logMinX = -2.) ! Pareto distribution.")
                    LogX(1) = getParetoLogQuan(logCDF = -3., Alpha = -2., logMinX = -2.) ! Pareto distribution.
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1:3) = getParetoLogQuan(logCDF = -[3., 4., 5.], Alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.")
                    LogX(1:3) = getParetoLogQuan(logCDF = -[3., 4., 5.], Alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Truncated Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1) = getParetoLogQuan(logCDF = -3., Alpha = -1., logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.")
                    LogX(1) = getParetoLogQuan(logCDF = -3., Alpha = -1., logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1:3) = getParetoLogQuan(logCDF = -[3., 1., 2.], Alpha = -[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.")
                    LogX(1:3) = getParetoLogQuan(logCDF = -[3., 1., 2.], Alpha = -[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example LogX array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        use pm_distPareto, only: getParetoLogCDF
        real :: Alpha(2), logMinX, logMaxX, LogX(4), LogCDF(2000)
        integer(IK) :: fileUnit, i
        call setLinSpace(LogCDF, x1 = log(.001), x2 = log(.999))
        Alpha = -[+0.5, +2.0]
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "getParetoLogQuan.RK.txt")
        do i = 1, size(LogCDF, 1, IK)
            LogX(1:2) = getParetoLogQuan(LogCDF(i), Alpha, logMinX)
            LogX(3:4) = getParetoLogQuan(LogCDF(i), Alpha, logMinX, logMaxX)
            write(fileUnit, "(5(g0,:,', '))") exp(LogCDF(i)), exp(LogX(1:4))
        end do
        close(fileUnit)
    end block

end program example