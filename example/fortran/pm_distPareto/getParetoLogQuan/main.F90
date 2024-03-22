program example

    use pm_kind, only: SK, IK, LK
    use pm_distPareto, only: getParetoLogQuan
    use pm_io, only: display_type

    implicit none

    real                    :: logx(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = getParetoLogQuan(logCDF = -3., alpha = -2., logMinX = -2.) ! Pareto distribution.")
                    logx(1) = getParetoLogQuan(logCDF = -3., alpha = -2., logMinX = -2.) ! Pareto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1:3) = getParetoLogQuan(logCDF = -[3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.")
                    logx(1:3) = getParetoLogQuan(logCDF = -[3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Truncated Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = getParetoLogQuan(logCDF = -3., alpha = -1., logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.")
                    logx(1) = getParetoLogQuan(logCDF = -3., alpha = -1., logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1:3) = getParetoLogQuan(logCDF = -[3., 1., 2.], alpha = -[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.")
                    logx(1:3) = getParetoLogQuan(logCDF = -[3., 1., 2.], alpha = -[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        use pm_distPareto, only: getParetoLogCDF
        real :: alpha(2), logMinX, logMaxX, logx(4), logCDF(2000)
        integer(IK) :: fileUnit, i
        call setLinSpace(logCDF, x1 = log(.001), x2 = log(.999))
        alpha = -[+0.5, +2.0]
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "getParetoLogQuan.RK.txt")
        do i = 1, size(logCDF, 1, IK)
            logx(1:2) = getParetoLogQuan(logCDF(i), alpha, logMinX)
            logx(3:4) = getParetoLogQuan(logCDF(i), alpha, logMinX, logMaxX)
            write(fileUnit, "(5(g0,:,', '))") exp(logCDF(i)), exp(logx(1:4))
        end do
        close(fileUnit)
    end block

end program example