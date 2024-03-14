program example

    use pm_kind, only: SK, IK, LK
    use pm_distPower, only: getPowerLogQuan
    use pm_io, only: display_type

    implicit none

    real                    :: LogX(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1) = getPowerLogQuan(logCDF = -3., Alpha = +2., logMaxX = -2.) ! Power distribution.")
                    LogX(1) = getPowerLogQuan(logCDF = -3., Alpha = +2., logMaxX = -2.) ! Power distribution.
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1:3) = getPowerLogQuan(logCDF = -[3., 4., 5.], Alpha = +[+2., +3., +4.], logMaxX = -2.) ! Power distribution.")
                    LogX(1:3) = getPowerLogQuan(logCDF = -[3., 4., 5.], Alpha = +[+2., +3., +4.], logMaxX = -2.) ! Power distribution.
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Truncated Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1) = getPowerLogQuan(logCDF = -3., Alpha = +1., logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.")
                    LogX(1) = getPowerLogQuan(logCDF = -3., Alpha = +1., logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1:3) = getPowerLogQuan(logCDF = -[3., 1., 2.], Alpha = +[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.")
                    LogX(1:3) = getPowerLogQuan(logCDF = -[3., 1., 2.], Alpha = +[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example LogX array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        use pm_distPower, only: getPowerLogCDF
        real :: Alpha(2), logMinX, logMaxX, LogX(4), LogCDF(2000)
        integer(IK) :: fileUnit, i
        call setLinSpace(LogCDF, x1 = log(.001), x2 = log(.999))
        Alpha = [+0.5, +2.0]
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "getPowerLogQuan.RK.txt")
        do i = 1, size(LogCDF, 1, IK)
            LogX(1:2) = getPowerLogQuan(LogCDF(i), Alpha, logMaxX)
            LogX(3:4) = getPowerLogQuan(LogCDF(i), Alpha, logMinX, logMaxX)
            write(fileUnit, "(5(g0,:,', '))") exp(LogCDF(i)), exp(LogX(1:4))
        end do
        close(fileUnit)
    end block

end program example