program example

    use pm_kind, only: SK, IK, LK
    use pm_distPower, only: getPowerLogCDFNF
    use pm_distPower, only: setPowerLogQuan
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
    call disp%show("call setPowerLogQuan(LogX(1), LogCDF = -3., alpha = +2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.")
                    call setPowerLogQuan(LogX(1), LogCDF = -3., alpha = +2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogQuan(LogX(1:3), LogCDF = -[3., 4., 5.], alpha = +[+2., +3., +4.], logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.")
                    call setPowerLogQuan(LogX(1:3), LogCDF = -[3., 4., 5.], alpha = +[+2., +3., +4.], logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Truncated Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogQuan(LogX(1), LogCDF = -3., alpha = +2., logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.")
                    call setPowerLogQuan(LogX(1), LogCDF = -3., alpha = +2., logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogQuan(LogX(1:3), LogCDF = -[3., 4., 5.], alpha = +[+2., +3., +4.], logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.")
                    call setPowerLogQuan(LogX(1:3), LogCDF = -[3., 4., 5.], alpha = +[+2., +3., +4.], logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example LogX array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        integer(IK) :: fileUnit, i
        real :: Alpha(2), LogMinX, LogMaxX, LogX(4), LogCDF(2000)
        call setLinSpace(LogCDF, x1 = log(0.001), x2 = log(.999))
        Alpha = [+0.5, +2.0]
        LogMinX = log(3.)
        LogMaxX = log(8.)
        open(newunit = fileUnit, file = "setPowerLogQuan.RK.txt")
        do i = 1, size(LogCDF, 1, IK)
            call setPowerLogQuan(LogX(1:2), LogCDF(i), Alpha, getPowerLogCDFNF(Alpha, LogMaxX))
            call setPowerLogQuan(LogX(3:4), LogCDF(i), Alpha, LogMinX, getPowerLogCDFNF(Alpha, LogMinX, LogMaxX))
            write(fileUnit, "(5(g0,:,', '))") exp(LogCDF(i)), exp(LogX)
        end do
        close(fileUnit)
    end block

end program example