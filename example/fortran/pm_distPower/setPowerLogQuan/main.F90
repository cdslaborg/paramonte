program example

    use pm_kind, only: SK, IK, LK
    use pm_distPower, only: getPowerLogCDFNF
    use pm_distPower, only: setPowerLogQuan
    use pm_io, only: display_type

    implicit none

    real                    :: logx(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogQuan(logx(1), logCDF = -3., alpha = +2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.")
                    call setPowerLogQuan(logx(1), logCDF = -3., alpha = +2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogQuan(logx(1:3), logCDF = -[3., 4., 5.], alpha = +[+2., +3., +4.], logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.")
                    call setPowerLogQuan(logx(1:3), logCDF = -[3., 4., 5.], alpha = +[+2., +3., +4.], logCDFNF = getPowerLogCDFNF(alpha = +2., logMaxX = -2.)) ! Power distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Truncated Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogQuan(logx(1), logCDF = -3., alpha = +2., logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.")
                    call setPowerLogQuan(logx(1), logCDF = -3., alpha = +2., logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowerLogQuan(logx(1:3), logCDF = -[3., 4., 5.], alpha = +[+2., +3., +4.], logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.")
                    call setPowerLogQuan(logx(1:3), logCDF = -[3., 4., 5.], alpha = +[+2., +3., +4.], logMinX = -2., logCDFNF = getPowerLogCDFNF(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Power distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        integer(IK) :: fileUnit, i
        real :: alpha(2), logMinX, logMaxX, logx(4), logCDF(2000)
        call setLinSpace(logCDF, x1 = log(0.001), x2 = log(.999))
        alpha = [+0.5, +2.0]
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "setPowerLogQuan.RK.txt")
        do i = 1, size(logCDF, 1, IK)
            call setPowerLogQuan(logx(1:2), logCDF(i), alpha, getPowerLogCDFNF(alpha, logMaxX))
            call setPowerLogQuan(logx(3:4), logCDF(i), alpha, logMinX, getPowerLogCDFNF(alpha, logMinX, logMaxX))
            write(fileUnit, "(*(g0,:,', '))") exp(logCDF(i)), exp(logx)
        end do
        close(fileUnit)
    end block

end program example