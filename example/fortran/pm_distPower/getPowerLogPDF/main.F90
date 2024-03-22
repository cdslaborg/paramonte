program example

    use pm_kind, only: SK, IK, LK
    use pm_distPower, only: getPowerLogPDF
    use pm_io, only: display_type

    implicit none

    real                    :: logPDF(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getPowerLogPDF(logx = 3., alpha = +2., logMaxX = 5.) ! Power distribution.")
                    logPDF(1) = getPowerLogPDF(logx = 3., alpha = +2., logMaxX = 5.) ! Power distribution.
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1:3) = getPowerLogPDF(logx = [3., 4., 5.], alpha = [+2., +3., +4.], logMaxX = 5.) ! Power distribution.")
                    logPDF(1:3) = getPowerLogPDF(logx = [3., 4., 5.], alpha = [+2., +3., +4.], logMaxX = 5.) ! Power distribution.
    call disp%show("logPDF(1:3)")
    call disp%show( logPDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Truncated Power distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getPowerLogPDF(logx = 3., alpha = 1., logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.")
                    logPDF(1) = getPowerLogPDF(logx = 3., alpha = 1., logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1:3) = getPowerLogPDF(logx = [3., 0., -2.], alpha = [+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.")
                    logPDF(1:3) = getPowerLogPDF(logx = [3., 0., -2.], alpha = [+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Power distribution.
    call disp%show("logPDF(1:3)")
    call disp%show( logPDF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example logPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: alpha(2), logMinX, logMaxX, logPDF(4), logx(2000)
        integer(IK) :: fileUnit, i
        call setLinSpace(logx, x1 = log(0.1), x2 = log(10.))
        alpha = [+0.5, +2.0]
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "getPowerLogPDF.RK.txt")
        do i = 1, size(logx, 1, IK)
            logPDF = -huge(0.)
            if (logx(i) <= logMaxX) then
                logPDF(1:2) = getPowerLogPDF(logx(i), alpha, logMaxX)
                if (logMinX <= logx(i)) logPDF(3:4) = getPowerLogPDF(logx(i), alpha, logMinX, logMaxX)
            end if
            write(fileUnit, "(*(g0,:,', '))") exp(logx(i)), exp(logPDF)
        end do
        close(fileUnit)
    end block

end program example