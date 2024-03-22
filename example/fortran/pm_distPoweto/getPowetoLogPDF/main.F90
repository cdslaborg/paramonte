program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLogSpace
    use pm_distPoweto, only: getPowetoLogPDF
    use pm_distPoweto, only: getPowetoLogPDFNF

    implicit none

    real :: logPDF(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getPowetoLogPDF(logx = 3., alpha = -2., logMinX = -2.) ! Poweto distribution.")
                    logPDF(1) = getPowetoLogPDF(logx = 3., alpha = -2., logMinX = -2.) ! Poweto distribution.
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1:3) = getPowetoLogPDF(logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2.) ! Poweto distribution.")
                    logPDF(1:3) = getPowetoLogPDF(logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2.) ! Poweto distribution.
    call disp%show("logPDF(1:3)")
    call disp%show( logPDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Truncated Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getPowetoLogPDF(logx = 3., alpha = -2., logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.")
                    logPDF(1) = getPowetoLogPDF(logx = 3., alpha = -2., logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1:3) = getPowetoLogPDF(logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.")
                    logPDF(1:3) = getPowetoLogPDF(logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.
    call disp%show("logPDF(1:3)")
    call disp%show( logPDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getPowetoLogPDF(logx = 3., alpha = +2., logMaxX = 5.) ! Poweto distribution.")
                    logPDF(1) = getPowetoLogPDF(logx = 3., alpha = +2., logMaxX = 5.) ! Poweto distribution.
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1:3) = getPowetoLogPDF(logx = [3., 4., 5.], alpha = [+2., +3., +4.], logMaxX = 5.) ! Poweto distribution.")
                    logPDF(1:3) = getPowetoLogPDF(logx = [3., 4., 5.], alpha = [+2., +3., +4.], logMaxX = 5.) ! Poweto distribution.
    call disp%show("logPDF(1:3)")
    call disp%show( logPDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Truncated Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getPowetoLogPDF(logx = 3., alpha = +2., logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.")
                    logPDF(1) = getPowetoLogPDF(logx = 3., alpha = +2., logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1:3) = getPowetoLogPDF(logx = [3., 4., 5.], alpha = [+2., +3., +4.], logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.")
                    logPDF(1:3) = getPowetoLogPDF(logx = [3., 4., 5.], alpha = [+2., +3., +4.], logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.
    call disp%show("logPDF(1:3)")
    call disp%show( logPDF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        use pm_arrayResize, only: setResized
        real, allocatable :: alpha(:), logPDF(:)
        real :: logMinX, logMaxX, logx(2000)
        integer(IK) :: fileUnit, i
        alpha = [-2.0, -1., 0., +1., +2.]
        call setLinSpace(logx, x1 = log(0.1), x2 = log(10.))
        call setResized(logPDF, size(alpha, 1, IK))
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "getPowetoLogPDF.RK.txt")
        do i = 1, size(logx, 1, IK)
            if (logMinX <= logx(i) .and. logx(i) <= logMaxX) then
                logPDF = getPowetoLogPDF(logx(i), alpha, logMinX, logMaxX)
            else
                logPDF = -huge(0.)
            end if
            write(fileUnit, "(*(g0,:,', '))") exp(logx(i)), exp(logPDF)
        end do
        close(fileUnit)
    end block

end program example