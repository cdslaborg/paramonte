program example

    use pm_kind, only: SK, IK, LK
    use pm_distPareto, only: getParetoLogPDFNF
    use pm_distPareto, only: setParetoLogPDF
    use pm_io, only: display_type

    implicit none

    real                    :: logPDF(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setParetoLogPDF(logPDF(1), logx = 3., alpha = -2., logPDFNF = getParetoLogPDFNF(alpha = -2., logMinX = -2.)) ! Pareto distribution.")
                    call setParetoLogPDF(logPDF(1), logx = 3., alpha = -2., logPDFNF = getParetoLogPDFNF(alpha = -2., logMinX = -2.)) ! Pareto distribution.
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setParetoLogPDF(logPDF(1:3), logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logPDFNF = getParetoLogPDFNF(alpha = -[+2., +3., +4.], logMinX = -2.)) ! Pareto distribution.")
                    call setParetoLogPDF(logPDF(1:3), logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logPDFNF = getParetoLogPDFNF(alpha = -[+2., +3., +4.], logMinX = -2.)) ! Pareto distribution.
    call disp%show("logPDF(1:3)")
    call disp%show( logPDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Truncated Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setParetoLogPDF(logPDF(1), logx = 3., alpha = -2., logPDFNF = getParetoLogPDFNF(alpha = -2., logMinX = -2., logMaxX = 5.)) ! Truncated Pareto distribution.")
                    call setParetoLogPDF(logPDF(1), logx = 3., alpha = -2., logPDFNF = getParetoLogPDFNF(alpha = -2., logMinX = -2., logMaxX = 5.)) ! Truncated Pareto distribution.
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setParetoLogPDF(logPDF(1:3), logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logPDFNF = getParetoLogPDFNF(alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Pareto distribution.")
                    call setParetoLogPDF(logPDF(1:3), logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logPDFNF = getParetoLogPDFNF(alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Pareto distribution.
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
        alpha = -[+0.5, +2.0]
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "setParetoLogPDF.RK.txt")
        do i = 1, size(logx, 1, IK)
            logPDF = -huge(0.)
            if (logMinX <= logx(i)) then
                call setParetoLogPDF(logPDF(1:2), logx(i), alpha, getParetoLogPDFNF(alpha, logMinX))
                if (logx(i) <= logMaxX) call setParetoLogPDF(logPDF(3:4), logx(i), alpha, getParetoLogPDFNF(alpha, logMinX, logMaxX))
            end if
            write(fileUnit, "(*(g0,:,', '))") exp(logx(i)), exp(logPDF)
        end do
        close(fileUnit)
    end block

end program example