program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distPiwiPoweto, only: setPiwiPowetoLogPDF
    use pm_distPiwiPoweto, only: getPiwiPowetoLogPDFNF
    use pm_arraySpace, only: setLinSpace
    use pm_arraySearch, only: getBin

    implicit none

    integer(IK), parameter  :: NP = 2000_IK
    real                    :: logx(NP), logPDF(NP)
    real, allocatable       :: logLimX(:), alpha(:)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("logx(1) = log(2.5)")
                    logx(1) = log(2.5)
    call disp%show("alpha = [1., 0., -1.5]")
                    alpha = [1., 0., -1.5]
    call disp%show("logLimX = log([0.5, 1., 1.5, huge(0.)])")
                    logLimX = log([0.5, 1., 1.5, huge(0.)])
    call disp%show("call setPiwiPowetoLogPDF(logPDF(1), logx(1), alpha, logLimX, logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX))")
                    call setPiwiPowetoLogPDF(logPDF(1), logx(1), alpha, logLimX, logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! If the target `x` value is known a prior to belong to a specific component of the PDF, specify it explicitly to expedite the PDF computation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = log(2.5)")
                    logx(1) = log(2.5)
    call disp%show("alpha = [1., 0., -1.5]")
                    alpha = [1., 0., -1.5]
    call disp%show("logLimX = log([0.5, 1., 1.5, huge(0.)])")
                    logLimX = log([0.5, 1., 1.5, huge(0.)])
    call disp%show("call setPiwiPowetoLogPDF(logPDF(1), logx = logx(1), alpha = alpha, bin = getBin(logLimX, logx(1)), logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX))")
                    call setPiwiPowetoLogPDF(logPDF(1), logx = logx(1), alpha = alpha, bin = getBin(logLimX, logx(1)), logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example logPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        real, parameter :: LOG_HUGE = log(huge(0.))
        integer(IK) :: fileUnit, i
        real :: logPDF(4)
        alpha = [3., 1., -1., -5.]
        logLimX = log([2., 5., 10., 15.])
        call setLinSpace(logx, x1 = log(0.001), x2 = log(20.), fopen = .true._LK, lopen = .true._LK)
        open(newunit = fileUnit, file = "setPiwiPowetoLogPDF.RK.txt")
        do i = 1, NP
            logPDF = -huge(0.)
            call setPiwiPowetoLogPDF(logPDF(1), logx(i), [alpha(1:2), 0., alpha(3:4)], [-LOG_HUGE, logLimX(1:4), LOG_HUGE], getPiwiPowetoLogPDFNF([alpha(1:2), 0., alpha(3:4)], [-LOG_HUGE, logLimX(1:4), LOG_HUGE])) ! PiwiPoweto
            if (logx(i) > logLimX(1)) call setPiwiPowetoLogPDF(logPDF(2), logx(i), alpha(2:4), [logLimX(1:3), LOG_HUGE], getPiwiPowetoLogPDFNF(alpha(2:4), [logLimX(1:3), LOG_HUGE])) ! left-truncated PiwiPoweto
            if (logx(i) < logLimX(4)) call setPiwiPowetoLogPDF(logPDF(3), logx(i), alpha(2:4), [-LOG_HUGE, logLimX(2:4)], getPiwiPowetoLogPDFNF(alpha(2:4), [-LOG_HUGE, logLimX(2:4)])) ! right-truncated PiwiPoweto
            if (logx(i) > logLimX(1) .and. logx(i) < logLimX(4)) call setPiwiPowetoLogPDF(logPDF(4), logx(i), alpha(1:3), logLimX(1:4), getPiwiPowetoLogPDFNF(alpha(1:3), logLimX(1:4))) ! doubly-truncated PiwiPoweto
            write(fileUnit,"(*(g0,:,', '))") exp(logx(i)), exp(logPDF)
        end do
        close(fileUnit)
    end block

end program example