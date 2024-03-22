program example

    use pm_kind, only: SK, IK, LK
    use pm_distPiwiPoweto, only: getPiwiPowetoLogPDF
    use pm_distPiwiPoweto, only: getPiwiPowetoLogPDFNF
    use pm_arraySpace, only: setLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999_IK
    real                    :: logX(NP), logPDF(NP)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the PiwiPoweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getPiwiPowetoLogPDF(logx = log(3.), alpha = [-2.], logLimX = [log(2.), log(huge(0.))])")
                    logPDF(1) = getPiwiPowetoLogPDF(logx = log(3.), alpha = [-2.], logLimX = [log(2.), log(huge(0.))])
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getPiwiPowetoLogPDF(logx = log(2.5), alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))])")
                    logPDF(1) = getPiwiPowetoLogPDF(logx = log(2.5), alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))])
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getPiwiPowetoLogPDF(logx = log(2.5), alpha = [1., 0.], logLimX = [0.5, 1., 1.5]) ! Truncated PiwiPoweto with the upper bound `1.5`.")
                    logPDF(1) = getPiwiPowetoLogPDF(logx = log(2.5), alpha = [1., 0.], logLimX = [0.5, 1., 1.5]) ! Truncated PiwiPoweto with the upper bound `1.5`.
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expedite repeated PDF computations by precomputing the normalization factors.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getPiwiPowetoLogPDF(logx = log(2.5), alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))], logPDFNF = getPiwiPowetoLogPDFNF(alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))]))")
                    logPDF(1) = getPiwiPowetoLogPDF(logx = log(2.5), alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))], logPDFNF = getPiwiPowetoLogPDFNF(alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))]))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example logPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        real, parameter :: LOG_HUGE = log(huge(0.))
        integer(IK) :: fileUnit, i
        real, allocatable :: logLimX(:), alpha(:)
        real :: logPDF(4)
        alpha = [3., 1., -1., -5.]
        logLimX = log([2., 5., 10., 15.])
        call setLinSpace(logX, x1 = log(0.001), x2 = log(20.), fopen = .true._LK, lopen = .true._LK)
        open(newunit = fileUnit, file = "getPiwiPowetoLogPDF.RK.txt")
        do i = 1, NP
            logPDF = -huge(0.)
            logPDF(1) = getPiwiPowetoLogPDF(logX(i), [alpha(1:2), 0., alpha(3:4)], [-LOG_HUGE, logLimX(1:4), LOG_HUGE]) ! PiwiPoweto
            if (logX(i) > logLimX(1)) logPDF(2) = getPiwiPowetoLogPDF(logX(i), alpha(2:4), [logLimX(1:3), LOG_HUGE]) ! left-truncated PiwiPoweto
            if (logX(i) < logLimX(4)) logPDF(3) = getPiwiPowetoLogPDF(logX(i), alpha(2:4), [-LOG_HUGE, logLimX(2:4)]) ! right-truncated PiwiPoweto
            if (logX(i) > logLimX(1) .and. logX(i) < logLimX(4)) logPDF(4) = getPiwiPowetoLogPDF(logX(i), alpha(1:3), logLimX(1:4)) ! doubly-truncated PiwiPoweto
            write(fileUnit,"(*(g0,:,', '))") exp(logX(i)), exp(logPDF)
        end do
        close(fileUnit)
    end block

end program example