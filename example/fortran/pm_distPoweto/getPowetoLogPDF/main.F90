program example

    use pm_kind, only: SK, IK, LK
    use pm_distPoweto, only: getPowetoLogPDF
    use pm_distPoweto, only: getPowetoLogPDFNF
    use pm_arraySpace, only: setLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999_IK
    real                    :: LogX(NP), LogPDF(NP)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1) = getPowetoLogPDF(logx = log(3.), Alpha = [-2.], LogLimX = [log(2.), log(huge(0.))])")
                    LogPDF(1) = getPowetoLogPDF(logx = log(3.), Alpha = [-2.], LogLimX = [log(2.), log(huge(0.))])
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1) = getPowetoLogPDF(logx = log(2.5), Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))])")
                    LogPDF(1) = getPowetoLogPDF(logx = log(2.5), Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))])
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1) = getPowetoLogPDF(logx = log(2.5), Alpha = [1., 0.], LogLimX = [0.5, 1., 1.5]) ! Truncated Poweto with the upper bound `1.5`.")
                    LogPDF(1) = getPowetoLogPDF(logx = log(2.5), Alpha = [1., 0.], LogLimX = [0.5, 1., 1.5]) ! Truncated Poweto with the upper bound `1.5`.
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expedite repeated PDF computations by precomputing the normalization factors.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDF(1) = getPowetoLogPDF(logx = log(2.5), Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))], LogNormFac = getPowetoLogPDFNF(Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))]))")
                    LogPDF(1) = getPowetoLogPDF(logx = log(2.5), Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))], LogNormFac = getPowetoLogPDFNF(Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))]))
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example LogPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        real, parameter :: LOG_HUGE = log(huge(0.))
        integer(IK) :: fileUnit, i
        real, allocatable :: LogLimX(:), Alpha(:)
        real :: LogPDF(4)
        Alpha = [3., 1., -1., -5.]
        LogLimX = log([2., 5., 10., 15.])
        call setLinSpace(LogX, x1 = log(0.001), x2 = log(20.), fopen = .true._LK, lopen = .true._LK)
        open(newunit = fileUnit, file = "getPowetoLogPDF.RK.txt")
        do i = 1, NP
            LogPDF = -huge(0.)
            LogPDF(1) = getPowetoLogPDF(LogX(i), [Alpha(1:2), 0., Alpha(3:4)], [-LOG_HUGE, LogLimX(1:4), LOG_HUGE]) ! Poweto
            if (LogX(i) > LogLimX(1)) LogPDF(2) = getPowetoLogPDF(LogX(i), Alpha(2:4), [LogLimX(1:3), LOG_HUGE]) ! left-truncated Poweto
            if (LogX(i) < LogLimX(4)) LogPDF(3) = getPowetoLogPDF(LogX(i), Alpha(2:4), [-LOG_HUGE, LogLimX(2:4)]) ! right-truncated Poweto
            if (LogX(i) > LogLimX(1) .and. LogX(i) < LogLimX(4)) LogPDF(4) = getPowetoLogPDF(LogX(i), Alpha(1:3), LogLimX(1:4)) ! doubly-truncated Poweto
            write(fileUnit,"(*(g0,:,', '))") exp(LogX(i)), exp(LogPDF)
        end do
        close(fileUnit)
    end block

end program example