program example

    use pm_kind, only: SK, IK, LK
    use pm_distPoweto, only: setPowetoLogPDF
    use pm_distPoweto, only: getPowetoLogPDFNF
    use pm_arraySpace, only: setLinSpace
    use pm_io, only: display_type
    use pm_arraySearch, only: getBin

    implicit none

    integer(IK), parameter  :: NP = 2000_IK
    real                    :: LogX(NP), LogPDF(NP)
    real, allocatable       :: LogLimX(:), Alpha(:)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("LogX(1) = log(2.5)")
                    LogX(1) = log(2.5)
    call disp%show("Alpha = [1., 0., -1.5]")
                    Alpha = [1., 0., -1.5]
    call disp%show("LogLimX = log([0.5, 1., 1.5, huge(0.)])")
                    LogLimX = log([0.5, 1., 1.5, huge(0.)])
    call disp%show("call setPowetoLogPDF(LogPDF(1), LogX(1), Alpha, LogLimX, LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX))")
                    call setPowetoLogPDF(LogPDF(1), LogX(1), Alpha, LogLimX, LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX))
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! If the target `x` value is known a prior to belong to a specific component of the PDF, specify it explicitly to expedite the PDF computation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1) = log(2.5)")
                    LogX(1) = log(2.5)
    call disp%show("Alpha = [1., 0., -1.5]")
                    Alpha = [1., 0., -1.5]
    call disp%show("LogLimX = log([0.5, 1., 1.5, huge(0.)])")
                    LogLimX = log([0.5, 1., 1.5, huge(0.)])
    call disp%show("call setPowetoLogPDF(LogPDF(1), logx = LogX(1), Alpha = Alpha, bin = getBin(LogLimX, LogX(1)), LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX))")
                    call setPowetoLogPDF(LogPDF(1), logx = LogX(1), Alpha = Alpha, bin = getBin(LogLimX, LogX(1)), LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX))
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example LogPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        real, parameter :: LOG_HUGE = log(huge(0.))
        integer(IK) :: fileUnit, i
        real :: LogPDF(4)
        Alpha = [3., 1., -1., -5.]
        LogLimX = log([2., 5., 10., 15.])
        call setLinSpace(LogX, x1 = log(0.001), x2 = log(20.), fopen = .true._LK, lopen = .true._LK)
        open(newunit = fileUnit, file = "setPowetoLogPDF.RK.txt")
        do i = 1, NP
            LogPDF = -huge(0.)
            call setPowetoLogPDF(LogPDF(1), LogX(i), [Alpha(1:2), 0., Alpha(3:4)], [-LOG_HUGE, LogLimX(1:4), LOG_HUGE], getPowetoLogPDFNF([Alpha(1:2), 0., Alpha(3:4)], [-LOG_HUGE, LogLimX(1:4), LOG_HUGE])) ! Poweto
            if (LogX(i) > LogLimX(1)) call setPowetoLogPDF(LogPDF(2), LogX(i), Alpha(2:4), [LogLimX(1:3), LOG_HUGE], getPowetoLogPDFNF(Alpha(2:4), [LogLimX(1:3), LOG_HUGE])) ! left-truncated Poweto
            if (LogX(i) < LogLimX(4)) call setPowetoLogPDF(LogPDF(3), LogX(i), Alpha(2:4), [-LOG_HUGE, LogLimX(2:4)], getPowetoLogPDFNF(Alpha(2:4), [-LOG_HUGE, LogLimX(2:4)])) ! right-truncated Poweto
            if (LogX(i) > LogLimX(1) .and. LogX(i) < LogLimX(4)) call setPowetoLogPDF(LogPDF(4), LogX(i), Alpha(1:3), LogLimX(1:4), getPowetoLogPDFNF(Alpha(1:3), LogLimX(1:4))) ! doubly-truncated Poweto
            write(fileUnit,"(*(g0,:,', '))") exp(LogX(i)), exp(LogPDF)
        end do
        close(fileUnit)
    end block

end program example