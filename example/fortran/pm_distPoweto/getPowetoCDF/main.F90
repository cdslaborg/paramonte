program example

    use pm_kind, only: SK, IK, LK
    use pm_distPoweto, only: getPowetoLogPDFNF
    use pm_distPoweto, only: getPowetoCDF
    use pm_arraySpace, only: setLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999_IK
    real                    :: LogX(NP), CDF(NP), CumSumArea(NP)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of the (truncated) Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("cdf = getPowetoCDF(logx = log(3.), Alpha = [-2.], LogLimX = [log(2.), log(huge(0.))]) ! left-truncated Poweto distribution at [2., +∞)")
                    cdf = getPowetoCDF(logx = log(3.), Alpha = [-2.], LogLimX = [log(2.), log(huge(0.))])
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()

    call disp%skip()
    call disp%show("cdf = getPowetoCDF(logx = log(2.5), Alpha = [1., 0.], LogLimX = [0.5, 1., 1.5]) ! doubly-truncated Poweto distribution at `[0.5, 1.5]`.")
                    cdf = getPowetoCDF(logx = log(2.5), Alpha = [1., 0.], LogLimX = [0.5, 1., 1.5]) ! doubly-truncated Poweto distribution at `[0.5, 1.5]`.
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()

    call disp%skip()
    call disp%show("cdf = getPowetoCDF(logx = log(2.5), Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))]) ! (truncated) Poweto with Semi-Infinite support.")
                    cdf = getPowetoCDF(logx = log(2.5), Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))]) ! (truncated) Poweto with Semi-Infinite support.
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expedite repeated PDF computations by precomputing the normalization factors.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("cdf = getPowetoCDF(logx = log(2.5), Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))], LogNormFac = getPowetoLogPDFNF(Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))], CumSumArea = CumSumArea(1:4)), CumSumArea = CumSumArea(1:4))")
                    cdf = getPowetoCDF(logx = log(2.5), Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))], LogNormFac = getPowetoLogPDFNF(Alpha = [1., 0., -1.5], LogLimX = [0.5, 1., 1.5, log(huge(0.))], CumSumArea = CumSumArea(1:4)), CumSumArea = CumSumArea(1:4))
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example CDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        real, parameter :: LOG_HUGE = log(huge(0.))
        real, allocatable :: CumSumArea(:), Alpha(:), LogLimX(:)
        integer(IK) :: fileUnit, i
        real :: CDF(4)
        Alpha = [3., 1., -1., -5.]
        LogLimX = log([2., 5., 10., 15.])
        call setLinSpace(LogX, x1 = log(0.001), x2 = log(20.), fopen = .true._LK, lopen = .true._LK)
        if (allocated(CumSumArea)) deallocate(CumSumArea); allocate(CumSumArea(size(LogLimX)+2))
        open(newunit = fileUnit, file = "getPowetoCDF.RK.txt")
        do i = 1, NP
            CDF(1) = getPowetoCDF(LogX(i), [Alpha(1:2), 0., Alpha(3:4)], [-LOG_HUGE, LogLimX(1:4), LOG_HUGE], getPowetoLogPDFNF([Alpha(1:2), 0., Alpha(3:4)], [-LOG_HUGE, LogLimX(1:4), LOG_HUGE], CumSumArea(1:6)), CumSumArea(1:6)) ! Poweto
            if (LogX(i) > LogLimX(1)) then
                CDF(2) = getPowetoCDF(LogX(i), Alpha(2:4), [LogLimX(1:3), LOG_HUGE], getPowetoLogPDFNF(Alpha(2:4), [LogLimX(1:3), LOG_HUGE], CumSumArea(1:4)), CumSumArea(1:4)) ! left-truncated Poweto
            else
                CDF(2) = 0.
            end if
            if (LogX(i) < LogLimX(4)) then
                CDF(3) = getPowetoCDF(LogX(i), Alpha(2:4), [-LOG_HUGE, LogLimX(2:4)], getPowetoLogPDFNF(Alpha(2:4), [-LOG_HUGE, LogLimX(2:4)], CumSumArea(1:4)), CumSumArea(1:4)) ! right-truncated Poweto
            else
                CDF(3) = 1.
            end if
            if (LogX(i) > LogLimX(1) .and. LogX(i) < LogLimX(4)) then
                CDF(4) = getPowetoCDF(LogX(i), Alpha(1:3), LogLimX(1:4), getPowetoLogPDFNF(Alpha(1:3), LogLimX(1:4), CumSumArea(1:4)), CumSumArea(1:4)) ! doubly-truncated Poweto
            elseif (LogX(i) <= LogLimX(1)) then
                CDF(4) = 0.
            else
                CDF(4) = 1.
            end if
            write(fileUnit,"(*(g0,:,', '))") exp(LogX(i)), CDF
        end do
        close(fileUnit)
    end block

end program example