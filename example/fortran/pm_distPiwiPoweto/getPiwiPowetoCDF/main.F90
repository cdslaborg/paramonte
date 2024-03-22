program example

    use pm_kind, only: SK, IK, LK
    use pm_distPiwiPoweto, only: getPiwiPowetoLogPDFNF
    use pm_distPiwiPoweto, only: getPiwiPowetoCDF
    use pm_arraySpace, only: setLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999_IK
    real                    :: logx(NP), CDF(NP), cumSumArea(NP)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of the (truncated) PiwiPoweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("cdf = getPiwiPowetoCDF(logx = log(3.), alpha = [-2.], logLimX = [log(2.), log(huge(0.))]) ! left-truncated PiwiPoweto distribution at [2., +∞)")
                    cdf = getPiwiPowetoCDF(logx = log(3.), alpha = [-2.], logLimX = [log(2.), log(huge(0.))])
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()

    call disp%skip()
    call disp%show("cdf = getPiwiPowetoCDF(logx = log(2.5), alpha = [1., 0.], logLimX = [0.5, 1., 1.5]) ! doubly-truncated PiwiPoweto distribution at `[0.5, 1.5]`.")
                    cdf = getPiwiPowetoCDF(logx = log(2.5), alpha = [1., 0.], logLimX = [0.5, 1., 1.5]) ! doubly-truncated PiwiPoweto distribution at `[0.5, 1.5]`.
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()

    call disp%skip()
    call disp%show("cdf = getPiwiPowetoCDF(logx = log(2.5), alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))]) ! (truncated) PiwiPoweto with Semi-Infinite support.")
                    cdf = getPiwiPowetoCDF(logx = log(2.5), alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))]) ! (truncated) PiwiPoweto with Semi-Infinite support.
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Expedite repeated PDF computations by precomputing the normalization factors.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("cdf = getPiwiPowetoCDF(logx = log(2.5), alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))], logPDFNF = getPiwiPowetoLogPDFNF(alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))], cumSumArea = cumSumArea(1:4)), cumSumArea = cumSumArea(1:4))")
                    cdf = getPiwiPowetoCDF(logx = log(2.5), alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))], logPDFNF = getPiwiPowetoLogPDFNF(alpha = [1., 0., -1.5], logLimX = [0.5, 1., 1.5, log(huge(0.))], cumSumArea = cumSumArea(1:4)), cumSumArea = cumSumArea(1:4))
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example CDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        real, parameter :: LOG_HUGE = log(huge(0.))
        real, allocatable :: cumSumArea(:), alpha(:), logLimX(:)
        integer(IK) :: fileUnit, i
        real :: CDF(4)
        alpha = [3., 1., -1., -5.]
        logLimX = log([2., 5., 10., 15.])
        call setLinSpace(logx, x1 = log(0.001), x2 = log(20.), fopen = .true._LK, lopen = .true._LK)
        if (allocated(cumSumArea)) deallocate(cumSumArea); allocate(cumSumArea(size(logLimX)+2))
        open(newunit = fileUnit, file = "getPiwiPowetoCDF.RK.txt")
        do i = 1, NP
            CDF(1) = getPiwiPowetoCDF(logx(i), [alpha(1:2), 0., alpha(3:4)], [-LOG_HUGE, logLimX(1:4), LOG_HUGE], getPiwiPowetoLogPDFNF([alpha(1:2), 0., alpha(3:4)], [-LOG_HUGE, logLimX(1:4), LOG_HUGE], cumSumArea(1:6)), cumSumArea(1:6)) ! PiwiPoweto
            if (logx(i) > logLimX(1)) then
                CDF(2) = getPiwiPowetoCDF(logx(i), alpha(2:4), [logLimX(1:3), LOG_HUGE], getPiwiPowetoLogPDFNF(alpha(2:4), [logLimX(1:3), LOG_HUGE], cumSumArea(1:4)), cumSumArea(1:4)) ! left-truncated PiwiPoweto
            else
                CDF(2) = 0.
            end if
            if (logx(i) < logLimX(4)) then
                CDF(3) = getPiwiPowetoCDF(logx(i), alpha(2:4), [-LOG_HUGE, logLimX(2:4)], getPiwiPowetoLogPDFNF(alpha(2:4), [-LOG_HUGE, logLimX(2:4)], cumSumArea(1:4)), cumSumArea(1:4)) ! right-truncated PiwiPoweto
            else
                CDF(3) = 1.
            end if
            if (logx(i) > logLimX(1) .and. logx(i) < logLimX(4)) then
                CDF(4) = getPiwiPowetoCDF(logx(i), alpha(1:3), logLimX(1:4), getPiwiPowetoLogPDFNF(alpha(1:3), logLimX(1:4), cumSumArea(1:4)), cumSumArea(1:4)) ! doubly-truncated PiwiPoweto
            elseif (logx(i) <= logLimX(1)) then
                CDF(4) = 0.
            else
                CDF(4) = 1.
            end if
            write(fileUnit,"(*(g0,:,', '))") exp(logx(i)), CDF
        end do
        close(fileUnit)
    end block

end program example