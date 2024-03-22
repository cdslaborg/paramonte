program example

    use pm_kind, only: SK, IK, LK
    use pm_distPiwiPoweto, only: getPiwiPowetoLogPDFNF
    use pm_distPiwiPoweto, only: setPiwiPowetoCDF
    use pm_arraySearch, only: getBin
    use pm_arraySpace, only: setLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999_IK
    real                    :: logx(NP), CDF(NP)
    real, allocatable       :: logLimX(:), alpha(:), cumSumArea(:)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("logx(1) = log(2.5)")
                    logx(1) = log(2.5)
    call disp%show("alpha = [1., 0., -1.5]")
                    alpha = [1., 0., -1.5]
    call disp%show("logLimX = log([0.5, 1., 1.5, huge(0.)])")
                    logLimX = log([0.5, 1., 1.5, huge(0.)])
    call disp%show("if (allocated(cumSumArea)) deallocate(cumSumArea); allocate(cumSumArea, mold = logLimX)")
                    if (allocated(cumSumArea)) deallocate(cumSumArea); allocate(cumSumArea, mold = logLimX)
    call disp%show("call setPiwiPowetoCDF(CDF(1), logx(1), alpha, logLimX, logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea), cumSumArea = cumSumArea)")
                    call setPiwiPowetoCDF(CDF(1), logx(1), alpha, logLimX, logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea), cumSumArea = cumSumArea)
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! If the target `x` value is known a prior to belong to a specific component of the CDF, specify it explicitly to expedite the CDF computation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = log(2.5)")
                    logx(1) = log(2.5)
    call disp%show("alpha = [1., 0., -1.5]")
                    alpha = [1., 0., -1.5]
    call disp%show("logLimX = log([0.5, 1., 1.5, huge(0.)])")
                    logLimX = log([0.5, 1., 1.5, huge(0.)])
    call disp%show("if (allocated(cumSumArea)) deallocate(cumSumArea); allocate(cumSumArea, mold = logLimX)")
                    if (allocated(cumSumArea)) deallocate(cumSumArea); allocate(cumSumArea, mold = logLimX)
    call disp%show("call setPiwiPowetoCDF(CDF(1), logx(1), alpha, logLimX, logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea), cumSumArea = cumSumArea, bin = getBin(logLimX, logx(1)))")
                    call setPiwiPowetoCDF(CDF(1), logx(1), alpha, logLimX, logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea), cumSumArea = cumSumArea, bin = getBin(logLimX, logx(1)))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example CDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        real, parameter :: LOG_HUGE = log(huge(0.))
        integer(IK) :: fileUnit, i
        real :: CDF(4)
        alpha = [3., 1., -1., -5.]
        logLimX = log([2., 5., 10., 15.])
        call setLinSpace(logx, x1 = log(0.001), x2 = log(20.), fopen = .true._LK, lopen = .true._LK)
        if (allocated(cumSumArea)) deallocate(cumSumArea); allocate(cumSumArea(size(logLimX)+2))
        open(newunit = fileUnit, file = "setPiwiPowetoCDF.RK.txt")
        do i = 1, NP
            call setPiwiPowetoCDF(CDF(1), logx(i), [alpha(1:2), 0., alpha(3:4)], [-LOG_HUGE, logLimX(1:4), LOG_HUGE], getPiwiPowetoLogPDFNF([alpha(1:2), 0., alpha(3:4)], [-LOG_HUGE, logLimX(1:4), LOG_HUGE], cumSumArea(1:6)), cumSumArea(1:6)) ! PiwiPoweto
            if (logx(i) > logLimX(1)) then
                call setPiwiPowetoCDF(CDF(2), logx(i), alpha(2:4), [logLimX(1:3), LOG_HUGE], getPiwiPowetoLogPDFNF(alpha(2:4), [logLimX(1:3), LOG_HUGE], cumSumArea(1:4)), cumSumArea(1:4)) ! left-truncated PiwiPoweto
            else
                CDF(2) = 0.
            end if
            if (logx(i) < logLimX(4)) then
                call setPiwiPowetoCDF(CDF(3), logx(i), alpha(2:4), [-LOG_HUGE, logLimX(2:4)], getPiwiPowetoLogPDFNF(alpha(2:4), [-LOG_HUGE, logLimX(2:4)], cumSumArea(1:4)), cumSumArea(1:4)) ! right-truncated PiwiPoweto
            else
                CDF(3) = 1.
            end if
            if (logx(i) > logLimX(1) .and. logx(i) < logLimX(4)) then
                call setPiwiPowetoCDF(CDF(4), logx(i), alpha(1:3), logLimX(1:4), getPiwiPowetoLogPDFNF(alpha(1:3), logLimX(1:4), cumSumArea(1:4)), cumSumArea(1:4)) ! doubly-truncated PiwiPoweto
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