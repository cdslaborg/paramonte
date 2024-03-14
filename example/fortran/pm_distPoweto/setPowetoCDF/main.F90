program example

    use pm_kind, only: SK, IK, LK
    use pm_distPoweto, only: getPowetoLogPDFNF
    use pm_distPoweto, only: setPowetoCDF
    use pm_arraySearch, only: getBin
    use pm_arraySpace, only: setLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999_IK
    real                    :: LogX(NP), CDF(NP)
    real, allocatable       :: LogLimX(:), Alpha(:), CumSumArea(:)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("LogX(1) = log(2.5)")
                    LogX(1) = log(2.5)
    call disp%show("Alpha = [1., 0., -1.5]")
                    Alpha = [1., 0., -1.5]
    call disp%show("LogLimX = log([0.5, 1., 1.5, huge(0.)])")
                    LogLimX = log([0.5, 1., 1.5, huge(0.)])
    call disp%show("if (allocated(CumSumArea)) deallocate(CumSumArea); allocate(CumSumArea, mold = LogLimX)")
                    if (allocated(CumSumArea)) deallocate(CumSumArea); allocate(CumSumArea, mold = LogLimX)
    call disp%show("call setPowetoCDF(CDF(1), LogX(1), Alpha, LogLimX, LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea), CumSumArea = CumSumArea)")
                    call setPowetoCDF(CDF(1), LogX(1), Alpha, LogLimX, LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea), CumSumArea = CumSumArea)
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! If the target `x` value is known a prior to belong to a specific component of the CDF, specify it explicitly to expedite the CDF computation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1) = log(2.5)")
                    LogX(1) = log(2.5)
    call disp%show("Alpha = [1., 0., -1.5]")
                    Alpha = [1., 0., -1.5]
    call disp%show("LogLimX = log([0.5, 1., 1.5, huge(0.)])")
                    LogLimX = log([0.5, 1., 1.5, huge(0.)])
    call disp%show("if (allocated(CumSumArea)) deallocate(CumSumArea); allocate(CumSumArea, mold = LogLimX)")
                    if (allocated(CumSumArea)) deallocate(CumSumArea); allocate(CumSumArea, mold = LogLimX)
    call disp%show("call setPowetoCDF(CDF(1), LogX(1), Alpha, LogLimX, LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea), CumSumArea = CumSumArea, bin = getBin(LogLimX, LogX(1)))")
                    call setPowetoCDF(CDF(1), LogX(1), Alpha, LogLimX, LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea), CumSumArea = CumSumArea, bin = getBin(LogLimX, LogX(1)))
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
        Alpha = [3., 1., -1., -5.]
        LogLimX = log([2., 5., 10., 15.])
        call setLinSpace(LogX, x1 = log(0.001), x2 = log(20.), fopen = .true._LK, lopen = .true._LK)
        if (allocated(CumSumArea)) deallocate(CumSumArea); allocate(CumSumArea(size(LogLimX)+2))
        open(newunit = fileUnit, file = "setPowetoCDF.RK.txt")
        do i = 1, NP
            call setPowetoCDF(CDF(1), LogX(i), [Alpha(1:2), 0., Alpha(3:4)], [-LOG_HUGE, LogLimX(1:4), LOG_HUGE], getPowetoLogPDFNF([Alpha(1:2), 0., Alpha(3:4)], [-LOG_HUGE, LogLimX(1:4), LOG_HUGE], CumSumArea(1:6)), CumSumArea(1:6)) ! Poweto
            if (LogX(i) > LogLimX(1)) then
                call setPowetoCDF(CDF(2), LogX(i), Alpha(2:4), [LogLimX(1:3), LOG_HUGE], getPowetoLogPDFNF(Alpha(2:4), [LogLimX(1:3), LOG_HUGE], CumSumArea(1:4)), CumSumArea(1:4)) ! left-truncated Poweto
            else
                CDF(2) = 0.
            end if
            if (LogX(i) < LogLimX(4)) then
                call setPowetoCDF(CDF(3), LogX(i), Alpha(2:4), [-LOG_HUGE, LogLimX(2:4)], getPowetoLogPDFNF(Alpha(2:4), [-LOG_HUGE, LogLimX(2:4)], CumSumArea(1:4)), CumSumArea(1:4)) ! right-truncated Poweto
            else
                CDF(3) = 1.
            end if
            if (LogX(i) > LogLimX(1) .and. LogX(i) < LogLimX(4)) then
                call setPowetoCDF(CDF(4), LogX(i), Alpha(1:3), LogLimX(1:4), getPowetoLogPDFNF(Alpha(1:3), LogLimX(1:4), CumSumArea(1:4)), CumSumArea(1:4)) ! doubly-truncated Poweto
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