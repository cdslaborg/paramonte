program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distUnif, only: getUnifRand
    use pm_distPoweto, only: getPowetoLogPDFNF

    implicit none

    integer(IK) , parameter :: NP = 4_IK
    real        , allocatable :: Alpha(:), LogNormFac(:), LogLimX(:), CumSumArea(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Alpha = [getUnifRand(-5., +5., NP - 1), -3.]
    LogLimX = [getLinSpace(log(0.001), log(20.), NP), log(huge(0.))]
    allocate(CumSumArea, mold = LogLimX)
    allocate(LogNormFac, mold = Alpha)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Alpha")
    call disp%show( Alpha )
    call disp%show("LogLimX")
    call disp%show( LogLimX )
    call disp%show("LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX)")
                    LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX)
    call disp%show("LogNormFac")
    call disp%show( LogNormFac )
    call disp%skip()

    call disp%skip()
    call disp%show("Alpha")
    call disp%show( Alpha )
    call disp%show("LogLimX")
    call disp%show( LogLimX )
    call disp%show("LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea)")
                    LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea)
    call disp%show("LogNormFac")
    call disp%show( LogNormFac )
    call disp%show("CumSumArea")
    call disp%show( CumSumArea )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the Truncated Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Alpha(1 : NP - 1)")
    call disp%show( Alpha(1 : NP - 1) )
    call disp%show("LogLimX(1 : NP)")
    call disp%show( LogLimX(1 : NP) )
    call disp%show("LogNormFac = getPowetoLogPDFNF(Alpha(1 : NP - 1), LogLimX(1 : NP), CumSumArea(1 : NP)) ! `LogLimX(NP)` is serves as the upper bound of the support of the distribution.")
                    LogNormFac = getPowetoLogPDFNF(Alpha(1 : NP - 1), LogLimX(1 : NP), CumSumArea(1 : NP))
    call disp%show("LogNormFac")
    call disp%show( LogNormFac )
    call disp%show("CumSumArea")
    call disp%show( CumSumArea )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example PDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        Alpha = getLinSpace(-5., 5., count = 100_IK)
        LogLimX = getLinSpace(log(0.001), log(20.), count = size(Alpha, 1, IK) + 1_IK)
        LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX)
        deallocate(CumSumArea); allocate(CumSumArea, mold = LogLimX)
        LogNormFac = getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea)
        open(newunit = fileUnit, file = "getPowetoLogPDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (exp(LogLimX(i)), CumSumArea(i), i = 1, size(CumSumArea))
        close(fileUnit)
    end block

end program example