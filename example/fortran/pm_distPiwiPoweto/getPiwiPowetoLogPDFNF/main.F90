program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distUnif, only: getUnifRand
    use pm_distPiwiPoweto, only: getPiwiPowetoLogPDFNF

    implicit none

    integer(IK) , parameter :: NP = 4_IK
    real        , allocatable :: alpha(:), logPDFNF(:), logLimX(:), cumSumArea(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    alpha = [getUnifRand(-5., +5., NP - 1), -3.]
    logLimX = [getLinSpace(log(0.001), log(20.), NP), log(huge(0.))]
    allocate(cumSumArea, mold = logLimX)
    allocate(logPDFNF, mold = alpha)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the PiwiPoweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alpha")
    call disp%show( alpha )
    call disp%show("logLimX")
    call disp%show( logLimX )
    call disp%show("logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX)")
                    logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX)
    call disp%show("logPDFNF")
    call disp%show( logPDFNF )
    call disp%skip()

    call disp%skip()
    call disp%show("alpha")
    call disp%show( alpha )
    call disp%show("logLimX")
    call disp%show( logLimX )
    call disp%show("logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea)")
                    logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea)
    call disp%show("logPDFNF")
    call disp%show( logPDFNF )
    call disp%show("cumSumArea")
    call disp%show( cumSumArea )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the Truncated PiwiPoweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("alpha(1 : NP - 1)")
    call disp%show( alpha(1 : NP - 1) )
    call disp%show("logLimX(1 : NP)")
    call disp%show( logLimX(1 : NP) )
    call disp%show("logPDFNF = getPiwiPowetoLogPDFNF(alpha(1 : NP - 1), logLimX(1 : NP), cumSumArea(1 : NP)) ! `logLimX(NP)` is serves as the upper bound of the support of the distribution.")
                    logPDFNF = getPiwiPowetoLogPDFNF(alpha(1 : NP - 1), logLimX(1 : NP), cumSumArea(1 : NP))
    call disp%show("logPDFNF")
    call disp%show( logPDFNF )
    call disp%show("cumSumArea")
    call disp%show( cumSumArea )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example PDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        alpha = getLinSpace(-5., 5., count = 100_IK)
        logLimX = getLinSpace(log(0.001), log(20.), count = size(alpha, 1, IK) + 1_IK)
        logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX)
        deallocate(cumSumArea); allocate(cumSumArea, mold = logLimX)
        logPDFNF = getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea)
        open(newunit = fileUnit, file = "getPiwiPowetoLogPDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (exp(logLimX(i)), cumSumArea(i), i = 1, size(cumSumArea))
        close(fileUnit)
    end block

end program example