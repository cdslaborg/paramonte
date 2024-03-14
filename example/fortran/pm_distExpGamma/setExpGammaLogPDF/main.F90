program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: RK => RK32 ! all other real kinds are also acceptable: RK, RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distExpGamma, only: setExpGammaLogPDF, getExpGammaLogPDFNF

    implicit none

    integer(IK) , parameter     :: NP = 1000_IK
    real(RK)    , allocatable   :: Point(:), LogPDF(:), Kappa(:), InvOmega(:), LogSigma(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLinSpace(+0.5_RK, +2._RK, count = NP)
    InvOmega = getLogSpace(logx1 = log(0.1_RK), logx2 = log(10._RK), count = NP)
    LogSigma = getLinSpace(-3._RK, +3._RK, count = NP)
    Point = getLinSpace(-10._RK, +10._RK, count = NP)
    allocate(LogPDF, mold = Point)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of ExpGamma distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setExpGammaLogPDF(LogPDF(1), x = 0.5_RK)")
                    call setExpGammaLogPDF(LogPDF(1), x = 0.5_RK)
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setExpGammaLogPDF(LogPDF(1), x = 0.5_RK, logNormFac = getExpGammaLogPDFNF(Kappa(1)), kappa = Kappa(1))")
                    call setExpGammaLogPDF(LogPDF(1), x = 0.5_RK, logNormFac = getExpGammaLogPDFNF(Kappa(1)), kappa = Kappa(1))
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("InvOmega(1)")
    call disp%show( InvOmega(1) )
    call disp%show("call setExpGammaLogPDF(LogPDF(1), x = 0.5_RK, logNormFac = getExpGammaLogPDFNF(Kappa(1)), kappa = Kappa(1))")
                    call setExpGammaLogPDF(LogPDF(1), x = 0.5_RK, logNormFac = getExpGammaLogPDFNF(Kappa(1)), kappa = Kappa(1))
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("InvOmega(1)")
    call disp%show( InvOmega(1) )
    call disp%show("LogSigma(1)")
    call disp%show( LogSigma(1) )
    call disp%show("call setExpGammaLogPDF(LogPDF(1), x = 0.5_RK, logNormFac = getExpGammaLogPDFNF(Kappa(1)), kappa = Kappa(1), logSigma = LogSigma(1))")
                    call setExpGammaLogPDF(LogPDF(1), x = 0.5_RK, logNormFac = getExpGammaLogPDFNF(Kappa(1)), kappa = Kappa(1), logSigma = LogSigma(1))
    call disp%show("LogPDF(1)")
    call disp%show( LogPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at an input vector real value with different parameter values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setExpGammaLogPDF(LogPDF(1:NP:NP/5), x = 0.5_RK, logNormFac = getExpGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))")
                    call setExpGammaLogPDF(LogPDF(1:NP:NP/5), x = 0.5_RK, logNormFac = getExpGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))
    call disp%show("LogPDF(1:NP:NP/5)")
    call disp%show( LogPDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setExpGammaLogPDF(LogPDF(1:NP:NP/5), x = Point(1:NP:NP/5), logNormFac = getExpGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))")
                    call setExpGammaLogPDF(LogPDF(1:NP:NP/5), x = Point(1:NP:NP/5), logNormFac = getExpGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))
    call disp%show("LogPDF(1:NP:NP/5)")
    call disp%show( LogPDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example PDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer     :: fileUnit, i
        real(RK)    :: LogPDF(NP,6)
        call setExpGammaLogPDF(LogPDF(:,1), Point, getExpGammaLogPDFNF(.5_RK), .5_RK, -.5)
        call setExpGammaLogPDF(LogPDF(:,2), Point, getExpGammaLogPDFNF(1._RK), 1._RK, -.5)
        call setExpGammaLogPDF(LogPDF(:,3), Point, getExpGammaLogPDFNF(2._RK), 2._RK, -.5)
        call setExpGammaLogPDF(LogPDF(:,4), Point, getExpGammaLogPDFNF(.5_RK), .5_RK, 2.0)
        call setExpGammaLogPDF(LogPDF(:,5), Point, getExpGammaLogPDFNF(1._RK), 1._RK, 2.0)
        call setExpGammaLogPDF(LogPDF(:,6), Point, getExpGammaLogPDFNF(2._RK), 2._RK, 2.0)
        open(newunit = fileUnit, file = "setExpGammaLogPDF.RK.txt")
        write(fileUnit,"(7(g0,:,' '))") (Point(i), exp(LogPDF(i,:)), i = 1, size(Point))
        close(fileUnit)
    end block

end program example