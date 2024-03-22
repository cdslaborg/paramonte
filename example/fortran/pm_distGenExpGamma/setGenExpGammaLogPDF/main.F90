program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: RK => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distGenExpGamma, only: setGenExpGammaLogPDF, getGenExpGammaLogPDFNF

    implicit none

    integer(IK) , parameter     :: NP = 1000_IK
    real(RK)    , allocatable   :: Point(:), logPDF(:), Kappa(:), invOmega(:), LogSigma(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLinSpace(+0.5_RK, +2._RK, count = NP)
    invOmega = getLogSpace(logx1 = log(0.1_RK), logx2 = log(10._RK), count = NP)
    LogSigma = getLinSpace(-3._RK, +3._RK, count = NP)
    Point = getLinSpace(-10._RK, +10._RK, count = NP)
    allocate(logPDF, mold = Point)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of GenExpGamma distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setGenExpGammaLogPDF(logPDF(1), x = 0.5_RK)")
                    call setGenExpGammaLogPDF(logPDF(1), x = 0.5_RK)
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setGenExpGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenExpGammaLogPDFNF(Kappa(1)), kappa = Kappa(1))")
                    call setGenExpGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenExpGammaLogPDFNF(Kappa(1)), kappa = Kappa(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("call setGenExpGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenExpGammaLogPDFNF(Kappa(1), invOmega(1)), kappa = Kappa(1), invOmega = invOmega(1))")
                    call setGenExpGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenExpGammaLogPDFNF(Kappa(1), invOmega(1)), kappa = Kappa(1), invOmega = invOmega(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("LogSigma(1)")
    call disp%show( LogSigma(1) )
    call disp%show("call setGenExpGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenExpGammaLogPDFNF(Kappa(1), invOmega(1)), kappa = Kappa(1), invOmega = invOmega(1), logSigma = LogSigma(1))")
                    call setGenExpGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenExpGammaLogPDFNF(Kappa(1), invOmega(1)), kappa = Kappa(1), invOmega = invOmega(1), logSigma = LogSigma(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at an input vector real value with different parameter values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setGenExpGammaLogPDF(logPDF(1:NP:NP/5), x = 0.5_RK, logPDFNF = getGenExpGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))")
                    call setGenExpGammaLogPDF(logPDF(1:NP:NP/5), x = 0.5_RK, logPDFNF = getGenExpGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setGenExpGammaLogPDF(logPDF(1:NP:NP/5), x = Point(1:NP:NP/5), logPDFNF = getGenExpGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))")
                    call setGenExpGammaLogPDF(logPDF(1:NP:NP/5), x = Point(1:NP:NP/5), logPDFNF = getGenExpGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer     :: fileUnit, i
        real(RK)    :: logPDF(NP,6)
        call setGenExpGammaLogPDF(logPDF(:,1), Point, getGenExpGammaLogPDFNF(.5_RK, 2._RK), .5_RK, 2._RK)
        call setGenExpGammaLogPDF(logPDF(:,2), Point, getGenExpGammaLogPDFNF(1._RK, 2._RK), 1._RK, 2._RK)
        call setGenExpGammaLogPDF(logPDF(:,3), Point, getGenExpGammaLogPDFNF(2._RK, 2._RK), 2._RK, 2._RK)
        call setGenExpGammaLogPDF(logPDF(:,4), Point, getGenExpGammaLogPDFNF(.5_RK, .5_RK), .5_RK, .5_RK)
        call setGenExpGammaLogPDF(logPDF(:,5), Point, getGenExpGammaLogPDFNF(1._RK, .5_RK), 1._RK, .5_RK)
        call setGenExpGammaLogPDF(logPDF(:,6), Point, getGenExpGammaLogPDFNF(2._RK, .5_RK), 2._RK, .5_RK)
        open(newunit = fileUnit, file = "setGenExpGammaLogPDF.RK.txt")
        write(fileUnit,"(7(g0,:,' '))") (Point(i), exp(logPDF(i,:)), i = 1, size(Point))
        close(fileUnit)
    end block

end program example