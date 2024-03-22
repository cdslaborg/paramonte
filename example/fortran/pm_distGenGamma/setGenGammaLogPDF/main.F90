program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: RK => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distGenGamma, only: setGenGammaLogPDF, getGenGammaLogPDFNF

    implicit none

    integer(IK) , parameter     :: NP = 1000_IK
    real(RK)    , allocatable   :: Point(:), logPDF(:), Kappa(:), invOmega(:), invSigma(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLogSpace(-2._RK, +2._RK, NP)
    invOmega = getLogSpace(-3._RK, +3._RK, NP)
    invSigma = getLogSpace(-3._RK, +3._RK, NP)
    Point = getLogSpace(log(0.01_RK), log(+10._RK), NP)
    allocate(logPDF, mold = Point)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of GenGamma distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setGenGammaLogPDF(logPDF(1), x = 0.5_RK)")
                    call setGenGammaLogPDF(logPDF(1), x = 0.5_RK)
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setGenGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenGammaLogPDFNF(Kappa(1)), kappa = Kappa(1))")
                    call setGenGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenGammaLogPDFNF(Kappa(1)), kappa = Kappa(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("call setGenGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenGammaLogPDFNF(Kappa(1), invOmega(1)), kappa = Kappa(1), invOmega = invOmega(1))")
                    call setGenGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenGammaLogPDFNF(Kappa(1), invOmega(1)), kappa = Kappa(1), invOmega = invOmega(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("invSigma(1)")
    call disp%show( invSigma(1) )
    call disp%show("call setGenGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenGammaLogPDFNF(Kappa(1), invOmega(1), invSigma(1)), kappa = Kappa(1), invOmega = invOmega(1), invSigma = invSigma(1))")
                    call setGenGammaLogPDF(logPDF(1), x = 0.5_RK, logPDFNF = getGenGammaLogPDFNF(Kappa(1), invOmega(1), invSigma(1)), kappa = Kappa(1), invOmega = invOmega(1), invSigma = invSigma(1))
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
    call disp%show("call setGenGammaLogPDF(logPDF(1:NP:NP/5), x = 0.5_RK, logPDFNF = getGenGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))")
                    call setGenGammaLogPDF(logPDF(1:NP:NP/5), x = 0.5_RK, logPDFNF = getGenGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setGenGammaLogPDF(logPDF(1:NP:NP/5), x = Point(1:NP:NP/5), logPDFNF = getGenGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))")
                    call setGenGammaLogPDF(logPDF(1:NP:NP/5), x = Point(1:NP:NP/5), logPDFNF = getGenGammaLogPDFNF(Kappa(1:NP:NP/5)), kappa = Kappa(1:NP:NP/5))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer     :: fileUnit, i
        real        :: logPDF(NP,6)
        call setGenGammaLogPDF(logPDF(:,1), Point, getGenGammaLogPDFNF(1.0, 0.5, 0.5), 1.0, 0.5, 0.5)
        call setGenGammaLogPDF(logPDF(:,2), Point, getGenGammaLogPDFNF(2.0, 0.5, 1.0), 2.0, 0.5, 1.0)
        call setGenGammaLogPDF(logPDF(:,3), Point, getGenGammaLogPDFNF(0.5, 2.0, 0.5), 0.5, 2.0, 0.5)
        call setGenGammaLogPDF(logPDF(:,4), Point, getGenGammaLogPDFNF(0.2, 5.0, 0.2), 0.2, 5.0, 0.2)
        call setGenGammaLogPDF(logPDF(:,5), Point, getGenGammaLogPDFNF(.14, 7.0, .14), .14, 7.0, .14)
        call setGenGammaLogPDF(logPDF(:,6), Point, getGenGammaLogPDFNF(2.0, 5.0, 0.3), 2.0, 5.0, 0.3)
        open(newunit = fileUnit, file = "setGenGammaLogPDF.RK.txt")
        write(fileUnit,"(7(g0,:,' '))") (Point(i), exp(logPDF(i,:)), i = 1, size(Point))
        close(fileUnit)
    end block

end program example