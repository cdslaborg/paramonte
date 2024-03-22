program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLinSpace
    use pm_distGamma, only: getGammaLogPDFNF

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real        , allocatable :: invSigma(:), logPDFNF(:), Kappa(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLinSpace(0.01, 10., count = NP)
    invSigma = getLinSpace(0.01, 10., count = NP)
    allocate(logPDFNF, mold = Kappa)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the Gamma distribution PDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("[Kappa(1), invSigma(1)]")
    call disp%show( [Kappa(1), invSigma(1)] )
    call disp%show("logPDFNF(1) = getGammaLogPDFNF(kappa = Kappa(1), invSigma = invSigma(1))")
                    logPDFNF(1) = getGammaLogPDFNF(kappa = Kappa(1), invSigma = invSigma(1))
    call disp%show("logPDFNF(1)")
    call disp%show( logPDFNF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at at a mix of scalar and vector input values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("[Kappa(1), invSigma(1)]")
    call disp%show( [Kappa(1), invSigma(1)] )
    call disp%show("logPDFNF(1) = getGammaLogPDFNF(kappa = Kappa(1), invSigma = invSigma(1))")
                    logPDFNF(1) = getGammaLogPDFNF(kappa = Kappa(1), invSigma = invSigma(1))
    call disp%show("logPDFNF(1)")
    call disp%show( logPDFNF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("logPDFNF(1:NP:NP/5) = getGammaLogPDFNF(kappa = Kappa(1:NP:NP/5), invSigma = 1.)")
                    logPDFNF(1:NP:NP/5) = getGammaLogPDFNF(kappa = Kappa(1:NP:NP/5), invSigma = 1.)
    call disp%show("logPDFNF(1:NP:NP/5)")
    call disp%show( logPDFNF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("logPDFNF(1:NP:NP/5) = getGammaLogPDFNF(kappa = 1., invSigma = invSigma(1:NP:NP/5))")
                    logPDFNF(1:NP:NP/5) = getGammaLogPDFNF(kappa = 1., invSigma = invSigma(1:NP:NP/5))
    call disp%show("logPDFNF(1:NP:NP/5)")
    call disp%show( logPDFNF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("logPDFNF(1:NP:NP/5) = getGammaLogPDFNF(kappa = Kappa(1:NP:NP/5), invSigma = invSigma(1:NP:NP/5))")
                    logPDFNF(1:NP:NP/5) = getGammaLogPDFNF(kappa = Kappa(1:NP:NP/5), invSigma = invSigma(1:NP:NP/5))
    call disp%show("logPDFNF(1:NP:NP/5)")
    call disp%show( logPDFNF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        logPDFNF = getGammaLogPDFNF(Kappa, invSigma)
        open(newunit = fileUnit, file = "getGammaLogPDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (Kappa(i), exp(logPDFNF(i)), i = 1, size(logPDFNF))
        close(fileUnit)
    end block

end program example