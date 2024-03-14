program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLinSpace
    use pm_distGamma, only: getGammaLogPDFNF

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real        , allocatable :: InvSigma(:), LogNormFac(:), Kappa(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLinSpace(0.01, 10., count = NP)
    InvSigma = getLinSpace(0.01, 10., count = NP)
    allocate(LogNormFac, mold = Kappa)

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
    call disp%show("[Kappa(1), InvSigma(1)]")
    call disp%show( [Kappa(1), InvSigma(1)] )
    call disp%show("LogNormFac(1) = getGammaLogPDFNF(kappa = Kappa(1), invSigma = InvSigma(1))")
                    LogNormFac(1) = getGammaLogPDFNF(kappa = Kappa(1), invSigma = InvSigma(1))
    call disp%show("LogNormFac(1)")
    call disp%show( LogNormFac(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at at a mix of scalar and vector input values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("[Kappa(1), InvSigma(1)]")
    call disp%show( [Kappa(1), InvSigma(1)] )
    call disp%show("LogNormFac(1) = getGammaLogPDFNF(kappa = Kappa(1), invSigma = InvSigma(1))")
                    LogNormFac(1) = getGammaLogPDFNF(kappa = Kappa(1), invSigma = InvSigma(1))
    call disp%show("LogNormFac(1)")
    call disp%show( LogNormFac(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("LogNormFac(1:NP:NP/5) = getGammaLogPDFNF(kappa = Kappa(1:NP:NP/5), invSigma = 1.)")
                    LogNormFac(1:NP:NP/5) = getGammaLogPDFNF(kappa = Kappa(1:NP:NP/5), invSigma = 1.)
    call disp%show("LogNormFac(1:NP:NP/5)")
    call disp%show( LogNormFac(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("LogNormFac(1:NP:NP/5) = getGammaLogPDFNF(kappa = 1., invSigma = InvSigma(1:NP:NP/5))")
                    LogNormFac(1:NP:NP/5) = getGammaLogPDFNF(kappa = 1., invSigma = InvSigma(1:NP:NP/5))
    call disp%show("LogNormFac(1:NP:NP/5)")
    call disp%show( LogNormFac(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("LogNormFac(1:NP:NP/5) = getGammaLogPDFNF(kappa = Kappa(1:NP:NP/5), invSigma = InvSigma(1:NP:NP/5))")
                    LogNormFac(1:NP:NP/5) = getGammaLogPDFNF(kappa = Kappa(1:NP:NP/5), invSigma = InvSigma(1:NP:NP/5))
    call disp%show("LogNormFac(1:NP:NP/5)")
    call disp%show( LogNormFac(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example PDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        LogNormFac = getGammaLogPDFNF(Kappa, InvSigma)
        open(newunit = fileUnit, file = "getGammaLogPDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (Kappa(i), exp(LogNormFac(i)), i = 1, size(LogNormFac))
        close(fileUnit)
    end block

end program example