program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_kind, only: RK => RK32 ! all other real kinds are also supported, e.g.: RK, RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLinSpace
    use pm_distGenGamma, only: getGenGammaLogPDFNF

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RK)    , allocatable :: InvSigma(:), InvOmega(:), Kappa(:), LogNormFac(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLinSpace(0.01_RK, 10._RK, count = NP)
    InvOmega = getLinSpace(0.1_RK, 10._RK, count = NP)
    InvSigma = getLinSpace(0.1_RK, 10._RK, count = NP)
    allocate(LogNormFac, mold = Kappa)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the GenGamma distribution PDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("LogNormFac(1) = getGenGammaLogPDFNF(kappa = Kappa(1))")
                    LogNormFac(1) = getGenGammaLogPDFNF(kappa = Kappa(1))
    call disp%show("LogNormFac(1)")
    call disp%show( LogNormFac(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("InvOmega(1)")
    call disp%show( InvOmega(1) )
    call disp%show("LogNormFac(1) = getGenGammaLogPDFNF(Kappa(1), InvOmega(1))")
                    LogNormFac(1) = getGenGammaLogPDFNF(Kappa(1), InvOmega(1))
    call disp%show("LogNormFac(1)")
    call disp%show( LogNormFac(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("InvOmega(1)")
    call disp%show( InvOmega(1) )
    call disp%show("InvSigma(1)")
    call disp%show( InvSigma(1) )
    call disp%show("LogNormFac(1) = getGenGammaLogPDFNF(Kappa(1), InvOmega(1), InvSigma(1))")
                    LogNormFac(1) = getGenGammaLogPDFNF(Kappa(1), InvOmega(1), InvSigma(1))
    call disp%show("LogNormFac(1)")
    call disp%show( LogNormFac(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at at a mix of scalar and vector input values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("LogNormFac(1:NP:NP/5) = getGenGammaLogPDFNF(Kappa(1:NP:NP/5))")
                    LogNormFac(1:NP:NP/5) = getGenGammaLogPDFNF(Kappa(1:NP:NP/5))
    call disp%show("LogNormFac(1:NP:NP/5)")
    call disp%show( LogNormFac(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("LogNormFac(1:NP:NP/5) = getGenGammaLogPDFNF(Kappa(1:NP:NP/5), 1._RK)")
                    LogNormFac(1:NP:NP/5) = getGenGammaLogPDFNF(Kappa(1:NP:NP/5), 1._RK)
    call disp%show("LogNormFac(1:NP:NP/5)")
    call disp%show( LogNormFac(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("InvOmega(1:NP:NP/5)")
    call disp%show( InvOmega(1:NP:NP/5) )
    call disp%show("LogNormFac(1:NP:NP/5) = getGenGammaLogPDFNF(1._RK, InvOmega(1:NP:NP/5))")
                    LogNormFac(1:NP:NP/5) = getGenGammaLogPDFNF(1._RK, InvOmega(1:NP:NP/5))
    call disp%show("LogNormFac(1:NP:NP/5)")
    call disp%show( LogNormFac(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("InvOmega(1:NP:NP/5)")
    call disp%show( InvOmega(1:NP:NP/5) )
    call disp%show("LogNormFac(1:NP:NP/5) = getGenGammaLogPDFNF(Kappa(1:NP:NP/5), InvOmega(1:NP:NP/5))")
                    LogNormFac(1:NP:NP/5) = getGenGammaLogPDFNF(Kappa(1:NP:NP/5), InvOmega(1:NP:NP/5))
    call disp%show("LogNormFac(1:NP:NP/5)")
    call disp%show( LogNormFac(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("InvOmega(1:NP:NP/5)")
    call disp%show( InvOmega(1:NP:NP/5) )
    call disp%show("InvSigma(1:NP:NP/5)")
    call disp%show( InvSigma(1:NP:NP/5) )
    call disp%show("LogNormFac(1:NP:NP/5) = getGenGammaLogPDFNF(Kappa(1:NP:NP/5), InvOmega(1:NP:NP/5), InvSigma(1:NP:NP/5))")
                    LogNormFac(1:NP:NP/5) = getGenGammaLogPDFNF(Kappa(1:NP:NP/5), InvOmega(1:NP:NP/5), InvSigma(1:NP:NP/5))
    call disp%show("LogNormFac(1:NP:NP/5)")
    call disp%show( LogNormFac(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example PDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        LogNormFac = getGenGammaLogPDFNF(Kappa)
        open(newunit = fileUnit, file = "getGenGammaLogPDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (Kappa(i), exp(LogNormFac(i)), i = 1, size(LogNormFac))
        close(fileUnit)
    end block

end program example