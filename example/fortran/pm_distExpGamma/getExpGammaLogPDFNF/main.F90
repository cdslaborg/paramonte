program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_kind, only: RK => RK32 ! all other real kinds are also supported, e.g.: RK, RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLinSpace
    use pm_distExpGamma, only: getExpGammaLogPDFNF

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RK)    , allocatable :: LogNormFac(:), Kappa(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLinSpace(0.01_RK, 10._RK, count = NP)
    allocate(LogNormFac, mold = Kappa)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the ExpGamma distribution PDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("LogNormFac(1:NP:NP/5) = getExpGammaLogPDFNF(1._RK)")
                    LogNormFac(1:NP:NP/5) = getExpGammaLogPDFNF(1._RK)
    call disp%show("LogNormFac(1:NP:NP/5)")
    call disp%show( LogNormFac(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("LogNormFac(1) = getExpGammaLogPDFNF(kappa = Kappa(1))")
                    LogNormFac(1) = getExpGammaLogPDFNF(kappa = Kappa(1))
    call disp%show("LogNormFac(1)")
    call disp%show( LogNormFac(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/5)")
    call disp%show( Kappa(1:NP:NP/5) )
    call disp%show("LogNormFac(1:NP:NP/5) = getExpGammaLogPDFNF(Kappa(1:NP:NP/5))")
                    LogNormFac(1:NP:NP/5) = getExpGammaLogPDFNF(Kappa(1:NP:NP/5))
    call disp%show("LogNormFac(1:NP:NP/5)")
    call disp%show( LogNormFac(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example PDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    LogNormFac = getExpGammaLogPDFNF(Kappa)

    block

        integer :: fileUnit, i

        open(newunit = fileUnit, file = "getExpGammaLogPDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (Kappa(i), exp(LogNormFac(i)), i = 1, size(LogNormFac))
        close(fileUnit)

    end block

end program example