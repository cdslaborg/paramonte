program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: RK => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distGenGamma, only: getGenGammaLogPDF

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RK)    , allocatable :: Point(:), logPDF(:), Kappa(:), invOmega(:), invSigma(:)

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
    call disp%show("logPDF(1) = getGenGammaLogPDF(0.5_RK)")
                    logPDF(1) = getGenGammaLogPDF(0.5_RK)
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("logPDF(1) = getGenGammaLogPDF(0.5_RK, Kappa(1))")
                    logPDF(1) = getGenGammaLogPDF(0.5_RK, Kappa(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("logPDF(1) = getGenGammaLogPDF(0.5_RK, Kappa(1), invOmega(1))")
                    logPDF(1) = getGenGammaLogPDF(0.5_RK, Kappa(1), invOmega(1))
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
    call disp%show("logPDF(1) = getGenGammaLogPDF(0.5_RK, Kappa(1), invOmega(1), invSigma(1))")
                    logPDF(1) = getGenGammaLogPDF(0.5_RK, Kappa(1), invOmega(1), invSigma(1))
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
    call disp%show("logPDF(1:NP:NP/5) = getGenGammaLogPDF(0.5_RK, Kappa(1:NP:NP/5))")
                    logPDF(1:NP:NP/5) = getGenGammaLogPDF(0.5_RK, Kappa(1:NP:NP/5))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("logPDF(1:NP:NP/5) = getGenGammaLogPDF(Point(1:NP:NP/5), Kappa(1:NP:NP/5))")
                    logPDF(1:NP:NP/5) = getGenGammaLogPDF(Point(1:NP:NP/5), Kappa(1:NP:NP/5))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer     :: fileUnit, i
        real(RK)    :: logPDF(NP, 6)
        logPDF(:,1) = getGenGammaLogPDF(Point, 1.0, 0.5, 0.5)
        logPDF(:,2) = getGenGammaLogPDF(Point, 2.0, 0.5, 1.0)
        logPDF(:,3) = getGenGammaLogPDF(Point, 0.5, 2.0, 0.5)
        logPDF(:,4) = getGenGammaLogPDF(Point, 0.2, 5.0, 0.2)
        logPDF(:,5) = getGenGammaLogPDF(Point, .14, 7.0, .14)
        logPDF(:,6) = getGenGammaLogPDF(Point, 2.0, 5.0, 0.3)
        open(newunit = fileUnit, file = "getGenGammaLogPDF.RK.txt")
        write(fileUnit,"(7(g0,:,' '))") (Point(i), exp(logPDF(i,:)), i = 1, size(Point))
        close(fileUnit)
    end block

end program example