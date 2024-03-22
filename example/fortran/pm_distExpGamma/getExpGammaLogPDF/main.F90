program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: RK => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distExpGamma, only: getExpGammaLogPDF

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RK)    , allocatable :: Point(:), logPDF(:), Kappa(:), LogSigma(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLinSpace(+0.5_RK, +2._RK, count = NP)
    LogSigma = getLinSpace(-3._RK, +3._RK, count = NP)
    Point = getLinSpace(-10._RK, +10._RK, count = NP)
    allocate(logPDF, mold = Point)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of ExpGamma distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the PDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDF(1) = getExpGammaLogPDF(0.5_RK)")
                    logPDF(1) = getExpGammaLogPDF(0.5_RK)
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("logPDF(1) = getExpGammaLogPDF(0.5_RK, Kappa(1))")
                    logPDF(1) = getExpGammaLogPDF(0.5_RK, Kappa(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("logPDF(1) = getExpGammaLogPDF(0.5_RK, Kappa(1), LogSigma(1))")
                    logPDF(1) = getExpGammaLogPDF(0.5_RK, Kappa(1), LogSigma(1))
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
    call disp%show("logPDF(1:NP:NP/5) = getExpGammaLogPDF(0.5_RK, Kappa(1:NP:NP/5))")
                    logPDF(1:NP:NP/5) = getExpGammaLogPDF(0.5_RK, Kappa(1:NP:NP/5))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("logPDF(1:NP:NP/5) = getExpGammaLogPDF(Point(1:NP:NP/5), Kappa(1:NP:NP/5))")
                    logPDF(1:NP:NP/5) = getExpGammaLogPDF(Point(1:NP:NP/5), Kappa(1:NP:NP/5))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer     :: fileUnit, i
        real        :: logPDF(NP,6)
        logPDF(:,1) = getExpGammaLogPDF(Point, .5, -.5)
        logPDF(:,2) = getExpGammaLogPDF(Point, 1., -.5)
        logPDF(:,3) = getExpGammaLogPDF(Point, 2., -.5)
        logPDF(:,4) = getExpGammaLogPDF(Point, .5, 2.0)
        logPDF(:,5) = getExpGammaLogPDF(Point, 1., 2.0)
        logPDF(:,6) = getExpGammaLogPDF(Point, 2., 2.0)
        open(newunit = fileUnit, file = "getExpGammaLogPDF.RK.txt")
        write(fileUnit,"(7(g0,:,' '))") (Point(i), exp(logPDF(i,:)), i = 1, size(Point))
        close(fileUnit)
    end block

end program example