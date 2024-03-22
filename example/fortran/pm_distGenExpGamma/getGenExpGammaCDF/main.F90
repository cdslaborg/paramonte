program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: RK => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distGenExpGamma, only: getGenExpGammaCDF

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RK)    , allocatable :: Point(:), CDF(:), Kappa(:), invOmega(:), LogSigma(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLinSpace(+0.5_RK, +2._RK, count = NP)
    invOmega = getLogSpace(logx1 = log(0.1_RK), logx2 = log(10._RK), count = NP)
    LogSigma = getLinSpace(-3._RK, +3._RK, count = NP)
    Point = getLinSpace(-10._RK, +10._RK, count = NP)
    allocate(CDF, mold = Point)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of GenExpGamma distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the CDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("CDF(1) = getGenExpGammaCDF(0.5_RK)")
                    CDF(1) = getGenExpGammaCDF(0.5_RK)
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("CDF(1) = getGenExpGammaCDF(0.5_RK, Kappa(1))")
                    CDF(1) = getGenExpGammaCDF(0.5_RK, Kappa(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("CDF(1) = getGenExpGammaCDF(0.5_RK, Kappa(1), invOmega(1))")
                    CDF(1) = getGenExpGammaCDF(0.5_RK, Kappa(1), invOmega(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("CDF(1) = getGenExpGammaCDF(0.5_RK, Kappa(1), invOmega(1), LogSigma(1))")
                    CDF(1) = getGenExpGammaCDF(0.5_RK, Kappa(1), invOmega(1), LogSigma(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the CDF at an input vector real value with different parameter values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("CDF(1:NP:NP/5) = getGenExpGammaCDF(0.5_RK, Kappa(1:NP:NP/5))")
                    CDF(1:NP:NP/5) = getGenExpGammaCDF(0.5_RK, Kappa(1:NP:NP/5))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("CDF(1:NP:NP/5) = getGenExpGammaCDF(Point(1:NP:NP/5), Kappa(1:NP:NP/5))")
                    CDF(1:NP:NP/5) = getGenExpGammaCDF(Point(1:NP:NP/5), Kappa(1:NP:NP/5))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer     :: fileUnit, i
        real(RK)    :: CDF(NP,6)
        CDF(:,1) = getGenExpGammaCDF(Point, .5_RK, 2._RK)
        CDF(:,2) = getGenExpGammaCDF(Point, 1._RK, 2._RK)
        CDF(:,3) = getGenExpGammaCDF(Point, 2._RK, 2._RK)
        CDF(:,4) = getGenExpGammaCDF(Point, .5_RK, .5_RK)
        CDF(:,5) = getGenExpGammaCDF(Point, 1._RK, .5_RK)
        CDF(:,6) = getGenExpGammaCDF(Point, 2._RK, .5_RK)
        open(newunit = fileUnit, file = "getGenExpGammaCDF.RK.txt")
        write(fileUnit,"(7(g0,:,' '))") (Point(i), CDF(i,:), i = 1, size(Point))
        close(fileUnit)
    end block

end program example