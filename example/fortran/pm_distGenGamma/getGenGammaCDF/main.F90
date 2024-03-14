program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: RK => RK32 ! all other real kinds are also acceptable: RK, RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distGenGamma, only: getGenGammaCDF

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RK)    , allocatable :: Point(:), CDF(:), Kappa(:), InvOmega(:), InvSigma(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLinSpace(+0.5_RK, +2._RK, count = NP)
    InvOmega = getLogSpace(logx1 = log(0.1_RK), logx2 = log(10._RK), count = NP)
    InvSigma = getLogSpace(-3._RK, +3._RK, count = NP)
    Point = getLogSpace(log(.01_RK), log(+10._RK), count = NP)
    allocate(CDF, mold = Point)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of GenGamma distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the CDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("CDF(1) = getGenGammaCDF(0.5_RK)")
                    CDF(1) = getGenGammaCDF(0.5_RK)
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("CDF(1) = getGenGammaCDF(0.5_RK, Kappa(1))")
                    CDF(1) = getGenGammaCDF(0.5_RK, Kappa(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("InvOmega(1)")
    call disp%show( InvOmega(1) )
    call disp%show("CDF(1) = getGenGammaCDF(0.5_RK, Kappa(1), InvOmega(1))")
                    CDF(1) = getGenGammaCDF(0.5_RK, Kappa(1), InvOmega(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("InvOmega(1)")
    call disp%show( InvOmega(1) )
    call disp%show("CDF(1) = getGenGammaCDF(0.5_RK, Kappa(1), InvOmega(1), InvSigma(1))")
                    CDF(1) = getGenGammaCDF(0.5_RK, Kappa(1), InvOmega(1), InvSigma(1))
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
    call disp%show("CDF(1:NP:NP/5) = getGenGammaCDF(0.5_RK, Kappa(1:NP:NP/5))")
                    CDF(1:NP:NP/5) = getGenGammaCDF(0.5_RK, Kappa(1:NP:NP/5))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("CDF(1:NP:NP/5) = getGenGammaCDF(Point(1:NP:NP/5), Kappa(1:NP:NP/5))")
                    CDF(1:NP:NP/5) = getGenGammaCDF(Point(1:NP:NP/5), Kappa(1:NP:NP/5))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example CDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer     :: fileUnit, i
        real(RK)    :: CDF(NP,6)
        CDF(:,1) = getGenGammaCDF(Point, 1.0, 0.5, 0.5)
        CDF(:,2) = getGenGammaCDF(Point, 2.0, 0.5, 1.0)
        CDF(:,3) = getGenGammaCDF(Point, 0.5, 2.0, 0.5)
        CDF(:,4) = getGenGammaCDF(Point, 0.2, 5.0, 0.2)
        CDF(:,5) = getGenGammaCDF(Point, .14, 7.0, .14)
        CDF(:,6) = getGenGammaCDF(Point, 2.0, 5.0, 0.3)
        open(newunit = fileUnit, file = "getGenGammaCDF.RK.txt")
        write(fileUnit,"(7(g0,:,' '))") (Point(i), CDF(i,:), i = 1, size(Point))
        close(fileUnit)
    end block

end program example