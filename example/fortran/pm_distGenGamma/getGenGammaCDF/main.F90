program example

    use pm_kind, only: SK
    use pm_kind, only: IK
    use pm_kind, only: RK => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distGenGamma, only: getGenGammaCDF

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RK)    , allocatable :: point(:), cdf(:), kappa(:), invOmega(:), invSigma(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    kappa = getLinSpace(+0.5_RK, +2._RK, count = NP)
    invOmega = getLogSpace(logx1 = log(0.1_RK), logx2 = log(10._RK), count = NP)
    invSigma = getLogSpace(-3._RK, +3._RK, count = NP)
    point = getLogSpace(log(.01_RK), log(+10._RK), count = NP)
    allocate(cdf, mold = point)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of GenGamma distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the cdf at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("cdf(1) = getGenGammaCDF(0.5_RK)")
                    cdf(1) = getGenGammaCDF(0.5_RK)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("cdf(1) = getGenGammaCDF(0.5_RK, kappa(1))")
                    cdf(1) = getGenGammaCDF(0.5_RK, kappa(1))
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("cdf(1) = getGenGammaCDF(0.5_RK, kappa(1), invOmega(1))")
                    cdf(1) = getGenGammaCDF(0.5_RK, kappa(1), invOmega(1))
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("cdf(1) = getGenGammaCDF(0.5_RK, kappa(1), invOmega(1), invSigma(1))")
                    cdf(1) = getGenGammaCDF(0.5_RK, kappa(1), invOmega(1), invSigma(1))
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the cdf at an input vector real value with different parameter values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("cdf(1:NP:NP/5) = getGenGammaCDF(0.5_RK, kappa(1:NP:NP/5))")
                    cdf(1:NP:NP/5) = getGenGammaCDF(0.5_RK, kappa(1:NP:NP/5))
    call disp%show("cdf(1:NP:NP/5)")
    call disp%show( cdf(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("cdf(1:NP:NP/5) = getGenGammaCDF(point(1:NP:NP/5), kappa(1:NP:NP/5))")
                    cdf(1:NP:NP/5) = getGenGammaCDF(point(1:NP:NP/5), kappa(1:NP:NP/5))
    call disp%show("cdf(1:NP:NP/5)")
    call disp%show( cdf(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer     :: fileUnit, i
        real(RK)    :: cdf(NP,6)
        cdf(:,1) = getGenGammaCDF(point, 1.0, 0.5, 0.5)
        cdf(:,2) = getGenGammaCDF(point, 2.0, 0.5, 1.0)
        cdf(:,3) = getGenGammaCDF(point, 0.5, 2.0, 0.5)
        cdf(:,4) = getGenGammaCDF(point, 0.2, 5.0, 0.2)
        cdf(:,5) = getGenGammaCDF(point, .14, 7.0, .14)
        cdf(:,6) = getGenGammaCDF(point, 2.0, 5.0, 0.3)
        open(newunit = fileUnit, file = "getGenGammaCDF.RK.txt")
        write(fileUnit,"(7(g0,:,' '))") (point(i), cdf(i,:), i = 1, size(point))
        close(fileUnit)
    end block

end program example