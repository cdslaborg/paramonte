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
    real(RK)    , allocatable :: point(:), cdf(:), kappa(:), invOmega(:), logSigma(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    kappa = getLinSpace(+0.5_RK, +2._RK, count = NP)
    invOmega = getLogSpace(logx1 = log(0.1_RK), logx2 = log(10._RK), count = NP)
    logSigma = getLinSpace(-3._RK, +3._RK, count = NP)
    point = getLinSpace(-10._RK, +10._RK, count = NP)
    allocate(cdf, mold = point)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (cdf) of GenExpGamma distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the cdf at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("cdf(1) = getGenExpGammaCDF(0.5_RK)")
                    cdf(1) = getGenExpGammaCDF(0.5_RK)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("cdf(1) = getGenExpGammaCDF(0.5_RK, kappa(1))")
                    cdf(1) = getGenExpGammaCDF(0.5_RK, kappa(1))
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("cdf(1) = getGenExpGammaCDF(0.5_RK, kappa(1), invOmega(1))")
                    cdf(1) = getGenExpGammaCDF(0.5_RK, kappa(1), invOmega(1))
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("cdf(1) = getGenExpGammaCDF(0.5_RK, kappa(1), invOmega(1), logSigma(1))")
                    cdf(1) = getGenExpGammaCDF(0.5_RK, kappa(1), invOmega(1), logSigma(1))
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
    call disp%show("cdf(1:NP:NP/5) = getGenExpGammaCDF(0.5_RK, kappa(1:NP:NP/5))")
                    cdf(1:NP:NP/5) = getGenExpGammaCDF(0.5_RK, kappa(1:NP:NP/5))
    call disp%show("cdf(1:NP:NP/5)")
    call disp%show( cdf(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("cdf(1:NP:NP/5) = getGenExpGammaCDF(point(1:NP:NP/5), kappa(1:NP:NP/5))")
                    cdf(1:NP:NP/5) = getGenExpGammaCDF(point(1:NP:NP/5), kappa(1:NP:NP/5))
    call disp%show("cdf(1:NP:NP/5)")
    call disp%show( cdf(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        real(RK) :: cdf(NP, 6)
        integer(IK) :: i, fileUnit
        cdf(:,1) = getGenExpGammaCDF(point, .5_RK, 2._RK)
        cdf(:,2) = getGenExpGammaCDF(point, 1._RK, 2._RK)
        cdf(:,3) = getGenExpGammaCDF(point, 2._RK, 2._RK)
        cdf(:,4) = getGenExpGammaCDF(point, .5_RK, .5_RK)
        cdf(:,5) = getGenExpGammaCDF(point, 1._RK, .5_RK)
        cdf(:,6) = getGenExpGammaCDF(point, 2._RK, .5_RK)
        open(newunit = fileUnit, file = "getGenExpGammaCDF.RK.txt")
        write(fileUnit,"(7(g0,:,' '))") (point(i), cdf(i,:), i = 1, size(point))
        close(fileUnit)
    end block

end program example