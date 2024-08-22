program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK
    use pm_kind, only: RK => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distGenGamma, only: setGenGammaCDF

    implicit none

    integer(IK) , parameter     :: NP = 1000_IK
    real(RK)    , dimension(NP) :: point, cdf, kappa, invOmega, invSigma
    integer(IK)                 :: info(NP)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    kappa = getLinSpace(+0.5_RK, +2._RK, count = NP)
    invOmega = getLogSpace(logx1 = log(0.1_RK), logx2 = log(10._RK), count = NP)
    invSigma = getLogSpace(-3._RK, +3._RK, count = NP)
    point = getLogSpace(log(.01_RK), log(+10._RK), count = NP)

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
    call disp%show("call setGenGammaCDF(cdf(1), 0.5_RK, info(1))")
                    call setGenGammaCDF(cdf(1), 0.5_RK, info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the cdf info.'")
                    if (info(1) < 0) error stop 'The computation of the cdf info.'
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("call setGenGammaCDF(cdf(1), 0.5_RK, kappa(1), info(1))")
                    call setGenGammaCDF(cdf(1), 0.5_RK, kappa(1), info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the cdf info.'")
                    if (info(1) < 0) error stop 'The computation of the cdf info.'
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("call setGenGammaCDF(cdf(1), 0.5_RK, kappa(1), invOmega(1), info(1))")
                    call setGenGammaCDF(cdf(1), 0.5_RK, kappa(1), invOmega(1), info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the cdf info.'")
                    if (info(1) < 0) error stop 'The computation of the cdf info.'
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("invOmega(1)")
    call disp%show( invOmega(1) )
    call disp%show("call setGenGammaCDF(cdf(1), 0.5_RK, kappa(1), invOmega(1), invSigma(1), info(1))")
                    call setGenGammaCDF(cdf(1), 0.5_RK, kappa(1), invOmega(1), invSigma(1), info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the cdf info.'")
                    if (info(1) < 0) error stop 'The computation of the cdf info.'
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
    call disp%show("call setGenGammaCDF(cdf(1:NP:NP/5), 0.5_RK, kappa(1:NP:NP/5), info(1:NP:NP/5))")
                    call setGenGammaCDF(cdf(1:NP:NP/5), 0.5_RK, kappa(1:NP:NP/5), info(1:NP:NP/5))
    call disp%show("if (info(1) < 0) error stop 'The computation of the cdf info.'")
                    if (info(1) < 0) error stop 'The computation of the cdf info.'
    call disp%show("cdf(1:NP:NP/5)")
    call disp%show( cdf(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("call setGenGammaCDF(cdf(1:NP:NP/5), point(1:NP:NP/5), kappa(1:NP:NP/5), info(1:NP:NP/5))")
                    call setGenGammaCDF(cdf(1:NP:NP/5), point(1:NP:NP/5), kappa(1:NP:NP/5), info(1:NP:NP/5))
    call disp%show("if (any(info(1:NP:NP/5) < 0)) error stop 'The computation of the cdf info.'")
                    if (any(info(1:NP:NP/5) < 0)) error stop 'The computation of the cdf info.'
    call disp%show("cdf(1:NP:NP/5)")
    call disp%show( cdf(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer     :: fileUnit, i
        real(RK)    :: cdf(NP,6)
        call setGenGammaCDF(cdf(:,1), point, 1.0, 0.5, 0.5, info)
        call setGenGammaCDF(cdf(:,2), point, 2.0, 0.5, 1.0, info)
        call setGenGammaCDF(cdf(:,3), point, 0.5, 2.0, 0.5, info)
        call setGenGammaCDF(cdf(:,4), point, 0.2, 5.0, 0.2, info)
        call setGenGammaCDF(cdf(:,5), point, .14, 7.0, .14, info)
        call setGenGammaCDF(cdf(:,6), point, 2.0, 5.0, 0.3, info)
        open(newunit = fileUnit, file = "setGenGammaCDF.RK.txt")
        if (any(info < 0)) error stop 'The computation of the cdf info.'
        write(fileUnit,"(7(g0,:,' '))") (point(i), cdf(i,:), i = 1, size(point))
        close(fileUnit)
    end block

end program example