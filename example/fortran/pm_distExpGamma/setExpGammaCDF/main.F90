program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK
    use pm_kind, only: RK => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distExpGamma, only: setExpGammaCDF

    implicit none

    integer(IK) , parameter     :: NP = 1000_IK
    real(RK)    , dimension(NP) :: point, CDF, kappa, logSigma
    integer(IK)                 :: info(NP)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    kappa = getLinSpace(+0.5_RK, +2._RK, count = NP)
    logSigma = getLinSpace(-3._RK, +3._RK, count = NP)
    point = getLinSpace(-10._RK, +10._RK, count = NP)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of ExpGamma distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the CDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setExpGammaCDF(CDF(1), 0.5_RK, info(1))")
                    call setExpGammaCDF(CDF(1), 0.5_RK, info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the CDF info.'")
                    if (info(1) < 0) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("call setExpGammaCDF(CDF(1), 0.5_RK, kappa(1), info(1))")
                    call setExpGammaCDF(CDF(1), 0.5_RK, kappa(1), info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the CDF info.'")
                    if (info(1) < 0) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("call setExpGammaCDF(CDF(1), 0.5_RK, kappa(1), info(1))")
                    call setExpGammaCDF(CDF(1), 0.5_RK, kappa(1), info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the CDF info.'")
                    if (info(1) < 0) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("call setExpGammaCDF(CDF(1), 0.5_RK, kappa(1), logSigma(1), info(1))")
                    call setExpGammaCDF(CDF(1), 0.5_RK, kappa(1), logSigma(1), info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the CDF info.'")
                    if (info(1) < 0) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the CDF at an input vector real value with different parameter values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("call setExpGammaCDF(CDF(1:NP:NP/5), 0.5_RK, kappa(1:NP:NP/5), info(1:NP:NP/5))")
                    call setExpGammaCDF(CDF(1:NP:NP/5), 0.5_RK, kappa(1:NP:NP/5), info(1:NP:NP/5))
    call disp%show("if (info(1) < 0) error stop 'The computation of the CDF info.'")
                    if (info(1) < 0) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1)")
    call disp%show( kappa(1) )
    call disp%show("call setExpGammaCDF(CDF(1:NP:NP/5), point(1:NP:NP/5), kappa(1:NP:NP/5), info(1:NP:NP/5))")
                    call setExpGammaCDF(CDF(1:NP:NP/5), point(1:NP:NP/5), kappa(1:NP:NP/5), info(1:NP:NP/5))
    call disp%show("if (any(info(1:NP:NP/5) < 0)) error stop 'The computation of the CDF info.'")
                    if (any(info(1:NP:NP/5) < 0)) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        real(RK) :: CDF(NP,6)
        integer :: fileUnit, i
        call setExpGammaCDF(CDF(:,1), point, .5_RK, 0.5, info)
        call setExpGammaCDF(CDF(:,2), point, 1._RK, 0.5, info)
        call setExpGammaCDF(CDF(:,3), point, 2._RK, 0.5, info)
        call setExpGammaCDF(CDF(:,4), point, .5_RK, 2.0, info)
        call setExpGammaCDF(CDF(:,5), point, 1._RK, 2.0, info)
        call setExpGammaCDF(CDF(:,6), point, 2._RK, 2.0, info)
        open(newunit = fileUnit, file = "setExpGammaCDF.RK.txt")
        if (any(info < 0)) error stop 'The computation of the CDF info.'
        write(fileUnit,"(7(g0,:,' '))") (point(i), CDF(i,:), i = 1, size(point))
        close(fileUnit)
    end block

end program example