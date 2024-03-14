program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK
    use pm_kind, only: RK => RK32 ! all other real kinds are also acceptable: RK, RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distGenGamma, only: setGenGammaCDF

    implicit none

    integer(IK) , parameter     :: NP = 1000_IK
    real(RK)    , dimension(NP) :: Point, CDF, Kappa, InvOmega, InvSigma
    integer(IK)                 :: info(NP)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    Kappa = getLinSpace(+0.5_RK, +2._RK, count = NP)
    InvOmega = getLogSpace(logx1 = log(0.1_RK), logx2 = log(10._RK), count = NP)
    InvSigma = getLogSpace(-3._RK, +3._RK, count = NP)
    Point = getLogSpace(log(.01_RK), log(+10._RK), count = NP)

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
    call disp%show("call setGenGammaCDF(CDF(1), 0.5_RK, info(1))")
                    call setGenGammaCDF(CDF(1), 0.5_RK, info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the CDF info.'")
                    if (info(1) < 0) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setGenGammaCDF(CDF(1), 0.5_RK, log_gamma(Kappa(1)), Kappa(1), info(1))")
                    call setGenGammaCDF(CDF(1), 0.5_RK, log_gamma(Kappa(1)), Kappa(1), info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the CDF info.'")
                    if (info(1) < 0) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("InvOmega(1)")
    call disp%show( InvOmega(1) )
    call disp%show("call setGenGammaCDF(CDF(1), 0.5_RK, log_gamma(Kappa(1)), Kappa(1), InvOmega(1), info(1))")
                    call setGenGammaCDF(CDF(1), 0.5_RK, log_gamma(Kappa(1)), Kappa(1), InvOmega(1), info(1))
    call disp%show("if (info(1) < 0) error stop 'The computation of the CDF info.'")
                    if (info(1) < 0) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("InvOmega(1)")
    call disp%show( InvOmega(1) )
    call disp%show("call setGenGammaCDF(CDF(1), 0.5_RK, log_gamma(Kappa(1)), Kappa(1), InvOmega(1), InvSigma(1), info(1))")
                    call setGenGammaCDF(CDF(1), 0.5_RK, log_gamma(Kappa(1)), Kappa(1), InvOmega(1), InvSigma(1), info(1))
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
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setGenGammaCDF(CDF(1:NP:NP/5), 0.5_RK, log_gamma(Kappa(1:NP:NP/5)), Kappa(1:NP:NP/5), info(1:NP:NP/5))")
                    call setGenGammaCDF(CDF(1:NP:NP/5), 0.5_RK, log_gamma(Kappa(1:NP:NP/5)), Kappa(1:NP:NP/5), info(1:NP:NP/5))
    call disp%show("if (info(1) < 0) error stop 'The computation of the CDF info.'")
                    if (info(1) < 0) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("call setGenGammaCDF(CDF(1:NP:NP/5), Point(1:NP:NP/5), log_gamma(Kappa(1:NP:NP/5)), Kappa(1:NP:NP/5), info(1:NP:NP/5))")
                    call setGenGammaCDF(CDF(1:NP:NP/5), Point(1:NP:NP/5), log_gamma(Kappa(1:NP:NP/5)), Kappa(1:NP:NP/5), info(1:NP:NP/5))
    call disp%show("if (any(info(1:NP:NP/5) < 0)) error stop 'The computation of the CDF info.'")
                    if (any(info(1:NP:NP/5) < 0)) error stop 'The computation of the CDF info.'
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example CDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer     :: fileUnit, i
        real(RK)    :: CDF(NP,6)
        call setGenGammaCDF(CDF(:,1), Point, log_gamma(1.0), 1.0, 0.5, 0.5, info)
        call setGenGammaCDF(CDF(:,2), Point, log_gamma(2.0), 2.0, 0.5, 1.0, info)
        call setGenGammaCDF(CDF(:,3), Point, log_gamma(0.5), 0.5, 2.0, 0.5, info)
        call setGenGammaCDF(CDF(:,4), Point, log_gamma(0.2), 0.2, 5.0, 0.2, info)
        call setGenGammaCDF(CDF(:,5), Point, log_gamma(.14), .14, 7.0, .14, info)
        call setGenGammaCDF(CDF(:,6), Point, log_gamma(2.0), 2.0, 5.0, 0.3, info)
        open(newunit = fileUnit, file = "setGenGammaCDF.RK.txt")
        if (any(info < 0)) error stop 'The computation of the CDF info.'
        write(fileUnit,"(7(g0,:,' '))") (Point(i), exp(CDF(i,:)), i = 1, size(Point))
        close(fileUnit)
    end block

end program example