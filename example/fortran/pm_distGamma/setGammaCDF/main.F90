program example

    use pm_kind, only: SK, IK, LK
    use pm_distGamma, only: setGammaCDF
    use pm_arraySpace, only: setLogSpace
    use pm_arraySpace, only: getLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999_IK
    real                    :: point(NP), CDF(NP)
    integer(IK)             :: info(NP)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call setLogSpace(point, logx1 = -5., logx2 = 5.)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of the Gamma distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("point(1)")
    call disp%show( point(1) )
    call disp%show("call setGammaCDF(CDF(1), point(1), info(1))")
                    call setGammaCDF(CDF(1), point(1), info(1))
    call disp%show("if (info(1) < 0) error stop 'setGammaCDF() failed.'")
                    if (info(1) < 0) error stop 'setGammaCDF() failed.'
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Accelerate the runtime performance for repeated calls when `kappa` and `invSigma` are fixed.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("point(1:NP:NP/4)")
    call disp%show( point(1:NP:NP/4) )
    call disp%show("call setGammaCDF(CDF(1:NP:NP/4), point(1:NP:NP/4), kappa = 2., invSigma = 2., info = info(1:NP:NP/4))")
                    call setGammaCDF(CDF(1:NP:NP/4), point(1:NP:NP/4), kappa = 2., invSigma = 2., info = info(1:NP:NP/4))
    call disp%show("if (any(info(1:NP:NP/4) < 0)) error stop 'setGammaCDF() failed.'")
                    if (any(info(1:NP:NP/4) < 0)) error stop 'setGammaCDF() failed.'
    call disp%show("CDF(1:NP:NP/4)")
    call disp%show( CDF(1:NP:NP/4) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at different points with the same CDF parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("point(1:NP:NP/4)")
    call disp%show( point(1:NP:NP/4) )
    call disp%show("call setGammaCDF(CDF(1:NP:NP/4), point(1:NP:NP/4), kappa = 0.5, invSigma = 5., info = info(1:NP:NP/4))")
                    call setGammaCDF(CDF(1:NP:NP/4), point(1:NP:NP/4), kappa = 0.5, invSigma = 5., info = info(1:NP:NP/4))
    call disp%show("if (any(info(1:NP:NP/4) < 0)) error stop 'setGammaCDF() failed.'")
                    if (any(info(1:NP:NP/4) < 0)) error stop 'setGammaCDF() failed.'
    call disp%show("CDF(1:NP:NP/4)")
    call disp%show( CDF(1:NP:NP/4) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at the same point but with different CDF parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("point(NP/4)")
    call disp%show( point(NP/4) )
    call disp%show("call setGammaCDF(CDF(1:NP:NP/4), point(NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5), info = info(1:NP:NP/4))")
                    call setGammaCDF(CDF(1:NP:NP/4), point(NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5), info = info(1:NP:NP/4))
    call disp%show("if (any(info(1:NP:NP/4) < 0)) error stop 'setGammaCDF() failed.'")
                    if (any(info(1:NP:NP/4) < 0)) error stop 'setGammaCDF() failed.'
    call disp%show("CDF(1:NP:NP/4)")
    call disp%show( CDF(1:NP:NP/4) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at different points with different CDF parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("point(1:NP:NP/4)")
    call disp%show( point(1:NP:NP/4) )
    call disp%show("call setGammaCDF(CDF(1:NP:NP/4), point(1:NP:NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5), info = info(1:NP:NP/4))")
                    call setGammaCDF(CDF(1:NP:NP/4), point(1:NP:NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5), info = info(1:NP:NP/4))
    call disp%show("if (any(info(1:NP:NP/4) < 0)) error stop 'setGammaCDF() failed.'")
                    if (any(info(1:NP:NP/4) < 0)) error stop 'setGammaCDF() failed.'
    call disp%show("CDF(1:NP:NP/4)")
    call disp%show( CDF(1:NP:NP/4) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        real, parameter :: kappa(*) = [0.5, 1.0, 2.0, 7.5]
        real, parameter :: invSigma(*) = [1.0, 0.5, 0.5, 1.0]
        open(newunit = fileUnit, file = "setGammaCDF.RK.txt")
        do i = 1, NP
            call setGammaCDF(CDF(1:4), point(i), kappa, invSigma, info(1:4))
            if (any(info(1:4) < 0)) error stop 'setGammaCDF() failed.'
            write(fileUnit,"(5(g0,:,' '))") point(i), CDF(1:4)
        end do
        close(fileUnit)
    end block

end program example