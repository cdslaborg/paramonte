program example

    use pm_kind, only: SK, IK, LK
    use pm_distGamma, only: getGammaCDF
    use pm_arraySpace, only: setLogSpace
    use pm_arraySpace, only: getLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999_IK
    real                    :: Point(NP), CDF(NP)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call setLogSpace(Point, logx1 = -5., logx2 = 5.)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of the Gamma distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(1)")
    call disp%show( Point(1) )
    call disp%show("CDF(1) = getGammaCDF(Point(1), 2., 2.)")
                    CDF(1) = getGammaCDF(Point(1), 2., 2.)
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Accelerate the runtime performance for repeated calls when `kappa` and `invSigma` are fixed.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(1:NP:NP/4)")
    call disp%show( Point(1:NP:NP/4) )
    call disp%show("CDF(1:NP:NP/4) = getGammaCDF(Point(1:NP:NP/4), 2., 2.)")
                    CDF(1:NP:NP/4) = getGammaCDF(Point(1:NP:NP/4), 2., 2.)
    call disp%show("CDF(1:NP:NP/4)")
    call disp%show( CDF(1:NP:NP/4) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at different points with the same CDF parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(1:NP:NP/4)")
    call disp%show( Point(1:NP:NP/4) )
    call disp%show("CDF(1:NP:NP/4) = getGammaCDF(Point(1:NP:NP/4), kappa = 0.5, invSigma = 5.)")
                    CDF(1:NP:NP/4) = getGammaCDF(Point(1:NP:NP/4), kappa = 0.5, invSigma = 5.)
    call disp%show("CDF(1:NP:NP/4)")
    call disp%show( CDF(1:NP:NP/4) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at the same point but with different CDF parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(NP/4)")
    call disp%show( Point(NP/4) )
    call disp%show("CDF(1:NP:NP/4) = getGammaCDF(Point(NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))")
                    CDF(1:NP:NP/4) = getGammaCDF(Point(NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))
    call disp%show("CDF(1:NP:NP/4)")
    call disp%show( CDF(1:NP:NP/4) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at different points with different CDF parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(1:NP:NP/4)")
    call disp%show( Point(1:NP:NP/4) )
    call disp%show("CDF(1:NP:NP/4) = getGammaCDF(Point(1:NP:NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))")
                    CDF(1:NP:NP/4) = getGammaCDF(Point(1:NP:NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))
    call disp%show("CDF(1:NP:NP/4)")
    call disp%show( CDF(1:NP:NP/4) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example CDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "getGammaCDF.RK.txt")
        write(fileUnit,"(5(g0,:,' '))") (Point(i), getGammaCDF(Point(i), kappa = [0.5, 1.0, 2.0, 7.5], invSigma = [1.0, 0.5, 0.5, 1.0]), i = 1, NP)
        close(fileUnit)
    end block

end program example