program example

    use pm_kind, only: SK, IK, LK
    use pm_distGamma, only: getGammaLogPDF
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: getLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999_IK
    real                    :: Point(NP), logPDF(NP)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(Point, x1 = 0., x2 = 20., fopen = .true._LK, lopen = .true._LK)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Gamma distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(1)")
    call disp%show( Point(1) )
    call disp%show("logPDF(1) = getGammaLogPDF(Point(1), 2., 2.)")
                    logPDF(1) = getGammaLogPDF(Point(1), 2., 2.)
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Accelerate the runtime performance for repeated calls when `kappa` and `invSigma` are fixed (i.e., the PDF normalization constant is fixed).")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(1:NP:NP/4)")
    call disp%show( Point(1:NP:NP/4) )
    call disp%show("logPDF(1:NP:NP/4) = getGammaLogPDF(Point(1:NP:NP/4), 2., 2.)")
                    logPDF(1:NP:NP/4) = getGammaLogPDF(Point(1:NP:NP/4), 2., 2.)
    call disp%show("logPDF(1:NP:NP/4)")
    call disp%show( logPDF(1:NP:NP/4) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of PDF at different points with the same PDF parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(1:NP:NP/4)")
    call disp%show( Point(1:NP:NP/4) )
    call disp%show("logPDF(1:NP:NP/4) = getGammaLogPDF(Point(1:NP:NP/4), kappa = 0.5, invSigma = 5.)")
                    logPDF(1:NP:NP/4) = getGammaLogPDF(Point(1:NP:NP/4), kappa = 0.5, invSigma = 5.)
    call disp%show("logPDF(1:NP:NP/4)")
    call disp%show( logPDF(1:NP:NP/4) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of PDF at the same point but with different PDF parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(NP/4)")
    call disp%show( Point(NP/4) )
    call disp%show("logPDF(1:NP:NP/4) = getGammaLogPDF(Point(NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))")
                    logPDF(1:NP:NP/4) = getGammaLogPDF(Point(NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))
    call disp%show("logPDF(1:NP:NP/4)")
    call disp%show( logPDF(1:NP:NP/4) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of PDF at different points with different PDF parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(1:NP:NP/4)")
    call disp%show( Point(1:NP:NP/4) )
    call disp%show("logPDF(1:NP:NP/4) = getGammaLogPDF(Point(1:NP:NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))")
                    logPDF(1:NP:NP/4) = getGammaLogPDF(Point(1:NP:NP/4), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))
    call disp%show("logPDF(1:NP:NP/4)")
    call disp%show( logPDF(1:NP:NP/4) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example logPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "getGammaLogPDF.RK.txt")
        write(fileUnit,"(5(g0,:,' '))") ( Point(i) &
                                        , exp(getGammaLogPDF(Point(i), kappa = 0.5, invSigma = 1.0)) &
                                        , exp(getGammaLogPDF(Point(i), kappa = 1.0, invSigma = 0.5)) &
                                        , exp(getGammaLogPDF(Point(i), kappa = 2.0, invSigma = 0.5)) &
                                        , exp(getGammaLogPDF(Point(i), kappa = 7.5, invSigma = 1.0)) &
                                        , i = 1, NP &
                                        )
        close(fileUnit)
    end block

end program example