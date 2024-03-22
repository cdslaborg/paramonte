program example

    use pm_kind, only: SK, IK
    use pm_distBeta, only: getBetaLogPDF
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: getLinSpace
    use pm_mathBeta, only: getLogBeta
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999_IK
    real                    :: Point(NP), logPDF(NP)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(Point, x1 = 0.001, x2 = .999)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the Beta distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(1)")
    call disp%show( Point(1) )
    call disp%show("logPDF(1) = getBetaLogPDF(Point(1), 2., 2.)")
                    logPDF(1) = getBetaLogPDF(Point(1), 2., 2.)
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Accelerate the runtime performance for repeated calls when `alpha` and `beta` are fixed (i.e., the PDF normalization constant is fixed).")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(NP/2)")
    call disp%show( Point(NP/2) )
    call disp%show("logPDF(NP/2) = getBetaLogPDF(Point(NP/2), 2., 2., getLogBeta(2., 2.))")
                    logPDF(NP/2) = getBetaLogPDF(Point(NP/2), 2., 2., getLogBeta(2., 2.))
    call disp%show("logPDF(NP/2)")
    call disp%show( logPDF(NP/2) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of PDF at different points with the same PDF parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(1:NP:NP/4)")
    call disp%show( Point(1:NP:NP/4) )
    call disp%show("logPDF(1:NP:NP/4) = getBetaLogPDF(Point(1:NP:NP/4), alpha = 0.5, beta = 5.)")
                    logPDF(1:NP:NP/4) = getBetaLogPDF(Point(1:NP:NP/4), alpha = 0.5, beta = 5.)
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
    call disp%show("logPDF(1:NP:NP/4) = getBetaLogPDF(Point(NP/4), alpha = getLinSpace(0.5, 5., 5), beta = getLinSpace(5., .5, 5))")
                    logPDF(1:NP:NP/4) = getBetaLogPDF(Point(NP/4), alpha = getLinSpace(0.5, 5., 5), beta = getLinSpace(5., .5, 5))
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
    call disp%show("logPDF(1:NP:NP/4) = getBetaLogPDF(Point(1:NP:NP/4), alpha = getLinSpace(0.5, 5., 5), beta = getLinSpace(5., .5, 5))")
                    logPDF(1:NP:NP/4) = getBetaLogPDF(Point(1:NP:NP/4), alpha = getLinSpace(0.5, 5., 5), beta = getLinSpace(5., .5, 5))
    call disp%show("logPDF(1:NP:NP/4)")
    call disp%show( logPDF(1:NP:NP/4) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example logPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "getBetaLogPDF.RK.txt")
        write(fileUnit,"(5(g0,:,' '))") ( Point(i) &
                                        , exp(getBetaLogPDF(Point(i), alpha = 0.5, beta = 0.5)) &
                                        , exp(getBetaLogPDF(Point(i), alpha = 2.0, beta = 2.0)) &
                                        , exp(getBetaLogPDF(Point(i), alpha = 2.0, beta = 5.0)) &
                                        , exp(getBetaLogPDF(Point(i), alpha = 5.0, beta = 2.0)) &
                                        , i = 1, NP &
                                        )
        close(fileUnit)
    end block

end program example