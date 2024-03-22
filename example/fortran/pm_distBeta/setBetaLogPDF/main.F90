program example

    use pm_kind, only: SK, IK
    use pm_distBeta, only: setBetaLogPDF
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
    call disp%show("call setBetaLogPDF(logPDF(1), Point(1), 2., 2.)")
                    call setBetaLogPDF(logPDF(1), Point(1), 2., 2.)
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
    call disp%show("call setBetaLogPDF(logPDF(NP/2), Point(NP/2), 2., 2., getLogBeta(2., 2.))")
                    call setBetaLogPDF(logPDF(NP/2), Point(NP/2), 2., 2., getLogBeta(2., 2.))
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
    call disp%show("call setBetaLogPDF(logPDF(1:NP:NP/4), Point(1:NP:NP/4), alpha = 0.5, beta = 5.)")
                    call setBetaLogPDF(logPDF(1:NP:NP/4), Point(1:NP:NP/4), alpha = 0.5, beta = 5.)
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
    call disp%show("call setBetaLogPDF(logPDF(1:NP:NP/4), Point(NP/4), alpha = getLinSpace(0.5, 5., 5), beta = getLinSpace(5., .5, 5))")
                    call setBetaLogPDF(logPDF(1:NP:NP/4), Point(NP/4), alpha = getLinSpace(0.5, 5., 5), beta = getLinSpace(5., .5, 5))
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
    call disp%show("call setBetaLogPDF(logPDF(1:NP:NP/4), Point(1:NP:NP/4), alpha = getLinSpace(0.5, 5., 5), beta = getLinSpace(5., .5, 5))")
                    call setBetaLogPDF(logPDF(1:NP:NP/4), Point(1:NP:NP/4), alpha = getLinSpace(0.5, 5., 5), beta = getLinSpace(5., .5, 5))
    call disp%show("logPDF(1:NP:NP/4)")
    call disp%show( logPDF(1:NP:NP/4) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example logPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "setBetaLogPDF.RK.txt")
        do i = 1, NP
            call setBetaLogPDF(logPDF(1:4), Point(i), alpha = [0.5, 2.0, 2.0, 5.0], beta = [0.5, 2.0, 5.0, 2.0])
            write(fileUnit,"(5(g0,:,' '))") Point(i), exp(logPDF(1:4))
        end do
        close(fileUnit)
    end block

end program example