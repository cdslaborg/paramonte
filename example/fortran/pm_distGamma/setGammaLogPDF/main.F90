program example

    use pm_kind, only: SK, IK, LK
    use pm_distGamma, only: getGammaLogPDFNF
    use pm_distGamma, only: setGammaLogPDF
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
    call disp%show("call setGammaLogPDF(logPDF(1), Point(1))")
                    call setGammaLogPDF(logPDF(1), Point(1))
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
    call disp%show("call setGammaLogPDF(logPDF(1:NP:NP/4), Point(1:NP:NP/4), getGammaLogPDFNF(2., 2.), kappa = 2., invSigma = 2.)")
                    call setGammaLogPDF(logPDF(1:NP:NP/4), Point(1:NP:NP/4), getGammaLogPDFNF(2., 2.), kappa = 2., invSigma = 2.)
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
    call disp%show("call setGammaLogPDF(logPDF(1:NP:NP/4), Point(1:NP:NP/4), getGammaLogPDFNF(0.5, 5.), kappa = 0.5, invSigma = 5.)")
                    call setGammaLogPDF(logPDF(1:NP:NP/4), Point(1:NP:NP/4), getGammaLogPDFNF(0.5, 5.), kappa = 0.5, invSigma = 5.)
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
    call disp%show("call setGammaLogPDF(logPDF(1:NP:NP/4), Point(NP/4), getGammaLogPDFNF(getLinSpace(0.5, 5., 5), getLinSpace(5., .5, 5)), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))")
                    call setGammaLogPDF(logPDF(1:NP:NP/4), Point(NP/4), getGammaLogPDFNF(getLinSpace(0.5, 5., 5), getLinSpace(5., .5, 5)), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))
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
    call disp%show("call setGammaLogPDF(logPDF(1:NP:NP/4), Point(1:NP:NP/4), getGammaLogPDFNF(getLinSpace(0.5, 5., 5), getLinSpace(5., .5, 5)), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))")
                    call setGammaLogPDF(logPDF(1:NP:NP/4), Point(1:NP:NP/4), getGammaLogPDFNF(getLinSpace(0.5, 5., 5), getLinSpace(5., .5, 5)), kappa = getLinSpace(0.5, 5., 5), invSigma = getLinSpace(5., .5, 5))
    call disp%show("logPDF(1:NP:NP/4)")
    call disp%show( logPDF(1:NP:NP/4) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example logPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        real, parameter :: Kappa(*) = [0.5, 1.0, 2.0, 7.5]
        real, parameter :: invSigma(*) = [1.0, 0.5, 0.5, 1.0]
        open(newunit = fileUnit, file = "setGammaLogPDF.RK.txt")
        do i = 1, NP
            call setGammaLogPDF(logPDF(1:4), Point(i), getGammaLogPDFNF(Kappa, invSigma), Kappa, invSigma)
            write(fileUnit,"(5(g0,:,' '))") Point(i), exp(logPDF(1:4))
        end do
        close(fileUnit)
    end block

end program example