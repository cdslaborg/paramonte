program example

    use pm_kind, only: SK
    use pm_kind, only: IK, RK ! all real kinds are supported.
    use pm_distLogNorm, only: getLogNormLogPDF
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    real(RK), dimension(NP) :: point, mu, Sigma, logPDF

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(mu, x1 = -5._RK, x2 = +5._RK)
    call setLogSpace(point, logx1 = log(0.0001_RK), logx2 = log(5._RK))
    call setLogSpace(Sigma, logx1 = log(0.1_RK), logx2 = log(10._RK))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the (Standard) LogNormal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("! Standard PDF.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("point(NP/2)")
    call disp%show( point(NP/2) )
    call disp%show("logPDF(1) = getLogNormLogPDF(point(NP/2))")
                    logPDF(1) = getLogNormLogPDF(point(NP/2))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("! PDF with a mean.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("point(1)")
    call disp%show( point(1) )
    call disp%show("logPDF(1) = getLogNormLogPDF(point(1), mu(1))")
                    logPDF(1) = getLogNormLogPDF(point(1), mu(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! PDF with a standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Sigma(1)")
    call disp%show( Sigma(1) )
    call disp%show("point(1)")
    call disp%show( point(1) )
    call disp%show("logPDF(1) = getLogNormLogPDF(point(1), sigma = Sigma(1))")
                    logPDF(1) = getLogNormLogPDF(point(1), sigma = Sigma(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! PDF with a mean and a standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("Sigma(1)")
    call disp%show( Sigma(1) )
    call disp%show("point(1)")
    call disp%show( point(1) )
    call disp%show("logPDF(1) = getLogNormLogPDF(point(1), mu(1), Sigma(1))")
                    logPDF(1) = getLogNormLogPDF(point(1), mu(1), Sigma(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of PDF at different points with the same mean and standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("Sigma(1)")
    call disp%show( Sigma(1) )
    call disp%show("point(1:NP:NP/5)")
    call disp%show( point(1:NP:NP/5) )
    call disp%show("logPDF(1:NP:NP/5) = getLogNormLogPDF(point(1:NP:NP/5), mu(1), Sigma(1))")
                    logPDF(1:NP:NP/5) = getLogNormLogPDF(point(1:NP:NP/5), mu(1), Sigma(1))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of PDF at the same point but with different means and standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1:NP:NP/5)")
    call disp%show( mu(1:NP:NP/5) )
    call disp%show("Sigma(1:NP:NP/5)")
    call disp%show( Sigma(1:NP:NP/5) )
    call disp%show("point(1)")
    call disp%show( point(1) )
    call disp%show("logPDF(1:NP:NP/5) = getLogNormLogPDF(point(1), mu(1:NP:NP/5), Sigma(1:NP:NP/5))")
                    logPDF(1:NP:NP/5) = getLogNormLogPDF(point(1), mu(1:NP:NP/5), Sigma(1:NP:NP/5))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of PDF at different points with different means and a standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1:NP:NP/5)")
    call disp%show( mu(1:NP:NP/5) )
    call disp%show("Sigma(1:NP:NP/5)")
    call disp%show( Sigma(1:NP:NP/5) )
    call disp%show("point(1:NP:NP/5)")
    call disp%show( point(1:NP:NP/5) )
    call disp%show("logPDF(1:NP:NP/5) = getLogNormLogPDF(point(1:NP:NP/5), mu(1:NP:NP/5), Sigma(1:NP:NP/5))")
                    logPDF(1:NP:NP/5) = getLogNormLogPDF(point(1:NP:NP/5), mu(1:NP:NP/5), Sigma(1:NP:NP/5))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example logPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "getLogNormLogPDF.RK.txt")
        do i = 1, NP
            write(fileUnit, "(*(E20.8e4,:,' '))") point(i) &
                                                , exp(getLogNormLogPDF(point(i), +0._RK, 2.00_RK)) &
                                                , exp(getLogNormLogPDF(point(i), +0._RK, 1.00_RK)) &
                                                , exp(getLogNormLogPDF(point(i), +0._RK, 0.50_RK)) &
                                                , exp(getLogNormLogPDF(point(i), +0._RK, 0.25_RK))
        end do
        close(fileUnit)
    end block

end program example