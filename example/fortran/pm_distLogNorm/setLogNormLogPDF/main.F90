program example

    use pm_kind, only: SK
    use pm_kind, only: IK, RK ! all real kinds are supported.
    use pm_distLogNorm, only: setLogNormLogPDF
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    real(RK), dimension(NP) :: logx, mu, invSigma, logPDF

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(mu, x1 = -5._RK, x2 = +5._RK)
    call setLinSpace(logx, x1 = log(0.0001_RK), x2 = log(5._RK))
    call setLogSpace(invSigma, logx1 = log(0.1_RK), logx2 = log(10._RK))

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
    call disp%show("logx(NP/2)")
    call disp%show( logx(NP/2) )
    call disp%show("call setLogNormLogPDF(logPDF(1), logx(NP/2))")
                    call setLogNormLogPDF(logPDF(1), logx(NP/2))
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
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%show("call setLogNormLogPDF(logPDF(1), logx(1), mu(1))")
                    call setLogNormLogPDF(logPDF(1), logx(1), mu(1))
    call disp%show("logPDF(1)")
    call disp%show( logPDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! PDF with a standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("invSigma(1)")
    call disp%show( invSigma(1) )
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%show("call setLogNormLogPDF(logPDF(1), logx(1), invSigma(1), log(invSigma(1)))")
                    call setLogNormLogPDF(logPDF(1), logx(1), invSigma(1), log(invSigma(1)))
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
    call disp%show("invSigma(1)")
    call disp%show( invSigma(1) )
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%show("call setLogNormLogPDF(logPDF(1), logx(1), mu(1), invSigma(1), log(invSigma(1)))")
                    call setLogNormLogPDF(logPDF(1), logx(1), mu(1), invSigma(1), log(invSigma(1)))
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
    call disp%show("invSigma(1)")
    call disp%show( invSigma(1) )
    call disp%show("logx(1:NP:NP/5)")
    call disp%show( logx(1:NP:NP/5) )
    call disp%show("call setLogNormLogPDF(logPDF(1:NP:NP/5), logx(1:NP:NP/5), mu(1), invSigma(1), log(invSigma(1)))")
                    call setLogNormLogPDF(logPDF(1:NP:NP/5), logx(1:NP:NP/5), mu(1), invSigma(1), log(invSigma(1)))
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
    call disp%show("invSigma(1:NP:NP/5)")
    call disp%show( invSigma(1:NP:NP/5) )
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%show("call setLogNormLogPDF(logPDF(1:NP:NP/5), logx(1), mu(1:NP:NP/5), invSigma(1:NP:NP/5), log(invSigma(1:NP:NP/5)))")
                    call setLogNormLogPDF(logPDF(1:NP:NP/5), logx(1), mu(1:NP:NP/5), invSigma(1:NP:NP/5), log(invSigma(1:NP:NP/5)))
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
    call disp%show("invSigma(1:NP:NP/5)")
    call disp%show( invSigma(1:NP:NP/5) )
    call disp%show("logx(1:NP:NP/5)")
    call disp%show( logx(1:NP:NP/5) )
    call disp%show("call setLogNormLogPDF(logPDF(1:NP:NP/5), logx(1:NP:NP/5), mu(1:NP:NP/5), invSigma(1:NP:NP/5), log(invSigma(1:NP:NP/5)))")
                    call setLogNormLogPDF(logPDF(1:NP:NP/5), logx(1:NP:NP/5), mu(1:NP:NP/5), invSigma(1:NP:NP/5), log(invSigma(1:NP:NP/5)))
    call disp%show("logPDF(1:NP:NP/5)")
    call disp%show( logPDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example logPDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        real(RK), parameter :: invSigma(*) = 1._RK / [2.00_RK, 1.00_RK, 0.50_RK, 0.25_RK]
        open(newunit = fileUnit, file = "setLogNormLogPDF.RK.txt")
        do i = 1, NP
            call setLogNormLogPDF(logPDF(1:4), logx(i), +0._RK, invSigma, log(invSigma))
            write(fileUnit, "(*(E20.8e4,:,' '))") exp(logx(i)), exp(logPDF(1:4))
        end do
        close(fileUnit)
    end block

end program example