program example

    use pm_kind, only: SK
    use pm_kind, only: IK, RK ! all real kinds are supported.
    use pm_distLogNorm, only: setLogNormCDF
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    real(RK), dimension(NP) :: LogPoint, mu, invSigma, CDF

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(mu, x1 = -5._RK, x2 = +5._RK)
    call setLinSpace(LogPoint, x1 = -10._RK, x2 = +3._RK)
    call setLogSpace(invSigma, logx1 = log(0.1_RK), logx2 = log(10._RK))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of the (Standard) LogNormal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("! Standard CDF.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogPoint(NP/2)")
    call disp%show( LogPoint(NP/2) )
    call disp%show("call setLogNormCDF(CDF(1), LogPoint(NP/2))")
                    call setLogNormCDF(CDF(1), LogPoint(NP/2))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("! CDF with a mean.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("LogPoint(1)")
    call disp%show( LogPoint(1) )
    call disp%show("call setLogNormCDF(CDF(1), LogPoint(1), mu(1))")
                    call setLogNormCDF(CDF(1), LogPoint(1), mu(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! CDF with a standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("invSigma(1)")
    call disp%show( invSigma(1) )
    call disp%show("LogPoint(1)")
    call disp%show( LogPoint(1) )
    call disp%show("call setLogNormCDF(CDF(1), LogPoint(1), invSigma(1))")
                    call setLogNormCDF(CDF(1), LogPoint(1), invSigma(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! CDF with a mean and a standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("invSigma(1)")
    call disp%show( invSigma(1) )
    call disp%show("LogPoint(1)")
    call disp%show( LogPoint(1) )
    call disp%show("call setLogNormCDF(CDF(1), LogPoint(1), mu(1), invSigma(1))")
                    call setLogNormCDF(CDF(1), LogPoint(1), mu(1), invSigma(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at different points with the same mean and standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1)")
    call disp%show( mu(1) )
    call disp%show("invSigma(1)")
    call disp%show( invSigma(1) )
    call disp%show("LogPoint(1:NP:NP/5)")
    call disp%show( LogPoint(1:NP:NP/5) )
    call disp%show("call setLogNormCDF(CDF(1:NP:NP/5), LogPoint(1:NP:NP/5), mu(1), invSigma(1))")
                    call setLogNormCDF(CDF(1:NP:NP/5), LogPoint(1:NP:NP/5), mu(1), invSigma(1))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at the same point but with different means and standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1:NP:NP/5)")
    call disp%show( mu(1:NP:NP/5) )
    call disp%show("invSigma(1:NP:NP/5)")
    call disp%show( invSigma(1:NP:NP/5) )
    call disp%show("LogPoint(1)")
    call disp%show( LogPoint(1) )
    call disp%show("call setLogNormCDF(CDF(1:NP:NP/5), LogPoint(1), mu(1:NP:NP/5), invSigma(1:NP:NP/5))")
                    call setLogNormCDF(CDF(1:NP:NP/5), LogPoint(1), mu(1:NP:NP/5), invSigma(1:NP:NP/5))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at different points with different means and a standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mu(1:NP:NP/5)")
    call disp%show( mu(1:NP:NP/5) )
    call disp%show("invSigma(1:NP:NP/5)")
    call disp%show( invSigma(1:NP:NP/5) )
    call disp%show("LogPoint(1:NP:NP/5)")
    call disp%show( LogPoint(1:NP:NP/5) )
    call disp%show("call setLogNormCDF(CDF(1:NP:NP/5), LogPoint(1:NP:NP/5), mu(1:NP:NP/5), invSigma(1:NP:NP/5))")
                    call setLogNormCDF(CDF(1:NP:NP/5), LogPoint(1:NP:NP/5), mu(1:NP:NP/5), invSigma(1:NP:NP/5))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        real(RK)    , parameter :: mu(*) = [0._RK, 0._RK, 0._RK, -2._RK]
        real(RK)    , parameter :: invSigma(*) = 1._RK / [3.0_RK, 1.0_RK, 0.3_RK, 1.0_RK]
        open(newunit = fileUnit, file = "setLogNormCDF.RK.txt")
        do i = 1, NP
            call setLogNormCDF(CDF(1:4), LogPoint(i), mu, invSigma)
            write(fileUnit,"(5(g0,:,' '))") exp(LogPoint(i)), CDF(1:4)
        end do
        close(fileUnit)
    end block

end program example