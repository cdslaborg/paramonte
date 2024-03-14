program example

    use pm_kind, only: SK
    use pm_kind, only: IK, RK ! all real kinds are supported.
    use pm_distNorm, only: setNormCDF
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    real(RK), dimension(NP) :: Point, Mu, InvSigma, CDF

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(Mu, x1 = -5._RK, x2 = +5._RK)
    call setLinSpace(Point, x1 = -10._RK, x2 = +10._RK)
    call setLogSpace(InvSigma, logx1 = log(0.1_RK), logx2 = log(10._RK))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of the (Standard) Normal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("! Standard CDF.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Point(NP/2)")
    call disp%show( Point(NP/2) )
    call disp%show("call setNormCDF(CDF(1), Point(NP/2))")
                    call setNormCDF(CDF(1), Point(NP/2))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("! CDF with a mean.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Mu(1)")
    call disp%show( Mu(1) )
    call disp%show("Point(1)")
    call disp%show( Point(1) )
    call disp%show("call setNormCDF(CDF(1), Point(1), Mu(1))")
                    call setNormCDF(CDF(1), Point(1), Mu(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! CDF with a standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("InvSigma(1)")
    call disp%show( InvSigma(1) )
    call disp%show("Point(1)")
    call disp%show( Point(1) )
    call disp%show("call setNormCDF(CDF(1), Point(1), InvSigma(1))")
                    call setNormCDF(CDF(1), Point(1), InvSigma(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! CDF with a mean and a standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Mu(1)")
    call disp%show( Mu(1) )
    call disp%show("InvSigma(1)")
    call disp%show( InvSigma(1) )
    call disp%show("Point(1)")
    call disp%show( Point(1) )
    call disp%show("call setNormCDF(CDF(1), Point(1), Mu(1), InvSigma(1))")
                    call setNormCDF(CDF(1), Point(1), Mu(1), InvSigma(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at different points with the same mean and standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Mu(1)")
    call disp%show( Mu(1) )
    call disp%show("InvSigma(1)")
    call disp%show( InvSigma(1) )
    call disp%show("Point(1:NP:NP/5)")
    call disp%show( Point(1:NP:NP/5) )
    call disp%show("call setNormCDF(CDF(1:NP:NP/5), Point(1:NP:NP/5), Mu(1), InvSigma(1))")
                    call setNormCDF(CDF(1:NP:NP/5), Point(1:NP:NP/5), Mu(1), InvSigma(1))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at the same point but with different means and standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Mu(1:NP:NP/5)")
    call disp%show( Mu(1:NP:NP/5) )
    call disp%show("InvSigma(1:NP:NP/5)")
    call disp%show( InvSigma(1:NP:NP/5) )
    call disp%show("Point(1)")
    call disp%show( Point(1) )
    call disp%show("call setNormCDF(CDF(1:NP:NP/5), Point(1), Mu(1:NP:NP/5), InvSigma(1:NP:NP/5))")
                    call setNormCDF(CDF(1:NP:NP/5), Point(1), Mu(1:NP:NP/5), InvSigma(1:NP:NP/5))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! A vector of CDF at different points with different means and a standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Mu(1:NP:NP/5)")
    call disp%show( Mu(1:NP:NP/5) )
    call disp%show("InvSigma(1:NP:NP/5)")
    call disp%show( InvSigma(1:NP:NP/5) )
    call disp%show("Point(1:NP:NP/5)")
    call disp%show( Point(1:NP:NP/5) )
    call disp%show("call setNormCDF(CDF(1:NP:NP/5), Point(1:NP:NP/5), Mu(1:NP:NP/5), InvSigma(1:NP:NP/5))")
                    call setNormCDF(CDF(1:NP:NP/5), Point(1:NP:NP/5), Mu(1:NP:NP/5), InvSigma(1:NP:NP/5))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example CDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        real(RK)    , parameter :: Mu(*) = [0._RK, 0._RK, 0._RK, -2._RK]
        real(RK)    , parameter :: InvSigma(*) = 1._RK / [3.0_RK, 1.0_RK, 0.3_RK, 1.0_RK]
        open(newunit = fileUnit, file = "setNormCDF.RK.txt")
        do i = 1, NP
            call setNormCDF(CDF(1:4), Point(i), Mu, InvSigma)
            write(fileUnit,"(5(g0,:,' '))") Point(i), CDF(1:4)
        end do
        close(fileUnit)
    end block

end program example