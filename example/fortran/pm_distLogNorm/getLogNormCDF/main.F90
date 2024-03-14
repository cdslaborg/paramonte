program example

    use pm_kind, only: SK
    use pm_kind, only: IK, RK ! all real kinds are supported.
    use pm_distLogNorm, only: getLogNormCDF
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    real(RK), dimension(NP) :: Point, Mu, Sigma, CDF

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(Mu, x1 = -5._RK, x2 = +5._RK)
    call setLogSpace(Point, logx1 = -10._RK, logx2 = +3._RK)
    call setLogSpace(Sigma, logx1 = log(0.1_RK), logx2 = log(10._RK))

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
    call disp%show("Point(NP/2)")
    call disp%show( Point(NP/2) )
    call disp%show("CDF(1) = getLogNormCDF(Point(NP/2))")
                    CDF(1) = getLogNormCDF(Point(NP/2))
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
    call disp%show("CDF(1) = getLogNormCDF(Point(1), Mu(1))")
                    CDF(1) = getLogNormCDF(Point(1), Mu(1))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! CDF with a standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Sigma(1)")
    call disp%show( Sigma(1) )
    call disp%show("Point(1)")
    call disp%show( Point(1) )
    call disp%show("CDF(1) = getLogNormCDF(Point(1), Sigma(1))")
                    CDF(1) = getLogNormCDF(Point(1), Sigma(1))
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
    call disp%show("Sigma(1)")
    call disp%show( Sigma(1) )
    call disp%show("Point(1)")
    call disp%show( Point(1) )
    call disp%show("CDF(1) = getLogNormCDF(Point(1), Mu(1), Sigma(1))")
                    CDF(1) = getLogNormCDF(Point(1), Mu(1), Sigma(1))
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
    call disp%show("Sigma(1)")
    call disp%show( Sigma(1) )
    call disp%show("Point(1:NP:NP/5)")
    call disp%show( Point(1:NP:NP/5) )
    call disp%show("CDF(1:NP:NP/5) = getLogNormCDF(Point(1:NP:NP/5), Mu(1), Sigma(1))")
                    CDF(1:NP:NP/5) = getLogNormCDF(Point(1:NP:NP/5), Mu(1), Sigma(1))
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
    call disp%show("Sigma(1:NP:NP/5)")
    call disp%show( Sigma(1:NP:NP/5) )
    call disp%show("Point(1)")
    call disp%show( Point(1) )
    call disp%show("CDF(1:NP:NP/5) = getLogNormCDF(Point(1), Mu(1:NP:NP/5), Sigma(1:NP:NP/5))")
                    CDF(1:NP:NP/5) = getLogNormCDF(Point(1), Mu(1:NP:NP/5), Sigma(1:NP:NP/5))
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
    call disp%show("Sigma(1:NP:NP/5)")
    call disp%show( Sigma(1:NP:NP/5) )
    call disp%show("Point(1:NP:NP/5)")
    call disp%show( Point(1:NP:NP/5) )
    call disp%show("CDF(1:NP:NP/5) = getLogNormCDF(Point(1:NP:NP/5), Mu(1:NP:NP/5), Sigma(1:NP:NP/5))")
                    CDF(1:NP:NP/5) = getLogNormCDF(Point(1:NP:NP/5), Mu(1:NP:NP/5), Sigma(1:NP:NP/5))
    call disp%show("CDF(1:NP:NP/5)")
    call disp%show( CDF(1:NP:NP/5) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example CDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "getLogNormCDF.RK.txt")
        write(fileUnit,"(5(g0,:,' '))") (Point(i), getLogNormCDF(Point(i), [0._RK, 0._RK, 0._RK, -2._RK], sigma = [3.0_RK, 1.0_RK, 0.3_RK, 1.0_RK]), i = 1, NP)
        close(fileUnit)
    end block

end program example