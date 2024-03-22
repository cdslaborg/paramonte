program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_distNegExp, only: getNegExpCDF

    implicit none

    real :: CDF(3)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of Negative Exponential distribution at the specified values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the CDF at an input scalar real value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("CDF(1) = getNegExpCDF(x = -.5)")
                    CDF(1) = getNegExpCDF(x = -.5)
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("CDF(1) = getNegExpCDF(x = -.5, mu = 0.)")
                    CDF(1) = getNegExpCDF(x = -.5, mu = 0.)
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("CDF(1) = getNegExpCDF(x = -.5, invSigma = 1.)")
                    CDF(1) = getNegExpCDF(x = -.5, invSigma = 1.)
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("CDF(1) = getNegExpCDF(x = -.5, mu = 0., invSigma = 1.)")
                    CDF(1) = getNegExpCDF(x = -.5, mu = 0., invSigma = 1.)
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the CDF at an input vector real value with different parameter values.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("CDF(1:3) = getNegExpCDF(x = -[.5, 1., 2.])")
                    CDF(1:3) = getNegExpCDF(x = -[.5, 1., 2.])
    call disp%show("CDF(1:3)")
    call disp%show( CDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("CDF(1:3) = getNegExpCDF(x = -[.5, 1., 2.], mu = 0.)")
                    CDF(1:3) = getNegExpCDF(x = -[.5, 1., 2.], mu = 0.)
    call disp%show("CDF(1:3)")
    call disp%show( CDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("CDF(1:3) = getNegExpCDF(x = -[.5, 1., 2.], invSigma = [2., 1., .5])")
                    CDF(1:3) = getNegExpCDF(x = -[.5, 1., 2.], invSigma = [2., 1., .5])
    call disp%show("CDF(1:3)")
    call disp%show( CDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("CDF(1:3) = getNegExpCDF(x = -[5., 1., 2.], mu = [-1., 0., +1.], invSigma = 1.)")
                    CDF(1:3) = getNegExpCDF(x = -[5., 1., 2.], mu = [-1., 0., +1.], invSigma = 1.)
    call disp%show("CDF(1:3)")
    call disp%show( CDF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: getLinSpace
        real, parameter :: invSigma(3) = [2., 1., .5]
        real, parameter :: mu(3) = [-2., 0., 2.]
        real :: Point(1000)
        integer :: fileUnit, i, j
        Point = -getLinSpace(-4., +8., size(Point,1,IK))
        open(newunit = fileUnit, file = "getNegExpCDF.RK.txt")
        do i = 1, size(Point)
            do j = 1, size(CDF)
                CDF(j) = 0.
                if (Point(i) <= mu(j)) CDF(j) = getNegExpCDF(Point(i), mu(j), invSigma(j))
            end do
            write(fileUnit,"(*(g0,:,', '))") Point(i), CDF
        end do
        close(fileUnit)
    end block

end program example