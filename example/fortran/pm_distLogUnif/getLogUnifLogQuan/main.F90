program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distLogUnif, only: getLogUnifLogQuan

    implicit none

    real :: LogX(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Inverse CDF (Quantile) Function of the LogUniform distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1) = getLogUnifLogQuan(cdf = 0.5, logMinX = log(2.), logMaxX = log(5.))")
                    LogX(1) = getLogUnifLogQuan(cdf = 0.5, logMinX = log(2.), logMaxX = log(5.))
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1:3) = getLogUnifLogQuan(cdf = [.3, .4, .25], logMinX = log(2.), logMaxX = log(5.))")
                    LogX(1:3) = getLogUnifLogQuan(cdf = [.3, .4, .25], logMinX = log(2.), logMaxX = log(5.))
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: LogMinX(2), LogMaxX(2), LogX(2), CDF(2000)
        integer(IK) :: fileUnit, i
        call setLinSpace(CDF, x1 = 0., x2 = 1.)
        LogMinX = log([3., 2.0])
        LogMaxX = log([7., 10.])
        open(newunit = fileUnit, file = "getLogUnifLogQuan.RK.txt")
        do i = 1, size(CDF)
            LogX = getLogUnifLogQuan(CDF(i), LogMinX, LogMaxX)
            write(fileUnit, "(*(g0,:,', '))") CDF(i), exp(LogX)
        end do
        close(fileUnit)
    end block

end program example