program example

    use pm_kind, only: SK, IK, LK
    use pm_distPoweto, only: getPowetoLogQuan
    use pm_io, only: display_type

    implicit none

    real                    :: logx(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = getPowetoLogQuan(logCDF = -3., alpha = -2., logMinX = -2.) ! Poweto distribution.")
                    logx(1) = getPowetoLogQuan(logCDF = -3., alpha = -2., logMinX = -2.) ! Poweto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1:3) = getPowetoLogQuan(logCDF = -[3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2.) ! Poweto distribution.")
                    logx(1:3) = getPowetoLogQuan(logCDF = -[3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2.) ! Poweto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Quantile Function of the Truncated Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = getPowetoLogQuan(logCDF = -3., alpha = -2., logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.")
                    logx(1) = getPowetoLogQuan(logCDF = -3., alpha = -2., logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1:3) = getPowetoLogQuan(logCDF = -[3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.")
                    logx(1:3) = getPowetoLogQuan(logCDF = -[3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        use pm_arrayResize, only: setResized
        real, allocatable :: alpha(:), logx(:)
        real :: logMinX, logMaxX, logCDF(2000)
        integer(IK) :: fileUnit, i
        alpha = [-2.0, -1., 0., +1., +2.]
        call setLinSpace(logCDF, x1 = log(0.001), x2 = log(1.))
        call setResized(logx, size(alpha, 1, IK))
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "getPowetoLogQuan.RK.txt")
        do i = 1, size(logCDF, 1, IK)
            logx = getPowetoLogQuan(logCDF(i), alpha, logMinX, logMaxX)
            write(fileUnit, "(*(g0,:,', '))") exp(logCDF(i)), exp(logx)
        end do
        close(fileUnit)
    end block

end program example