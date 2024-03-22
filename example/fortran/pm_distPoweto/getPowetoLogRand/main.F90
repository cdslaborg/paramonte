program example

    use pm_kind, only: SK, IK, LK
    use pm_distPoweto, only: getPowetoLogCDFNF
    use pm_distPoweto, only: getPowetoLogRand
    use pm_io, only: display_type

    implicit none

    real                    :: logx(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = getPowetoLogRand(alpha = -2., logMinX = -2.) ! Poweto distribution.")
                    logx(1) = getPowetoLogRand(alpha = -2., logMinX = -2.) ! Poweto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1:3) = getPowetoLogRand(alpha = -[+2., +3., +4.], logMinX = -2.) ! Poweto distribution.")
                    logx(1:3) = getPowetoLogRand(alpha = -[+2., +3., +4.], logMinX = -2.) ! Poweto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = getPowetoLogRand(alpha = +2., logMaxX = -2.) ! Poweto distribution.")
                    logx(1) = getPowetoLogRand(alpha = +2., logMaxX = -2.) ! Poweto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1:3) = getPowetoLogRand(alpha = [+2., +3., +4.], logMaxX = -2.) ! Poweto distribution.")
                    logx(1:3) = getPowetoLogRand(alpha = [+2., +3., +4.], logMaxX = -2.) ! Poweto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Truncated Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = getPowetoLogRand(alpha = -2., logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.")
                    logx(1) = getPowetoLogRand(alpha = -2., logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1:3) = getPowetoLogRand(alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.")
                    logx(1:3) = getPowetoLogRand(alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = getPowetoLogRand(alpha = +2., logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.")
                    logx(1) = getPowetoLogRand(alpha = +2., logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1:3) = getPowetoLogRand(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.")
                    logx(1:3) = getPowetoLogRand(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.) ! Truncated Poweto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        use pm_arrayResize, only: setResized
        use pm_arrayRefill, only: setRefilled
        real, allocatable :: alpha(:), sigma(:)
        real :: logMinX, logMaxX
        integer(IK) :: fileUnit, i
        alpha = [-2.0, -1., 0., +1., +2.]
        call setRefilled(sigma, 1., size(alpha, 1, IK))
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "getPowetoLogRand.RK.txt")
        do i = 1, 2000
            write(fileUnit, "(*(g0,:,', '))") exp(getPowetoLogRand(alpha, logMinX, logMaxX))
        end do
        close(fileUnit)
    end block

end program example