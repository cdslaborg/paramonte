program example

    use pm_kind, only: SK, IK, LK
    use pm_distPareto, only: getParetoLogRand
    use pm_io, only: display_type

    implicit none

    real                    :: logx(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = getParetoLogRand(alpha = -2., logMinX = -2.) ! Pareto distribution.")
                    logx(1) = getParetoLogRand(alpha = -2., logMinX = -2.) ! Pareto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1:3) = getParetoLogRand(alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.")
                    logx(1:3) = getParetoLogRand(alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Truncated Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1) = getParetoLogRand(alpha = -1., logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.")
                    logx(1) = getParetoLogRand(alpha = -1., logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logx(1:3) = getParetoLogRand(alpha = -[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.")
                    logx(1:3) = getParetoLogRand(alpha = -[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: alpha(3), logMinX(3), logMaxX(3), logx(3)
        integer(IK) :: fileUnit, i
        alpha = -[.5, 1., 10.]
        logMinX = log([3., 1., 2.])
        logMaxX = log([5., 4., huge(0.)])
        open(newunit = fileUnit, file = "getParetoLogRand.RK.txt")
        do i = 1, 2000_IK
            logx(1:2) = getParetoLogRand(alpha(1:2), logMinX(1:2), logMaxX(1:2))
            logx(3) = getParetoLogRand(alpha(3), logMinX(3))
            write(fileUnit, "(*(g0,:,', '))") exp(logx)
        end do
        close(fileUnit)
    end block

end program example