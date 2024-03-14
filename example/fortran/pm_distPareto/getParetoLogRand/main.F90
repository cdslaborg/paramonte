program example

    use pm_kind, only: SK, IK, LK
    use pm_distPareto, only: getParetoLogRand
    use pm_io, only: display_type

    implicit none

    real                    :: LogX(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1) = getParetoLogRand(Alpha = -2., logMinX = -2.) ! Pareto distribution.")
                    LogX(1) = getParetoLogRand(Alpha = -2., logMinX = -2.) ! Pareto distribution.
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1:3) = getParetoLogRand(Alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.")
                    LogX(1:3) = getParetoLogRand(Alpha = -[+2., +3., +4.], logMinX = -2.) ! Pareto distribution.
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Truncated Pareto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1) = getParetoLogRand(Alpha = -1., logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.")
                    LogX(1) = getParetoLogRand(Alpha = -1., logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.
    call disp%show("LogX(1)")
    call disp%show( LogX(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogX(1:3) = getParetoLogRand(Alpha = -[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.")
                    LogX(1:3) = getParetoLogRand(Alpha = -[+1., +2., +3.], logMinX = -2., logMaxX = 5.) ! Truncated Pareto distribution.
    call disp%show("LogX(1:3)")
    call disp%show( LogX(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example LogX array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: Alpha(3), LogMinX(3), LogMaxX(3), LogX(3)
        integer(IK) :: fileUnit, i
        Alpha = -[.5, 1., 10.]
        LogMinX = log([3., 1., 2.])
        LogMaxX = log([5., 4., huge(0.)])
        open(newunit = fileUnit, file = "getParetoLogRand.RK.txt")
        do i = 1, 2000_IK
            LogX(1:2) = getParetoLogRand(Alpha(1:2), LogMinX(1:2), LogMaxX(1:2))
            LogX(3) = getParetoLogRand(Alpha(3), LogMinX(3))
            write(fileUnit, "(*(g0,:,', '))") exp(LogX)
        end do
        close(fileUnit)
    end block

end program example