program example

    use pm_kind, only: SK, IK, LK
    use pm_distPoweto, only: getPowetoLogCDFNF
    use pm_distPoweto, only: setPowetoLogRand
    use pm_distNegExp, only: getNegExpRand
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
    call disp%show("call setPowetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = -2., logMinX = -2.) ! Poweto distribution.")
                    call setPowetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = -2., logMinX = -2.) ! Poweto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = -[+2., +3., +4.], logMinX = -2.) ! Poweto distribution.")
                    call setPowetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = -[+2., +3., +4.], logMinX = -2.) ! Poweto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = +2., logCDFNF = getPowetoLogCDFNF(alpha = +2., logMaxX = -2.)) ! Poweto distribution.")
                    call setPowetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = +2., logCDFNF = getPowetoLogCDFNF(alpha = +2., logMaxX = -2.)) ! Poweto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = [+2., +3., +4.], logCDFNF = getPowetoLogCDFNF(alpha = [+2., +3., +4.], logMaxX = -2.)) ! Poweto distribution.")
                    call setPowetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = [+2., +3., +4.], logCDFNF = getPowetoLogCDFNF(alpha = [+2., +3., +4.], logMaxX = -2.)) ! Poweto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute random value(s) from the Truncated Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = -2., logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = -2., logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.")
                    call setPowetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = -2., logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = -2., logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = -[+2., +3., +4.], logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.")
                    call setPowetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = -[+2., +3., +4.], logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = +2., logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.")
                    call setPowetoLogRand(logx(1), negExpRand = getNegExpRand(1.), alpha = +2., logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = +[+2., +3., +4.], logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.")
                    call setPowetoLogRand(logx(1:3), negExpRand = getNegExpRand([1., 1., 1.]), alpha = +[+2., +3., +4.], logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = +[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.
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
        real, allocatable :: alpha(:), logx(:), sigma(:)
        real :: logMinX, logMaxX
        integer(IK) :: fileUnit, i
        alpha = [-2.0, -1., 0., +1., +2.]
        call setRefilled(sigma, 1., size(alpha, 1, IK))
        call setResized(logx, size(alpha, 1, IK))
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "setPowetoLogRand.RK.txt")
        do i = 1, 2000
            call setPowetoLogRand(logx, getNegExpRand(1.), alpha, logMinX, getPowetoLogCDFNF(alpha, logMinX, logMaxX))
            write(fileUnit, "(*(g0,:,', '))") exp(logx)
        end do
        close(fileUnit)
    end block

end program example