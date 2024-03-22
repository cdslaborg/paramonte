program example

    use pm_kind, only: SK, IK, LK
    use pm_distPoweto, only: getPowetoLogCDFNF
    use pm_distPoweto, only: setPowetoLogCDF
    use pm_io, only: display_type

    implicit none

    real                    :: logCDF(3)

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of the Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogCDF(logCDF(1), logx = 3., alpha = -2., logMinX = -2.) ! Poweto distribution.")
                    call setPowetoLogCDF(logCDF(1), logx = 3., alpha = -2., logMinX = -2.) ! Poweto distribution.
    call disp%show("logCDF(1)")
    call disp%show( logCDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogCDF(logCDF(1:3), logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2.) ! Poweto distribution.")
                    call setPowetoLogCDF(logCDF(1:3), logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2.) ! Poweto distribution.
    call disp%show("logCDF(1:3)")
    call disp%show( logCDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogCDF(logCDF(1), logx = 3., alpha = +2., logCDFNF = getPowetoLogCDFNF(alpha = +2., logMaxX = 5.)) ! Poweto distribution.")
                    call setPowetoLogCDF(logCDF(1), logx = 3., alpha = +2., logCDFNF = getPowetoLogCDFNF(alpha = +2., logMaxX = 5.)) ! Poweto distribution.
    call disp%show("logCDF(1)")
    call disp%show( logCDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogCDF(logCDF(1:3), logx = [3., 4., 5.], alpha = [+2., +3., +4.], logCDFNF = getPowetoLogCDFNF(alpha = [+2., +3., +4.], logMaxX = 5.)) ! Poweto distribution.")
                    call setPowetoLogCDF(logCDF(1:3), logx = [3., 4., 5.], alpha = [+2., +3., +4.], logCDFNF = getPowetoLogCDFNF(alpha = [+2., +3., +4.], logMaxX = 5.)) ! Poweto distribution.
    call disp%show("logCDF(1:3)")
    call disp%show( logCDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of the Truncated Poweto distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogCDF(logCDF(1), logx = 3., alpha = -2., logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = -2., logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.")
                    call setPowetoLogCDF(logCDF(1), logx = 3., alpha = -2., logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = -2., logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.
    call disp%show("logCDF(1)")
    call disp%show( logCDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogCDF(logCDF(1:3), logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.")
                    call setPowetoLogCDF(logCDF(1:3), logx = [3., 4., 5.], alpha = -[+2., +3., +4.], logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = -[+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.
    call disp%show("logCDF(1:3)")
    call disp%show( logCDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogCDF(logCDF(1), logx = 3., alpha = +2., logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.")
                    call setPowetoLogCDF(logCDF(1), logx = 3., alpha = +2., logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = +2., logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.
    call disp%show("logCDF(1)")
    call disp%show( logCDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPowetoLogCDF(logCDF(1:3), logx = [3., 4., 5.], alpha = [+2., +3., +4.], logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = [+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.")
                    call setPowetoLogCDF(logCDF(1:3), logx = [3., 4., 5.], alpha = [+2., +3., +4.], logMinX = -2., logCDFNF = getPowetoLogCDFNF(alpha = [+2., +3., +4.], logMinX = -2., logMaxX = 5.)) ! Truncated Poweto distribution.
    call disp%show("logCDF(1:3)")
    call disp%show( logCDF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        use pm_arrayResize, only: setResized
        real, allocatable :: alpha(:), logCDF(:)
        real :: logMinX, logMaxX, logx(2000)
        integer(IK) :: fileUnit, i
        alpha = [-2.0, -1., 0., +1., +2.]
        call setLinSpace(logx, x1 = log(0.1), x2 = log(10.))
        call setResized(logCDF, size(alpha, 1, IK))
        logMinX = log(3.)
        logMaxX = log(8.)
        open(newunit = fileUnit, file = "setPowetoLogCDF.RK.txt")
        do i = 1, size(logx, 1, IK)
            if (logx(i) < logMinX) then
                logCDF = -huge(0.)
            elseif (logx(i) < logMaxX) then
                call setPowetoLogCDF(logCDF, logx(i), alpha, logMinX, getPowetoLogCDFNF(alpha, logMinX, logMaxX))
            else
                logCDF = 0.
            end if
            write(fileUnit, "(*(g0,:,', '))") exp(logx(i)), exp(logCDF)
        end do
        close(fileUnit)
    end block

end program example