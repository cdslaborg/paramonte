program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLogSpace
    use pm_distPower, only: getPowerLogCDFNF

    implicit none

    real :: logCDFNF(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the (Truncated) Power distribution CDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logCDFNF(1) = getPowerLogCDFNF(alpha = 1., logMaxX = 1.)")
                    logCDFNF(1) = getPowerLogCDFNF(alpha = 1., logMaxX = 1.)
    call disp%show("logCDFNF(1)")
    call disp%show( logCDFNF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logCDFNF(1:3) = getPowerLogCDFNF(alpha = 2., logMaxX = log([1., 2., 3.]))")
                    logCDFNF(1:3) = getPowerLogCDFNF(alpha = 2., logMaxX = log([1., 2., 3.]))
    call disp%show("logCDFNF(1:3)")
    call disp%show( logCDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logCDFNF(1:3) = getPowerLogCDFNF(alpha = [+2., +3., +4.], logMaxX = log([1., 2., 3.]))")
                    logCDFNF(1:3) = getPowerLogCDFNF(alpha = [+2., +3., +4.], logMaxX = log([1., 2., 3.]))
    call disp%show("logCDFNF(1:3)")
    call disp%show( logCDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logCDFNF(1:3) = getPowerLogCDFNF(alpha = [+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))")
                    logCDFNF(1:3) = getPowerLogCDFNF(alpha = [+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))
    call disp%show("logCDFNF(1:3)")
    call disp%show( logCDFNF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        real, allocatable :: alpha(:), logCDFNF(:)
        alpha = getLogSpace(-10., +10., count = 500_IK)
        logCDFNF = getPowerLogCDFNF(alpha, logMinX = -1., logMaxX = +1.)
        open(newunit = fileUnit, file = "getPowerLogCDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (alpha(i), exp(logCDFNF(i)), i = 1, size(logCDFNF))
        close(fileUnit)
    end block

end program example