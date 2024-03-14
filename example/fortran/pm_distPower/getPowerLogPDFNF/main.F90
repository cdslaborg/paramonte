program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLogSpace
    use pm_distPower, only: getPowerLogPDFNF

    implicit none

    real :: LogPDFNF(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the (Truncated) Power distribution PDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDFNF(1) = getPowerLogPDFNF(alpha = 1., logMaxX = 1.)")
                    LogPDFNF(1) = getPowerLogPDFNF(alpha = 1., logMaxX = 1.)
    call disp%show("LogPDFNF(1)")
    call disp%show( LogPDFNF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDFNF(1:3) = getPowerLogPDFNF(alpha = 2., logMaxX = log([1., 2., 3.]))")
                    LogPDFNF(1:3) = getPowerLogPDFNF(alpha = 2., logMaxX = log([1., 2., 3.]))
    call disp%show("LogPDFNF(1:3)")
    call disp%show( LogPDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDFNF(1:3) = getPowerLogPDFNF(alpha = [+2., +3., +4.], logMaxX = log([1., 2., 3.]))")
                    LogPDFNF(1:3) = getPowerLogPDFNF(alpha = [+2., +3., +4.], logMaxX = log([1., 2., 3.]))
    call disp%show("LogPDFNF(1:3)")
    call disp%show( LogPDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("LogPDFNF(1:3) = getPowerLogPDFNF(alpha = [+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))")
                    LogPDFNF(1:3) = getPowerLogPDFNF(alpha = [+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))
    call disp%show("LogPDFNF(1:3)")
    call disp%show( LogPDFNF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example PDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        real, allocatable :: Alpha(:), LogPDFNF(:)
        Alpha = getLogSpace(-10., +10., count = 500_IK)
        LogPDFNF = getPowerLogPDFNF(Alpha, logMinX = -1., logMaxX = +1.)
        open(newunit = fileUnit, file = "getPowerLogPDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (Alpha(i), exp(LogPDFNF(i)), i = 1, size(LogPDFNF))
        close(fileUnit)
    end block

end program example