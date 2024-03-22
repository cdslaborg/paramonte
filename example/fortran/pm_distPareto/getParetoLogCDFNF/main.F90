program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLogSpace
    use pm_distPareto, only: getParetoLogCDFNF

    implicit none

    real :: logCDFNF(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the (Truncated) Pareto distribution CDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logCDFNF(1:3) = getParetoLogCDFNF(alpha = -[+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))")
                    logCDFNF(1:3) = getParetoLogCDFNF(alpha = -[+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))
    call disp%show("logCDFNF(1:3)")
    call disp%show( logCDFNF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        real, allocatable :: alpha(:), logCDFNF(:)
        alpha = -getLogSpace(-10., +10., count = 500_IK)
        logCDFNF = getParetoLogCDFNF(alpha, logMinX = -1., logMaxX = +1.)
        open(newunit = fileUnit, file = "getParetoLogCDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (-alpha(i), exp(logCDFNF(i)), i = 1, size(logCDFNF))
        close(fileUnit)
    end block

end program example