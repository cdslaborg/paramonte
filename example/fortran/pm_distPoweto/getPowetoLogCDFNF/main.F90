program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_distPoweto, only: getPowetoLogCDFNF

    implicit none

    real :: logCDFNF(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the (Truncated) Poweto distribution CDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logCDFNF(1:3) = getPowetoLogCDFNF(alpha = -[+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))")
                    logCDFNF(1:3) = getPowetoLogCDFNF(alpha = -[+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))
    call disp%show("logCDFNF(1:3)")
    call disp%show( logCDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logCDFNF(1) = getPowetoLogCDFNF(alpha = 1., logMaxX = 1.)")
                    logCDFNF(1) = getPowetoLogCDFNF(alpha = 1., logMaxX = 1.)
    call disp%show("logCDFNF(1)")
    call disp%show( logCDFNF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logCDFNF(1:3) = getPowetoLogCDFNF(alpha = 2., logMaxX = log([1., 2., 3.]))")
                    logCDFNF(1:3) = getPowetoLogCDFNF(alpha = 2., logMaxX = log([1., 2., 3.]))
    call disp%show("logCDFNF(1:3)")
    call disp%show( logCDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logCDFNF(1:3) = getPowetoLogCDFNF(alpha = [+2., +3., +4.], logMaxX = log([1., 2., 3.]))")
                    logCDFNF(1:3) = getPowetoLogCDFNF(alpha = [+2., +3., +4.], logMaxX = log([1., 2., 3.]))
    call disp%show("logCDFNF(1:3)")
    call disp%show( logCDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logCDFNF(1:3) = getPowetoLogCDFNF(alpha = [+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))")
                    logCDFNF(1:3) = getPowetoLogCDFNF(alpha = [+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))
    call disp%show("logCDFNF(1:3)")
    call disp%show( logCDFNF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        real, allocatable :: alpha(:), logCDFNF(:)
        alpha = 5 * getLinSpace(-1., +1., count = 500_IK); !alpha(size(alpha) + 1) = 0.
        logCDFNF = getPowetoLogCDFNF(alpha, logMinX = -1., logMaxX = +1.)
        open(newunit = fileUnit, file = "getPowetoLogCDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (alpha(i), exp(logCDFNF(i)), i = 1, size(logCDFNF))
        close(fileUnit)
    end block

end program example