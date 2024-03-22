program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_distUnif, only: getUnifRand
    use pm_distPoweto, only: getPowetoLogPDFNF

    implicit none

    real :: logPDFNF(3)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the natural logarithm of the normalization factor of the (Truncated) Poweto distribution PDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF(1) = getPowetoLogPDFNF(alpha = -1., logMinX = 1.)")
                    logPDFNF(1) = getPowetoLogPDFNF(alpha = -1., logMinX = 1.)
    call disp%show("logPDFNF(1)")
    call disp%show( logPDFNF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF(1:3) = getPowetoLogPDFNF(alpha = -2., logMinX = log([1., 2., 3.]))")
                    logPDFNF(1:3) = getPowetoLogPDFNF(alpha = -2., logMinX = log([1., 2., 3.]))
    call disp%show("logPDFNF(1:3)")
    call disp%show( logPDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF(1:3) = getPowetoLogPDFNF(alpha = -[+2., +3., +4.], logMinX = log([1., 2., 3.]))")
                    logPDFNF(1:3) = getPowetoLogPDFNF(alpha = -[+2., +3., +4.], logMinX = log([1., 2., 3.]))
    call disp%show("logPDFNF(1:3)")
    call disp%show( logPDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF(1:3) = getPowetoLogPDFNF(alpha = -[+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))")
                    logPDFNF(1:3) = getPowetoLogPDFNF(alpha = -[+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))
    call disp%show("logPDFNF(1:3)")
    call disp%show( logPDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF(1) = getPowetoLogPDFNF(alpha = 1., logMaxX = 1.)")
                    logPDFNF(1) = getPowetoLogPDFNF(alpha = 1., logMaxX = 1.)
    call disp%show("logPDFNF(1)")
    call disp%show( logPDFNF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF(1:3) = getPowetoLogPDFNF(alpha = 2., logMaxX = log([1., 2., 3.]))")
                    logPDFNF(1:3) = getPowetoLogPDFNF(alpha = 2., logMaxX = log([1., 2., 3.]))
    call disp%show("logPDFNF(1:3)")
    call disp%show( logPDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF(1:3) = getPowetoLogPDFNF(alpha = [+2., +3., +4.], logMaxX = log([1., 2., 3.]))")
                    logPDFNF(1:3) = getPowetoLogPDFNF(alpha = [+2., +3., +4.], logMaxX = log([1., 2., 3.]))
    call disp%show("logPDFNF(1:3)")
    call disp%show( logPDFNF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF(1:3) = getPowetoLogPDFNF(alpha = [+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))")
                    logPDFNF(1:3) = getPowetoLogPDFNF(alpha = [+2., +3., +4.], logMinX = log([1., 2., 3.]), logMaxX = log(20.))
    call disp%show("logPDFNF(1:3)")
    call disp%show( logPDFNF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        real, allocatable :: alpha(:), logPDFNF(:)
        alpha = getLinSpace(-10., +10., count = 500_IK)
        logPDFNF = getPowetoLogPDFNF(alpha, logMinX = -1., logMaxX = +1.)
        open(newunit = fileUnit, file = "getPowetoLogPDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (alpha(i), exp(logPDFNF(i)), i = 1, size(logPDFNF))
        close(fileUnit)
    end block

end program example