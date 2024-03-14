program example

    use pm_kind, only: IK
    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_distLogUnif, only: getLogUnifPDFNF

    implicit none

    real :: PDFNF(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the normalization factor of the LogUniform distribution PDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("PDFNF(1) = getLogUnifPDFNF(logMinX = 1., logMaxX = 3.)")
                    PDFNF(1) = getLogUnifPDFNF(logMinX = 1., logMaxX = 3.)
    call disp%show("PDFNF(1)")
    call disp%show( PDFNF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("PDFNF(1:3) = getLogUnifPDFNF(logMinX = log([1., 2., 3.]), logMaxX = log(20.))")
                    PDFNF(1:3) = getLogUnifPDFNF(logMinX = log([1., 2., 3.]), logMaxX = log(20.))
    call disp%show("PDFNF(1:3)")
    call disp%show( PDFNF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example PDF array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer :: fileUnit, i
        real, allocatable :: LogMaxX(:), PDFNF(:)
        LogMaxX = getLinSpace(0.01, 10., count = 500_IK)
        PDFNF = getLogUnifPDFNF(logMinX = 0., logMaxX = LogMaxX)
        open(newunit = fileUnit, file = "getLogUnifPDFNF.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (exp(LogMaxX(i)), PDFNF(i), i = 1, size(PDFNF))
        close(fileUnit)
    end block

end program example