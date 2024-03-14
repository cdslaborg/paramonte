program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distLogUnif, only: setLogUnifPDF, getLogUnifPDFNF

    implicit none

    real :: PDF(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Probability Density Function (PDF) of the LogUniform distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogUnifPDF(PDF(1), x = 3., pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))")
                    call setLogUnifPDF(PDF(1), x = 3., pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))
    call disp%show("PDF(1)")
    call disp%show( PDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogUnifPDF(PDF(1:3), x = [3., 4., 2.5], pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))")
                    call setLogUnifPDF(PDF(1:3), x = [3., 4., 2.5], pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))
    call disp%show("PDF(1:3)")
    call disp%show( PDF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLogSpace
        real :: MinX(2), MaxX(2), X(2000), PDF(2), PDFNF(2)
        integer(IK) :: fileUnit, i
        call setLogSpace(X, logx1 = log(0.1), logx2 = log(13.))
        MinX = [3., 2.0]
        MaxX = [7., 10.]
        PDFNF = getLogUnifPDFNF(log(MinX), log(MaxX))
        open(newunit = fileUnit, file = "setLogUnifPDF.RK.txt")
        do i = 1, size(X)
            PDF = 0.
            if (MinX(1) <= X(i) .and. X(i) <= MaxX(1)) call setLogUnifPDF(PDF(1), X(i), PDFNF(1))
            if (MinX(2) <= X(i) .and. X(i) <= MaxX(2)) call setLogUnifPDF(PDF(2), X(i), PDFNF(2))
            write(fileUnit, "(*(g0,:,', '))") X(i), PDF
        end do
        close(fileUnit)
    end block

end program example