program example

    use pm_kind, only: SK, IK
    use pm_distKolm, only: getKolmPDF
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: x(:), pdf(:)
        call disp%skip()
        call disp%show("x = [0._TKC, epsilon(0._TKC), 1._TKC, 2._TKC, 5._TKC]")
                        x = [0._TKC, epsilon(0._TKC), 1._TKC, 2._TKC, 5._TKC]
        call disp%show("x")
        call disp%show( x )
        call disp%show("pdf = getKolmPDF(x)")
                        pdf = getKolmPDF(x)
        call disp%show("pdf")
        call disp%show( pdf )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKH
        real(TKC), allocatable :: x(:), pdf(:)
        call disp%skip()
        call disp%show("x = [0._TKC, epsilon(0._TKC), 1._TKC, 2._TKC, 5._TKC]")
                        x = [0._TKC, epsilon(0._TKC), 1._TKC, 2._TKC, 5._TKC]
        call disp%show("x")
        call disp%show( x )
        call disp%show("pdf = getKolmPDF(x)")
                        pdf = getKolmPDF(x)
        call disp%show("pdf")
        call disp%show( pdf )
        call disp%skip()
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example pdf array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_io, only: getErrTableWrite
        use pm_kind, only: TKC => RKH
        real(TKC) :: xpdf(1000, 2)
        xpdf(:, 1) = getLinSpace(0._TKC, +3._TKC, size(xpdf, 1, IK))
        xpdf(:, 2) = getKolmPDF(xpdf(:, 1))
        if (0 /= getErrTableWrite(SK_"getKolmPDF.RK.txt", xpdf, header = SK_"x,pdf")) error stop "table output failed."
    end block

end program example