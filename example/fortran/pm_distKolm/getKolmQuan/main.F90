program example

    use pm_kind, only: SK, IK
    use pm_distKolm, only: getKolmCDF
    use pm_distKolm, only: getKolmQuan
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: x(:), cdf(:)
        call disp%skip()
        call disp%show("x = [0._TKC, epsilon(0._TKC), 1._TKC, 2._TKC]")
                        x = [0._TKC, epsilon(0._TKC), 1._TKC, 2._TKC]
        call disp%show("x")
        call disp%show( x )
        call disp%show("cdf = getKolmCDF(x)")
                        cdf = getKolmCDF(x)
        call disp%show("cdf")
        call disp%show( cdf )
        call disp%show("x = getKolmQuan(cdf)")
                        x = getKolmQuan(cdf)
        call disp%show("x")
        call disp%show( x )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKH
        real(TKC), allocatable :: x(:), cdf(:)
        call disp%skip()
        call disp%show("x = [0._TKC, epsilon(0._TKC), 1._TKC, 2._TKC]")
                        x = [0._TKC, epsilon(0._TKC), 1._TKC, 2._TKC]
        call disp%show("x")
        call disp%show( x )
        call disp%show("cdf = getKolmCDF(x)")
                        cdf = getKolmCDF(x)
        call disp%show("cdf")
        call disp%show( cdf )
        call disp%show("x = getKolmQuan(cdf)")
                        x = getKolmQuan(cdf)
        call disp%show("x")
        call disp%show( x )
        call disp%skip()
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_io, only: getErrTableWrite
        use pm_kind, only: TKC => RKH
        real(TKC) :: cdfx(1500, 2)
        cdfx(:, 1) = getLinSpace(0._TKC, +.9999_TKC, size(cdfx, 1, IK))
        cdfx(:, 2) = getKolmQuan(cdfx(:, 1))
        if (0 /= getErrTableWrite(SK_"getKolmQuan.RK.txt", cdfx, header = SK_"cdf,quantile")) error stop "table output failed."
    end block

end program example