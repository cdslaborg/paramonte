program example

    use pm_kind, only: SK, IK
    use pm_distKolm, only: setKolmCDF
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: getLogSpace
    use pm_arrayResize, only: setResized
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: x(:), cdf(:)
        call disp%skip()
        call disp%show("x = [0._TKG, epsilon(0._TKG), 1._TKG, 2._TKG, 5._TKG]")
                        x = [0._TKG, epsilon(0._TKG), 1._TKG, 2._TKG, 5._TKG]
        call disp%show("x")
        call disp%show( x )
        call disp%show("call setResized(cdf, size(x, 1, IK))")
                        call setResized(cdf, size(x, 1, IK))
        call disp%show("call setKolmCDF(cdf, x)")
                        call setKolmCDF(cdf, x)
        call disp%show("cdf")
        call disp%show( cdf )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKH
        real(TKG), allocatable :: x(:), cdf(:)
        call disp%skip()
        call disp%show("x = [0._TKG, epsilon(0._TKG), 1._TKG, 2._TKG, 5._TKG]")
                        x = [0._TKG, epsilon(0._TKG), 1._TKG, 2._TKG, 5._TKG]
        call disp%show("x")
        call disp%show( x )
        call disp%show("call setResized(cdf, size(x, 1, IK))")
                        call setResized(cdf, size(x, 1, IK))
        call disp%show("call setKolmCDF(cdf, x)")
                        call setKolmCDF(cdf, x)
        call disp%show("cdf")
        call disp%show( cdf )
        call disp%skip()
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example cdf array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_io, only: getErrTableWrite
        use pm_kind, only: TKG => RKH
        real(TKG) :: xcdf(1000, 2)
        xcdf(:, 1) = getLinSpace(0._TKG, +3._TKG, size(xcdf, 1, IK))
        call setKolmCDF(xcdf(:, 2), xcdf(:, 1))
        if (0 /= getErrTableWrite(SK_"setKolmCDF.RK.txt", xcdf, header = SK_"x,cdf")) error stop "table output failed."
    end block

end program example