program example

    use pm_kind, only: SK, IK, LK
    use pm_distNorm, only: getNormQuan
    use pm_arraySpace, only: getLinSpace
    use pm_io, only: display_type

    implicit none

    real, allocatable :: cdf(:), quantile(:), mu(:), sigma(:) ! all real kinds are supported.

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("cdf = [0., 0.05, 0.158655, 0.5, 0.658655, 0.95, 1.]")
                    cdf = [0., 0.05, 0.158655, 0.5, 0.658655, 0.95, 1.]
    call disp%show("quantile = getNormQuan(cdf)")
                    quantile = getNormQuan(cdf)
    call disp%show("quantile")
    call disp%show( quantile )
    call disp%skip()

    call disp%skip()
    call disp%show("cdf = [0., 0.05, 0.158655, 0.5, 0.658655, 0.95, 1.]")
                    cdf = [0., 0.05, 0.158655, 0.5, 0.658655, 0.95, 1.]
    call disp%show("quantile = getNormQuan(cdf, mu = 1.)")
                    quantile = getNormQuan(cdf, mu = 1.)
    call disp%show("quantile")
    call disp%show( quantile )
    call disp%skip()

    call disp%skip()
    call disp%show("cdf = [0., 0.05, 0.158655, 0.5, 0.658655, 0.95, 1.]")
                    cdf = [0., 0.05, 0.158655, 0.5, 0.658655, 0.95, 1.]
    call disp%show("quantile = getNormQuan(cdf, mu = 1., sigma = 10.)")
                    quantile = getNormQuan(cdf, mu = 1., sigma = 10.)
    call disp%show("quantile")
    call disp%show( quantile )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example quantile array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "getNormQuan.RK.txt")
        cdf = getLinSpace(0., 1., count = 1000_IK, lopen = .true._LK, fopen = .true._LK)
        do i = 1, size(cdf)
            write(fileUnit, "(5(g0,:,','))") cdf(i), getNormQuan(cdf(i), [0., 0., 0., -2.], sigma = [3., 1., 0.3, 1.])
        end do
        close(fileUnit)
    end block

end program example