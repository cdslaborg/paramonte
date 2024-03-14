program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKC => RK ! all real kinds are supported.
    use pm_distNormShell, only: getNormShellLogUDF
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 5_IK

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getNormShellLogUDF([0._RKC]) ! 1D")
    call disp%show( getNormShellLogUDF([0._RKC]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormShellLogUDF([0._RKC], [0._RKC], reshape([1._RKC], [1,1]), 5._RKC, 2._RKC) ! 1D")
    call disp%show( getNormShellLogUDF([0._RKC], [0._RKC], reshape([1._RKC], [1,1]), 5._RKC, 2._RKC) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormShellLogUDF([real(RKC) :: 0, 1], center = [-1._RKC, 1._RKC], invCov = reshape([1._RKC, .5_RKC, .5_RKC, 1._RKC], [2, 2]), width = 0.5_RKC, radius = 2._RKC) ! 2D")
    call disp%show( getNormShellLogUDF([real(RKC) :: 0, 1], center = [-1._RKC, 1._RKC], invCov = reshape([1._RKC, .5_RKC, .5_RKC, 1._RKC], [2, 2]), width = 0.5_RKC, radius = 2._RKC) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example density array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i, j
        integer(IK), parameter :: nmix = 4
        real(RKC) , parameter :: signif = 2
        real(RKC) :: point(1000), width(nmix), radius(nmix), center(1, nmix), invCov(1, 1, nmix)
        center = reshape([-1._RKC, +0._RKC, +1._RKC, +2._RKC], shape(center))
        invCov = 1._RKC / reshape([+1._RKC, +2._RKC, +5._RKC, 0.5_RKC], shape(invCov))
        radius = [1._RKC, 5._RKC, 2._RKC, 0.5_RKC]
        width = .5_RKC
        call setLinSpace( point &
                        , x1 = -10._RKC &
                        , x2 = +10._RKC &
                        )
        open(newunit = fileUnit, file = "getNormShellLogUDF.D1.RK.txt")
        do i = 1, size(point)
            write(fileUnit,"(*(f0.8,:,','))") point(i), exp(getNormShellLogUDF(point(i:i), center, invCov, width, radius))
        end do
        close(fileUnit)
    end block

    block
        use pm_mathSubAdd, only: operator(.subadd.)
        use pm_arrayMembership, only: operator(.inrange.)
        use pm_matrixInit, only: uppLowDia, getMatInit
        use pm_mathLogSumExp, only: getLogSumExp
        integer(IK) :: fileUnit, i, j
        real(RKC)   , parameter :: signif = 2_IK
        integer(IK) , parameter :: ndim = 2_IK, nmix = 2_IK, npnt = 700_IK
        real(RKC) :: grid(ndim, npnt, npnt), center(ndim, nmix), invCov(ndim, ndim, nmix), width(nmix), radius(nmix)
        center = 2 * reshape([-1._RKC, -1._RKC, 1._RKC, 1._RKC], shape(center))
        invCov = spread(getMatInit([ndim, ndim], uppLowDia, 0._RKC, 0._RKC, 1._RKC), 3, nmix)
        invCov(1,2,1) = 0.5
        invCov(2,1,1) = 0.5
        radius = [2._RKC, 2._RKC]
        width = .5_RKC * [1._RKC, 1._RKC]
        do i = 1, ndim
            grid(i, :, :) = spread(getLinSpace(x1 = -6._RKC, x2 = +6._RKC, count = npnt) , 3 - i, npnt)
        end do
        open(newunit = fileUnit, file = "getNormShellLogUDF.D2.RK.txt")
        do i = 1, size(grid, 2)
            do j = 1, size(grid, 3)
                write(fileUnit,"(*(f0.8,:,','))") grid(:, i, j), exp(getLogSumExp(getNormShellLogUDF(grid(:, i, j), center, invCov, width, radius)))
            end do
        end do
        close(fileUnit)
    end block

end program example