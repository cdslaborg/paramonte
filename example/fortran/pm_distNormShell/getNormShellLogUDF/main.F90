program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RK ! all real kinds are supported.
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
    call disp%show("getNormShellLogUDF([0._RKG]) ! 1D")
    call disp%show( getNormShellLogUDF([0._RKG]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormShellLogUDF([0._RKG], [0._RKG], reshape([1._RKG], [1,1]), 5._RKG, 2._RKG) ! 1D")
    call disp%show( getNormShellLogUDF([0._RKG], [0._RKG], reshape([1._RKG], [1,1]), 5._RKG, 2._RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormShellLogUDF([real(RKG) :: 0, 1], center = [-1._RKG, 1._RKG], invCov = reshape([1._RKG, .5_RKG, .5_RKG, 1._RKG], [2, 2]), width = 0.5_RKG, radius = 2._RKG) ! 2D")
    call disp%show( getNormShellLogUDF([real(RKG) :: 0, 1], center = [-1._RKG, 1._RKG], invCov = reshape([1._RKG, .5_RKG, .5_RKG, 1._RKG], [2, 2]), width = 0.5_RKG, radius = 2._RKG) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example density array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i, j
        integer(IK), parameter :: nmix = 4
        real(RKG) , parameter :: signif = 2
        real(RKG) :: point(1000), width(nmix), radius(nmix), center(1, nmix), invCov(1, 1, nmix)
        center = reshape([-1._RKG, +0._RKG, +1._RKG, +2._RKG], shape(center))
        invCov = 1._RKG / reshape([+1._RKG, +2._RKG, +5._RKG, 0.5_RKG], shape(invCov))
        radius = [1._RKG, 5._RKG, 2._RKG, 0.5_RKG]
        width = .5_RKG
        call setLinSpace( point &
                        , x1 = -10._RKG &
                        , x2 = +10._RKG &
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
        real(RKG)   , parameter :: signif = 2_IK
        integer(IK) , parameter :: ndim = 2_IK, nmix = 2_IK, npnt = 700_IK
        real(RKG) :: grid(ndim, npnt, npnt), center(ndim, nmix), invCov(ndim, ndim, nmix), width(nmix), radius(nmix)
        center = 2 * reshape([-1._RKG, -1._RKG, 1._RKG, 1._RKG], shape(center))
        invCov = spread(getMatInit([ndim, ndim], uppLowDia, 0._RKG, 0._RKG, 1._RKG), 3, nmix)
        invCov(1,2,1) = 0.5
        invCov(2,1,1) = 0.5
        radius = [2._RKG, 2._RKG]
        width = .5_RKG * [1._RKG, 1._RKG]
        do i = 1, ndim
            grid(i, :, :) = spread(getLinSpace(x1 = -6._RKG, x2 = +6._RKG, count = npnt) , 3 - i, npnt)
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