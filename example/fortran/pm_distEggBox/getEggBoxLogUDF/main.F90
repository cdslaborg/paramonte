program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKC => RK ! all real kinds are supported.
    use pm_distEggBox, only: getEggBoxLogUDF
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 5_IK

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getEggBoxLogUDF([0._RKC]) ! 1D eggbox")
    call disp%show( getEggBoxLogUDF([0._RKC]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getEggBoxLogUDF([0._RKC], [0._RKC], [1._RKC], 5._RKC, 2._RKC) ! 1D eggbox")
    call disp%show( getEggBoxLogUDF([0._RKC], [0._RKC], [1._RKC], 5._RKC, 2._RKC) )
    call disp%skip()

    call disp%skip()
    call disp%show("getEggBoxLogUDF([real(RKC) :: 0, 1], mu = [-1._RKC, 1._RKC], sigma = [2._RKC, .5_RKC]) ! 2D EggBox")
    call disp%show( getEggBoxLogUDF([real(RKC) :: 0, 1], mu = [-1._RKC, 1._RKC], sigma = [2._RKC, .5_RKC]) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example density array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_mathSubAdd, only: operator(.subadd.)
        use pm_arrayMembership, only: operator(.inrange.)
        integer(IK) :: fileUnit, i, j
        real(RKC) , parameter :: signif = 2
        real(RKC) :: point(1000), density(4), mu(4), sigma(4)
        mu = [+0._RKC, +0._RKC, +0._RKC, -2._RKC]
        sigma = [+3._RKC, +1._RKC, +.3_RKC, 1._RKC]
        call setLinSpace( point &
                        , x1 = minval(mu) - signif * sigma(minloc(mu, 1)) &
                        , x2 = maxval(mu) + signif * sigma(maxloc(mu, 1)) &
                        )
        open(newunit = fileUnit, file = "getEggBoxLogUDF.D1.RK.txt")
        do i = 1, size(point)
            do j = 1, size(density)
                density(j) = getEggBoxLogUDF(point(i:i), mu(j:j), sigma(j:j))
            end do
            write(fileUnit,"(5(g0,:,','))") point(i), density
        end do
        close(fileUnit)
    end block

    block
        use pm_mathSubAdd, only: operator(.subadd.)
        use pm_arrayMembership, only: operator(.inrange.)
        integer(IK) :: fileUnit, i, j
        real(RKC)   , parameter :: signif = 2
        integer(IK) , parameter :: ndim = 2, npnt = 500
        real(RKC) :: grid(ndim, npnt, npnt), mu(ndim), sigma(ndim)
        mu = [+0._RKC, -2._RKC]
        sigma = [+3._RKC, +1._RKC]
        do i = 1, ndim
            grid(i, :, :) = spread  ( getLinSpace   ( x1 = mu(i) - signif * sigma(i) &
                                                    , x2 = mu(i) + signif * sigma(i) &
                                                    , count = npnt &
                                                    ) , 3 - i, npnt)
        end do
        open(newunit = fileUnit, file = "getEggBoxLogUDF.D2.RK.txt")
        do i = 1, size(grid, 2)
            do j = 1, size(grid, 3)
                write(fileUnit,"(5(g0,:,','))") grid(:, i, j), getEggBoxLogUDF(grid(:, i, j), mu, sigma)
            end do
        end do
        close(fileUnit)
    end block

end program example