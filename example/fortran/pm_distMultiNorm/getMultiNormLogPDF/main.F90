program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKC => RK ! all real kinds are supported.
    use pm_distMultiNorm, only: getMultiNormLogPDF
    use pm_arraySpace, only: getLinSpace
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 5_IK

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getMultiNormLogPDF([0._RKC]) ! 1D Norm")
    call disp%show( getMultiNormLogPDF([0._RKC]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getMultiNormLogPDF([0._RKC], mean = [10._RKC], invCov = reshape([1._RKC], [1, 1])) ! 1D Norm")
    call disp%show( getMultiNormLogPDF([0._RKC], mean = [10._RKC], invCov = reshape([1._RKC], [1, 1])) )
    call disp%skip()

    call disp%skip()
    call disp%show("getMultiNormLogPDF([real(RKC) :: 0, 1], mean = [-1._RKC, 1._RKC], invCov = reshape([+1.3_RKC, -.66_RKC, -.66_RKC, 1.3_RKC], [2, 2])) ! 2D Norm")
    call disp%show( getMultiNormLogPDF([real(RKC) :: 0, 1], mean = [-1._RKC, 1._RKC], invCov = reshape([+1.3_RKC, -.66_RKC, -.66_RKC, 1.3_RKC], [2, 2])) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example density array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i, j
        real(RKC)   , parameter :: signif = 5
        integer(IK) , parameter :: ndim = 2, npnt = 500
        real(RKC) :: grid(ndim, npnt, npnt), mean(ndim), invCov(ndim, ndim)
        mean = [+5._RKC, -5._RKC]
        invCov = reshape([+1.3_RKC, -.66_RKC, -.66_RKC, 1.3_RKC], [ndim, ndim])
        do i = 1, ndim
            grid(i, :, :) = spread  ( getLinSpace   ( x1 = mean(i) - signif &
                                                    , x2 = mean(i) + signif &
                                                    , count = npnt &
                                                    ) , 3 - i, npnt)
        end do
        open(newunit = fileUnit, file = "getMultiNormLogPDF.D2.RK.txt")
        do i = 1, size(grid, 2)
            do j = 1, size(grid, 3)
                write(fileUnit,"(5(g0,:,','))") grid(:, i, j), exp(getMultiNormLogPDF(grid(:, i, j), mean, invCov))
            end do
        end do
        close(fileUnit)
    end block

end program example