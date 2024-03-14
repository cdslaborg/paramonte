program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK, RK ! all real kinds are supported.
    use pm_io, only: display_type
    use pm_matrixChol, only: getMatChol, uppDia, lowDia
    use pm_distUnifEll, only: getUnifEllRand

    implicit none

    integer(IK), parameter  :: NDIM = 2_IK  ! The number of dimensions of the MVN distribution.
    real(RK)                :: mean(NDIM), gramian(NDIM, NDIM), choLow(NDIM, NDIM), rand(NDIM)
    integer(IK)             :: info

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random numbers from the (Standard) Multivariate Uniform Ellipsoidal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Multivariate Uniform Ellipsoidal random vector with a particular mean and Identity Gramian matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean = [-5., -5.]")
                    mean = [-5., -5.]
    call disp%show("rand = getUnifEllRand(mean)")
                    rand = getUnifEllRand(mean)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Multivariate Uniform Ellipsoidal random vector with zero mean and Gramian matrix specified via the Cholesky Lower Triangle.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("gramian = reshape([1., 1., 1., 4.], shape = [NDIM, NDIM])")
                    gramian = reshape([1., 1., 1., 4.], shape = [NDIM, NDIM])
    call disp%show("choLow = getMatChol(gramian, lowDia)")
                    choLow = getMatChol(gramian, lowDia)
    call disp%show("choLow")
    call disp%show( choLow )
    call disp%show("rand = getUnifEllRand(choLow, lowDia)")
                    rand = getUnifEllRand(choLow, lowDia)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Multivariate Uniform Ellipsoidal random vector with given mean and Gramian matrix specified via the Cholesky Lower Triangle.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean")
    call disp%show( mean )
    call disp%show("choLow")
    call disp%show( choLow )
    call disp%show("rand = getUnifEllRand(mean, choLow, lowDia)")
                    rand = getUnifEllRand(mean, choLow, lowDia)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        integer(IK), parameter :: NVEC = 5000_IK
        open(newunit = fileUnit, file = "getUnifEllRandMean.RK.txt")
            do i = 1, NVEC
                write(fileUnit,"(*(g0,:,','))") getUnifEllRand(mean)
            end do
        close(fileUnit)
        open(newunit = fileUnit, file = "getUnifEllRandChol.RK.txt")
            do i = 1, NVEC
                write(fileUnit,"(*(g0,:,','))") getUnifEllRand(choLow, lowDia)
            end do
        close(fileUnit)
        open(newunit = fileUnit, file = "getUnifEllRandMeanChol.RK.txt")
            do i = 1, NVEC
                write(fileUnit,"(*(g0,:,','))") getUnifEllRand(mean, choLow, lowDia)
            end do
        close(fileUnit)
    end block

end program example