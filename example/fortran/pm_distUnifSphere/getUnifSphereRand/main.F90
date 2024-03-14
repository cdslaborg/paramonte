program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK, RK ! all real kinds are supported.
    use pm_io, only: display_type
    use pm_matrixChol, only: getMatChol, uppDia, lowDia
    use pm_distUnifSphere, only: getUnifSphereRand
    use pm_distCov, only: getCovRand

    implicit none

    integer(IK), parameter  :: NDIM = 2_IK  ! The number of dimensions of the MVN distribution.
    real(RK)                :: mean(NDIM), gramian(NDIM, NDIM), chol(NDIM, NDIM), rand(NDIM)
    integer(IK)             :: info

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random numbers on the surface an n-sphere or its affine transformation embedded in ndim-dimensions.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean = [-5., -5.]")
                    mean = [-5., -5.]
    call disp%show("rand = getUnifSphereRand(mean)")
                    rand = getUnifSphereRand(mean)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    call disp%skip()
    call disp%show("gramian = getCovRand(mold = 1., ndim = ndim)")
                    gramian = getCovRand(mold = 1., ndim = ndim)
    call disp%show("chol = getMatChol(gramian, lowDia)")
                    chol = getMatChol(gramian, lowDia)
    call disp%show("chol")
    call disp%show( chol )
    call disp%show("rand = getUnifSphereRand(chol, lowDia)")
                    rand = getUnifSphereRand(chol, lowDia)
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
    call disp%show("chol")
    call disp%show( chol )
    call disp%show("rand = getUnifSphereRand(mean, chol, lowDia)")
                    rand = getUnifSphereRand(mean, chol, lowDia)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        integer(IK), parameter :: NVEC = 100_IK
        open(newunit = fileUnit, file = "getUnifSphereRandMean.RK.txt")
            do i = 1, NVEC
                write(fileUnit,"(*(g0,:,','))") getUnifSphereRand(mean)
            end do
        close(fileUnit)
        open(newunit = fileUnit, file = "getUnifSphereRandChol.RK.txt")
            chol = getMatChol(getCovRand(mold = 1., ndim = ndim), lowDia)
            do i = 1, NVEC
                write(fileUnit,"(*(g0,:,','))") getUnifSphereRand(chol, lowDia)
            end do
        close(fileUnit)
        open(newunit = fileUnit, file = "getUnifSphereRandMeanChol.RK.txt")
            chol = getMatChol(getCovRand(ndim = ndim, mold = 1.), lowDia)
            do i = 1, NVEC
                write(fileUnit,"(*(g0,:,','))") getUnifSphereRand(mean, chol, lowDia)
            end do
        close(fileUnit)
    end block

end program example