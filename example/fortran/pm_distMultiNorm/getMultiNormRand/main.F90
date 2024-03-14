program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK, RKC => RK32 ! all real kinds are supported.
    use pm_io, only: display_type
    use pm_matrixChol, only: getMatChol, uppDia
    use pm_distMultiNorm, only: getMultiNormRand
    use pm_distCov, only: getCovRand

    implicit none

    real(RKC), allocatable :: mean(:), chol(:,:), rand(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random numbers from the (Standard) Multivariate Normal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Multivariate Normal random vector with a particular mean and Identity covariance matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean = [-5., -5.]")
                    mean = [-5., -5.]
    call disp%show("rand = getMultiNormRand(mean)")
                    rand = getMultiNormRand(mean)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Multivariate Normal random vector with zero mean and Covariance matrix specified via the Cholesky Lower Triangle.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("chol = getMatChol(getCovRand(mold = 1._RKC, scale = [real(RKC) :: 1, 2]), uppDia)")
                    chol = getMatChol(getCovRand(mold = 1._RKC, scale = [real(RKC) :: 1, 2]), uppDia)
    call disp%show("chol")
    call disp%show( chol )
    call disp%show("rand = getMultiNormRand(chol, uppDia)")
                    rand = getMultiNormRand(chol, uppDia)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Multivariate Normal random vector with given mean and Covariance matrix specified via the Cholesky Lower Triangle.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean")
    call disp%show( mean )
    call disp%show("chol")
    call disp%show( chol )
    call disp%show("rand = getMultiNormRand(mean, chol, uppDia)")
                    rand = getMultiNormRand(mean, chol, uppDia)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_io, only: getErrTableWrite, trans
        real(RKC) :: rand(2, 5000)

        rand = getMultiNormRand(mean, size(rand, 2, IK))
        if (0 /= getErrTableWrite("getMultiNormRandMean.RK.txt", rand, trans)) error stop 'table write failed.'

        rand = getMultiNormRand(chol, uppDia, size(rand, 2, IK))
        if (0 /= getErrTableWrite("getMultiNormRandChol.RK.txt", rand, trans)) error stop 'table write failed.'

        rand = getMultiNormRand(chol, uppDia, size(rand, 2, IK))
        if (0 /= getErrTableWrite("getMultiNormRandMeanChol.RK.txt", rand, trans)) error stop 'table write failed.'

    end block

end program example