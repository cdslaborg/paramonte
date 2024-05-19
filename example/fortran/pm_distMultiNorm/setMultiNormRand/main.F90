program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK, RKG => RKS ! all real kinds are supported.
    use pm_io, only: display_type
    use pm_arrayResize, only: setResized
    use pm_matrixChol, only: getMatChol
    use pm_distMultiNorm, only: setMultiNormRand, uppDia, lowDia

    implicit none

    real(RKG), allocatable :: mean(:), cov(:,:), chol(:,:), rand(:)
    integer(IK) :: info

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
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Multivariate Normal random vector from a Standard distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setResized(rand, 3_IK)")
                    call setResized(rand, 3_IK)
    call disp%show("call setMultiNormRand(rand)")
                    call setMultiNormRand(rand)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Multivariate Normal random vector with a particular mean and Identity covariance matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean = [-5., -5.]")
                    mean = [-5., -5.]
    call disp%show("mean")
    call disp%show( mean )
    call disp%show("call setResized(rand, size(mean, 1, IK))")
                    call setResized(rand, size(mean, 1, IK))
    call disp%show("call setMultiNormRand(rand, mean)")
                    call setMultiNormRand(rand, mean)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Multivariate Normal random vector with zero mean and Covariance matrix specified via the Cholesky Lower Triangle.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("cov = reshape([1., 1., 1., 4.], shape = [2, 2])")
                    cov = reshape([1., 1., 1., 4.], shape = [2, 2])
    call disp%show("cov")
    call disp%show( cov )
    call disp%show("chol = getMatChol(cov, uppDia)")
                    chol = getMatChol(cov, uppDia)
    call disp%show("call setResized(rand, size(chol, 1, IK))")
                    call setResized(rand, size(chol, 1, IK))
    call disp%show("call setMultiNormRand(rand, chol, uppDia)")
                    call setMultiNormRand(rand, chol, uppDia)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Multivariate Normal random vector with given mean and Covariance matrix specified via the Cholesky Lower Triangle.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean = [-5., -5.]")
                    mean = [-5., -5.]
    call disp%show("mean")
    call disp%show( mean )
    call disp%show("cov = reshape([1., 1., 1., 4.], shape = [2, 2])")
                    cov = reshape([1., 1., 1., 4.], shape = [2, 2])
    call disp%show("cov")
    call disp%show( cov )
    call disp%show("chol = getMatChol(cov, uppDia)")
                    chol = getMatChol(cov, uppDia)
    call disp%show("call setResized(rand, size(chol, 1, IK))")
                    call setResized(rand, size(chol, 1, IK))
    call disp%show("call setMultiNormRand(rand, mean, chol, uppDia)")
                    call setMultiNormRand(rand, mean, chol, uppDia)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_io, only: getErrTableWrite, trans
        real(RKG) :: rand(2, 5000)

        call setMultiNormRand(rand)
        if (0 /= getErrTableWrite("setMultiNormRand.RK.txt", rand, trans)) error stop 'table write failed.'

        call setMultiNormRand(rand, mean)
        if (0 /= getErrTableWrite("setMultiNormRandMean.RK.txt", rand, trans)) error stop 'table write failed.'

        call setMultiNormRand(rand, chol, uppDia)
        if (0 /= getErrTableWrite("setMultiNormRandChol.RK.txt", rand, trans)) error stop 'table write failed.'

        call setMultiNormRand(rand, chol, uppDia)
        if (0 /= getErrTableWrite("setMultiNormRandMeanChol.RK.txt", rand, trans)) error stop 'table write failed.'

    end block

end program example