program example

    use pm_kind, only: SK, IK, LK
    use pm_sampleMean, only: getMean
    use pm_sampleMean, only: setMean
    use pm_sampleShift, only: getShifted
    use pm_arraySpace, only: getLinSpace
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_arrayFill, only: getFilled
    use pm_sampleCov, only: uppDia
    use pm_sampleCov, only: lowDia
    use pm_sampleCov, only: getCov
    use pm_sampleCor, only: getCor
    use pm_io, only: display_type
    use pm_io, only: getFormat

    implicit none

    integer(IK) :: itry, ntry = 10
    type(display_type) :: disp
    character(:), allocatable :: format
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Convert correlation matrix and standard deviation to covariance matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported, e.g., RK32, RK64, RK128.
        use pm_matrixCopy, only: setMatCopy, rdpack
        use pm_distCov, only: getCovRand
        integer(IK) :: ndim
        real(TKC), allocatable :: cov(:,:), cor(:,:), std(:)
        format = getFormat(mold = [0._TKC], ed = SK_"es", signed = .true._LK)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("ndim = getUnifRand(1, 7)")
                            ndim = getUnifRand(1, 7)
            call disp%show("ndim")
            call disp%show( ndim )
            call disp%show("std = getUnifRand(1, 10, ndim)")
                            std = getUnifRand(1, 10, ndim)
            call disp%show("std")
            call disp%show( std , format = format )
            call disp%show("call setResized(cov, [ndim, ndim])")
                            call setResized(cov, [ndim, ndim])
            call disp%show("cor = getCovRand(1., ndim)")
                            cor = getCovRand(1., ndim)
            call disp%show("cor")
            call disp%show( cor , format = format )
            call disp%show("cov = getCov(cor, uppDia, std) ! convert upper correlation matrix to full covariance matrix.")
                            cov = getCov(cor, uppDia, std)
            call disp%show("cov")
            call disp%show( cov , format = format )
            call disp%show("cov = getCov(cor, lowDia, std) ! convert upper correlation matrix to full covariance matrix.")
                            cov = getCov(cor, lowDia, std)
            call disp%show("cov")
            call disp%show( cov , format = format )
            call disp%show("getCor(getCov(cor, lowDia, std), lowDia) ! reconstruct the original correlation matrix.")
            call disp%show( getCor(getCov(cor, lowDia, std), lowDia) , format = format )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the covariance matrix of a 2-D sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported, e.g., RK32, RK64, RK128.
        real(TKC), allocatable :: sample(:,:), cov(:,:), mean(:)
        integer(IK) :: ndim, nsam
        call disp%show("ndim = 2; nsam = 10")
                        ndim = 2; nsam = 10
        call disp%show("sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample )
        call disp%skip()
        call disp%show("Compute the sample covariance along the second dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(sample, dim = 2_IK)")
                        cov = getCov(sample, dim = 2_IK)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the sample covariance along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(transpose(sample), dim = 1_IK)")
                        cov = getCov(transpose(sample), dim = 1_IK)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the full sample covariance for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(sample(1,:), sample(2,:))")
                        cov = getCov(sample(1,:), sample(2,:))
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported, e.g., RK32, RK64, RK128.
        complex(TKC), allocatable :: sample(:,:), cov(:,:), mean(:)
        integer(IK) :: ndim, nsam
        call disp%show("ndim = 2; nsam = 10")
                        ndim = 2; nsam = 10
        call disp%show("sample = reshape(cmplx(getUnifRand(1, 20, ndim * nsam), -getUnifRand(1, 20, ndim * nsam), TKC), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(cmplx(getUnifRand(1, 20, ndim * nsam), -getUnifRand(1, 20, ndim * nsam), TKC), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample )
        call disp%skip()
        call disp%show("Compute the sample covariance along the second dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(sample, dim = 2_IK)")
                        cov = getCov(sample, dim = 2_IK)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the sample covariance along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(transpose(sample), dim = 1_IK)")
                        cov = getCov(transpose(sample), dim = 1_IK)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the full sample covariance for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(sample(1,:), sample(2,:))")
                        cov = getCov(sample(1,:), sample(2,:))
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased covariance matrix of a weighted 2-D sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported, e.g., RK32, RK64, RK128.
        real(TKC) :: rweisum
        integer(IK) :: iweisum
        real(TKC), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        real(TKC), allocatable :: sample(:,:), cov(:,:), mean(:)
        integer(IK) :: ndim, nsam
        call disp%show("ndim = 2; nsam = 10")
                        ndim = 2; nsam = 10
        call disp%show("sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("call setResized(mean, ndim)")
                        call setResized(mean, ndim)
        call disp%show("iweight = getUnifRand(1, 10, nsam) ! integer-valued weights.")
                        iweight = getUnifRand(1, 10, nsam) ! integer-valued weights.
        call disp%show("iweight")
        call disp%show( iweight )
        call disp%show("call setMean(mean, sample, 2_IK, iweight, iweisum)")
                        call setMean(mean, sample, 2_IK, iweight, iweisum)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("iweisum")
        call disp%show( iweisum )
        call disp%show("rweight = iweight ! or real-valued weights.")
                        rweight = iweight ! or real-valued weights.
        call disp%show("iweight")
        call disp%show( iweight )
        call disp%show("call setMean(mean, sample, 2_IK, rweight, rweisum)")
                        call setMean(mean, sample, 2_IK, rweight, rweisum)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("rweisum")
        call disp%show( rweisum )

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Compute the covariance matrix with integer weights.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%show("cov = getCov(sample, 2_IK, iweight)")
                        cov = getCov(sample, 2_IK, iweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the sample covariance along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(transpose(sample), 1_IK, iweight)")
                        cov = getCov(transpose(sample), 1_IK, iweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the full sample covariance for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(sample(1,:), sample(2,:), weight = iweight)")
                        cov = getCov(sample(1,:), sample(2,:), weight = iweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Compute the covariance matrix with    real weights.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%show("cov = getCov(sample, 2_IK, rweight)")
                        cov = getCov(sample, 2_IK, rweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the sample covariance along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(transpose(sample), 1_IK, rweight)")
                        cov = getCov(transpose(sample), 1_IK, rweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the full sample covariance for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(sample(1,:), sample(2,:), weight = rweight)")
                        cov = getCov(sample(1,:), sample(2,:), weight = rweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported, e.g., RK32, RK64, RK128.
        real(TKC) :: rweisum
        integer(IK) :: iweisum
        real(TKC), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        complex(TKC), allocatable :: sample(:,:), cov(:,:), mean(:)
        integer(IK) :: ndim, nsam
        call disp%show("ndim = 2; nsam = 10")
                        ndim = 2; nsam = 10
        call disp%show("sample = reshape(cmplx(getUnifRand(1, 20, ndim * nsam), -getUnifRand(1, 20, ndim * nsam), TKC), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(cmplx(getUnifRand(1, 20, ndim * nsam), -getUnifRand(1, 20, ndim * nsam), TKC), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("call setResized(mean, ndim)")
                        call setResized(mean, ndim)
        call disp%show("iweight = getUnifRand(1, 10, nsam) ! integer-valued weights.")
                        iweight = getUnifRand(1, 10, nsam) ! integer-valued weights.
        call disp%show("iweight")
        call disp%show( iweight )
        call disp%show("call setMean(mean, sample, 2_IK, iweight, iweisum)")
                        call setMean(mean, sample, 2_IK, iweight, iweisum)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("iweisum")
        call disp%show( iweisum )
        call disp%show("rweight = iweight ! or real-valued weights.")
                        rweight = iweight ! or real-valued weights.
        call disp%show("iweight")
        call disp%show( iweight )
        call disp%show("call setMean(mean, sample, 2_IK, rweight, rweisum)")
                        call setMean(mean, sample, 2_IK, rweight, rweisum)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("rweisum")
        call disp%show( rweisum )

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Compute the covariance matrix with integer weights.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%show("cov = getCov(sample, 2_IK, iweight)")
                        cov = getCov(sample, 2_IK, iweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the sample covariance along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(transpose(sample), 1_IK, iweight)")
                        cov = getCov(transpose(sample), 1_IK, iweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the full sample covariance for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(sample(1,:), sample(2,:), weight = iweight)")
                        cov = getCov(sample(1,:), sample(2,:), weight = iweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Compute the covariance matrix with    real weights.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%show("cov = getCov(sample, 2_IK, rweight)")
                        cov = getCov(sample, 2_IK, rweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the sample covariance along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(transpose(sample), 1_IK, rweight)")
                        cov = getCov(transpose(sample), 1_IK, rweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
        call disp%show("Compute the full sample covariance for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cov = getCov(sample(1,:), sample(2,:), weight = rweight)")
                        cov = getCov(sample(1,:), sample(2,:), weight = rweight)
        call disp%show("cov")
        call disp%show( cov )
        call disp%skip()
    end block

end program example