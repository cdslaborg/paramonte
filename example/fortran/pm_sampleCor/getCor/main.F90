program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_distCov, only: getCovRand
    use pm_arrayRange, only: getRange
    use pm_sampleCor, only: getCor
    use pm_sampleCor, only: uppDia
    use pm_sampleCor, only: lowDia
    use pm_sampleCor, only: upp
    use pm_sampleCor, only: low
    use pm_sampleCov, only: getCov
    use pm_sampleMean, only: getMean
    use pm_sampleMean, only: setMean
    use pm_sampleShift, only: getShifted
    use pm_arraySpace, only: getLinSpace
    use pm_arrayResize, only: setResized
    use pm_arrayFill, only: getFilled
    use pm_io, only: getFormat

    implicit none

    type(display_type) :: disp
    integer(IK) :: itry, ntry = 10
    character(:), allocatable :: format
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the correlation matrix from covariance matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS ! All other real types are also supported.
        real(TKG), allocatable :: cov(:,:), cor(:,:), std(:)
        integer(IK) :: ndim
        do itry = 1, 10
            call disp%skip()
            call disp%show("ndim = getUnifRand(0, 7)")
                            ndim = getUnifRand(0, 7)
            call disp%show("ndim")
            call disp%show( ndim )
            call disp%show("std = getUnifRand(1, ndim, ndim)")
                            std = getUnifRand(1, ndim, ndim)
            call disp%show("std")
            call disp%show( std )
            call disp%show("cov = getCovRand(1._TKG, std)")
                            cov = getCovRand(1._TKG, std)
            call disp%show("cov")
            call disp%show( cov )
            call disp%show("cor = getCor(cov, uppDia)")
                            cor = getCor(cov, uppDia)
            call disp%show("cor")
            call disp%show( cor )
            call disp%show("cor = getCor(cov, upp, stdinv = 1 / std)")
                            cor = getCor(cov, upp, stdinv = 1 / std)
            call disp%show("cor")
            call disp%show( cor )
            call disp%show("cor = getCor(cov, lowDia)")
                            cor = getCor(cov, lowDia)
            call disp%show("cor")
            call disp%show( cor )
            call disp%show("cor = getCor(cov, low, stdinv = 1 / std)")
                            cor = getCor(cov, low, stdinv = 1 / std)
            call disp%show("cor")
            call disp%show( cor )
            call disp%show("getCov(getCor(cov, lowDia), lowDia, std) ! reconstruct the original covariance matrix.")
            call disp%show( getCov(getCor(cov, lowDia), lowDia, std) , format = format )
            call disp%show("getCov(getCor(cov, uppDia), uppDia, std) ! reconstruct the original covariance matrix.")
            call disp%show( getCov(getCor(cov, uppDia), uppDia, std) , format = format )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Pearson correlation matrix for a pair of time series.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS ! All other real types are also supported.
        integer(IK) :: ndim, nsam, dim
        real(TKG), allocatable :: sample(:,:), cor(:,:), mean(:)
        format = getFormat(mold = [0._TKG], ed = SK_"es", signed = .true._LK)
        call disp%show("ndim = 2; nsam = 10; dim = 2")
                        ndim = 2; nsam = 10; dim = 2
        call disp%show("sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample , format = format )
        call disp%show("mean = getMean(sample, dim)")
                        mean = getMean(sample, dim)
        call disp%show("mean")
        call disp%show( mean , format = format )
        call disp%show("cor = getCor(sample, dim) ! same result as above.")
                        cor = getCor(sample, dim)
        call disp%show("cor")
        call disp%show( cor , format = format )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("dim = 1")
                        dim = 1
        call disp%show("cor = getCor(transpose(sample), dim) ! same result as above.")
                        cor = getCor(transpose(sample), dim)
        call disp%show("cor")
        call disp%show( cor , format = format )
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cor(1,1) = getCor(sample(1,:), sample(2,:)) ! same result as above.")
                        cor(1,1) = getCor(sample(1,:), sample(2,:))
        call disp%show("cor(1,1)")
        call disp%show( cor(1,1) , format = format )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Pearson correlation matrix for a weighted pair of time series.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_arrayVerbose, only: getVerbose
        use pm_kind, only: TKG => RKS ! All other real types are also supported.
        integer(IK), allocatable :: iweight(:)
        real(TKG), allocatable :: rweight(:)
        real(TKG) :: rweisum
        integer(IK) :: iweisum
        integer(IK) :: ndim, nsam, dim
        real(TKG), allocatable :: sample(:,:), cor(:,:), mean(:)
        format = getFormat(mold = [0._TKG], ed = SK_"es", signed = .true._LK)
        call disp%show("ndim = 2; nsam = 10; dim = 2")
                        ndim = 2; nsam = 10; dim = 2
        call disp%show("sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample , format = format )
        call disp%show("call setResized(mean, ndim)")
                        call setResized(mean, ndim)
        call disp%show("iweight = getUnifRand(1, 10, nsam) ! integer-valued weights.")
                        iweight = getUnifRand(1, 10, nsam) ! integer-valued weights.
        call disp%show("iweight")
        call disp%show( iweight )
        call disp%show("call setMean(mean, sample, dim, iweight, iweisum)")
                        call setMean(mean, sample, dim, iweight, iweisum)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("iweisum")
        call disp%show( iweisum )
        call disp%show("rweight = iweight ! or real-valued weights.")
                        rweight = iweight ! or real-valued weights.
        call disp%show("iweight")
        call disp%show( iweight )
        call disp%show("call setMean(mean, sample, dim, rweight, rweisum)")
                        call setMean(mean, sample, dim, rweight, rweisum)
        call disp%show("mean")
        call disp%show( mean , format = format )
        call disp%show("rweisum")
        call disp%show( rweisum , format = format )

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Compute the correlation matrix with integer weights.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%show("cor = getCor(sample, dim, iweight) ! same result as above.")
                        cor = getCor(sample, dim, iweight)
        call disp%show("cor")
        call disp%show( cor , format = format )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("dim = 1")
                        dim = 1
        call disp%show("cor = getCor(transpose(sample), dim, iweight) ! same result as above.")
                        cor = getCor(transpose(sample), dim, iweight)
        call disp%show("cor")
        call disp%show( cor , format = format )
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cor(1,1) = getCor(sample(1,:), sample(2,:), iweight) ! same result as above.")
                        cor(1,1) = getCor(sample(1,:), sample(2,:), iweight)
        call disp%show("cor(1,1)")
        call disp%show( cor(1,1) , format = format )
        call disp%skip()

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Compute the correlation matrix with    real weights.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%show("dim = 2")
                        dim = 2
        call disp%show("cor = getCor(sample, dim, rweight) ! same result as above.")
                        cor = getCor(sample, dim, rweight)
        call disp%show("cor")
        call disp%show( cor , format = format )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("dim = 1")
                        dim = 1
        call disp%show("cor = getCor(transpose(sample), dim, rweight) ! same result as above.")
                        cor = getCor(transpose(sample), dim, rweight)
        call disp%show("cor")
        call disp%show( cor , format = format )
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("cor(1,1) = getCor(sample(1,:), sample(2,:), rweight) ! same result as above.")
                        cor(1,1) = getCor(sample(1,:), sample(2,:), rweight)
        call disp%show("cor(1,1)")
        call disp%show( cor(1,1) , format = format )
        call disp%skip()
    end block

end program example