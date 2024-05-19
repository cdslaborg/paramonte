program example

    use pm_kind, only: SK, IK
    use pm_kind, only: TKG => RKS ! All other real types are also supported.
    use pm_sampleCov, only: getCov
    use pm_sampleMean, only: getMean
    use pm_sampleCov, only: uppDia, lowDia
    use pm_sampleCov, only: setCovUpdated
    use pm_arrayRebind, only: setRebound
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_arrayRange, only: getRange
    use pm_io, only: display_type

    implicit none

    integer(IK) :: isam, ndim
    integer(IK) , allocatable :: iweight(:)
    integer(IK) :: dim, idim, nsamA, nsamB, itry, ntry = 10
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased merged covariance of a multivariate sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG), allocatable :: cov(:,:), covA(:,:)
        real(TKG), allocatable :: mean(:), meanA(:), meanB(:)
        real(TKG), allocatable :: sample(:,:), sampleA(:,:), sampleB(:,:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("dim = 2; ndim = getUnifRand(1, 5); nsamA = getUnifRand(ndim + 1_IK, ndim * 2_IK); nsamB = getUnifRand(1_IK, nsamA)")
                            dim = 2; ndim = getUnifRand(1, 5); nsamA = getUnifRand(ndim + 1_IK, ndim * 2_IK); nsamB = getUnifRand(1_IK, nsamA)
            call disp%show("[ndim, nsamA]")
            call disp%show( [ndim, nsamA] )
            call disp%show("sampleA = getUnifRand(-1., +1., ndim, nsamA) ! Generate a non-singular sample.")
                            sampleA = getUnifRand(-1., +1., ndim, nsamA)
            call disp%show("sampleA")
            call disp%show( sampleA )
            call disp%show("covA = getCov(sampleA, dim)")
                            covA = getCov(sampleA, dim)
            call disp%show("covA")
            call disp%show( covA )
            call disp%show("meanA = getMean(sampleA, dim)")
                            meanA = getMean(sampleA, dim)
            call disp%show("meanA")
            call disp%show( meanA )
            call disp%show("sampleB = spread([(idim, idim = 1, ndim)], dim, nsamB) ! Generate a singular sample.")
                            sampleB = spread([(idim, idim = 1, ndim)], dim, nsamB)
            call disp%show("sampleB")
            call disp%show( sampleB )
            call disp%show("meanB = getMean(sampleB, dim)")
                            meanB = getMean(sampleB, dim)
            call disp%show("meanB")
            call disp%show( meanB )
            call disp%show("sample = reshape([sampleA, sampleB], [ndim, nsamA + nsamB])")
                            sample = reshape([sampleA, sampleB], [ndim, nsamA + nsamB])
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("cov = getCov(sample, dim)")
                            cov = getCov(sample, dim)
            call disp%show("cov")
            call disp%show( cov )
            call disp%show("mean = getMean(sample, dim)")
                            mean = getMean(sample, dim)
            call disp%show("mean")
            call disp%show( mean )
            call disp%show("call setCovUpdated(covA, meanA - meanB, real(nsamA, TKG) / real(nsamA + nsamB, TKG), uppDia)")
                            call setCovUpdated(covA, meanA - meanB, real(nsamA, TKG) / real(nsamA + nsamB, TKG), uppDia)
            call disp%show("covA")
            call disp%show( covA )
            call disp%show("cov ! reference")
            call disp%show( cov )
            call disp%skip()
        end do
    end block

end program example