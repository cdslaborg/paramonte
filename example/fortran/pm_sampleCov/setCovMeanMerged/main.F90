program example

    use pm_kind, only: SK, IK
    use pm_kind, only: TKC => RKS ! All other real types are also supported.
    use pm_sampleCov, only: getCov
    use pm_sampleMean, only: getMean
    use pm_sampleCov, only: uppDia, lowDia
    use pm_sampleCov, only: setCovMeanMerged
    use pm_arrayRebind, only: setRebound
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_arrayRange, only: getRange
    use pm_io, only: display_type

    implicit none

    integer(IK) , parameter :: nsam = 2
    integer(IK) , allocatable :: iweight(:)
    integer(IK) :: isam, ndim, lb(nsam), ub(nsam)
    integer(IK) :: dim, itry, ntry = 10
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased merged covariance of a multivariate sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKC), allocatable :: mean(:,:), cov(:,:,:), meanMerged(:), covMerged(:,:)
        real(TKC), allocatable :: sample(:,:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("dim = 2; lb(1) = 1; ub(1) = getUnifRand(2, 7)")
                            dim = 2; lb(1) = 1; ub(1) = getUnifRand(2, 7)
            call disp%show("do isam = 2, nsam")
            call disp%show("    lb(isam) = ub(isam - 1) + 1")
            call disp%show("    ub(isam) = ub(isam - 1) + getUnifRand(2, 7)")
            call disp%show("end do")
                            do isam = 2, nsam
                                lb(isam) = ub(isam - 1) + 1
                                ub(isam) = ub(isam - 1) + getUnifRand(2, 7)
                            end do
            call disp%show("lb")
            call disp%show( lb )
            call disp%show("ub")
            call disp%show( ub )
            call disp%show("ndim = getUnifRand(1, minval(ub - lb + 1, 1))")
                            ndim = getUnifRand(1, minval(ub - lb + 1, 1))
            call disp%show("call setRebound(cov, [1_IK, 1_IK, 0_IK], [ndim, ndim, nsam])")
                            call setRebound(cov, [1_IK, 1_IK, 0_IK], [ndim, ndim, nsam])
            call disp%show("call setRebound(mean, [1_IK, 0_IK], [ndim, nsam])")
                            call setRebound(mean, [1_IK, 0_IK], [ndim, nsam])
            call disp%show("call setResized(covMerged, [ndim, ndim])")
                            call setResized(covMerged, [ndim, ndim])
            call disp%show("call setResized(meanMerged, ndim)")
                            call setResized(meanMerged, ndim)
            call disp%show("sample = getUnifRand(-1., +1., ndim, ub(nsam))")
                            sample = getUnifRand(-1., +1., ndim, ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("cov(:,:,0) = getCov(sample, dim)")
                            cov(:,:,0) = getCov(sample, dim)
            call disp%show("cov(:,:,0) ! reference")
            call disp%show( cov(:,:,0) )
            call disp%show("mean(:,0) = getMean(sample, dim)")
                            mean(:,0) = getMean(sample, dim)
            call disp%show("mean(:,0) ! reference")
            call disp%show( mean(:,0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), dim)")
            call disp%show("    mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim)")
            call disp%show("end do")
                            do isam = 1, nsam
                                cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), dim)
                                mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim)
                            end do
            call disp%show("call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(ub(1), TKC) / real(ub(2), TKC), uppDia)")
                            call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(ub(1), TKC) / real(ub(2), TKC), uppDia)
            call disp%show("covMerged")
            call disp%show( covMerged )
            call disp%show("call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(ub(1), TKC) / real(ub(2), TKC), lowDia)")
                            call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(ub(1), TKC) / real(ub(2), TKC), lowDia)
            call disp%show("covMerged")
            call disp%show( covMerged )
            call disp%show("call setCovMeanMerged(cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(ub(1), TKC) / real(ub(2), TKC), uppDia)")
                            call setCovMeanMerged(cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(ub(1), TKC) / real(ub(2), TKC), uppDia)
            call disp%show("cov(:,:,2)")
            call disp%show( cov(:,:,2) )
            call disp%show("cov(:,:,0) ! reference")
            call disp%show( cov(:,:,0) )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased merged covariance of a frequency weighted multivariate sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKC), allocatable :: mean(:,:), cov(:,:,:), meanMerged(:), covMerged(:,:)
        real(TKC), allocatable :: sample(:,:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("dim = 2; lb(1) = 1; ub(1) = getUnifRand(2, 7)")
                            dim = 2; lb(1) = 1; ub(1) = getUnifRand(2, 7)
            call disp%show("do isam = 2, nsam")
            call disp%show("    lb(isam) = ub(isam - 1) + 1")
            call disp%show("    ub(isam) = ub(isam - 1) + getUnifRand(2, 7)")
            call disp%show("end do")
                            do isam = 2, nsam
                                lb(isam) = ub(isam - 1) + 1
                                ub(isam) = ub(isam - 1) + getUnifRand(2, 7)
                            end do
            call disp%show("lb")
            call disp%show( lb )
            call disp%show("ub")
            call disp%show( ub )
            call disp%show("ndim = getUnifRand(1, minval(ub - lb + 1, 1))")
                            ndim = getUnifRand(1, minval(ub - lb + 1, 1))
            call disp%show("call setRebound(cov, [1_IK, 1_IK, 0_IK], [ndim, ndim, nsam])")
                            call setRebound(cov, [1_IK, 1_IK, 0_IK], [ndim, ndim, nsam])
            call disp%show("call setRebound(mean, [1_IK, 0_IK], [ndim, nsam])")
                            call setRebound(mean, [1_IK, 0_IK], [ndim, nsam])
            call disp%show("call setResized(covMerged, [ndim, ndim])")
                            call setResized(covMerged, [ndim, ndim])
            call disp%show("call setResized(meanMerged, ndim)")
                            call setResized(meanMerged, ndim)
            call disp%show("sample = getUnifRand(-1., +1., ndim, ub(nsam))")
                            sample = getUnifRand(-1., +1., ndim, ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("iweight = getUnifRand(1, 10, size(sample, dim, IK))")
                            iweight = getUnifRand(1, 10, size(sample, dim, IK))
            call disp%show("iweight")
            call disp%show( iweight )
            call disp%show("cov(:,:,0) = getCov(sample, dim, iweight)")
                            cov(:,:,0) = getCov(sample, dim, iweight)
            call disp%show("cov(:,:,0) ! reference")
            call disp%show( cov(:,:,0) )
            call disp%show("mean(:,0) = getMean(sample, dim, iweight)")
                            mean(:,0) = getMean(sample, dim, iweight)
            call disp%show("mean(:,0) ! reference")
            call disp%show( mean(:,0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), dim, iweight(lb(isam):ub(isam)))")
            call disp%show("    mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim, iweight(lb(isam):ub(isam)))")
            call disp%show("end do")
                            do isam = 1, nsam
                                cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), dim, iweight(lb(isam):ub(isam)))
                                mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim, iweight(lb(isam):ub(isam)))
                            end do
            call disp%show("call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(iweight(:ub(1))), TKC) / real(sum(iweight), TKC), uppDia)")
                            call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(iweight(:ub(1))), TKC) / real(sum(iweight), TKC), uppDia)
            call disp%show("covMerged")
            call disp%show( covMerged )
            call disp%show("call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(iweight(:ub(1))), TKC) / real(sum(iweight), TKC), lowDia)")
                            call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(iweight(:ub(1))), TKC) / real(sum(iweight), TKC), lowDia)
            call disp%show("covMerged")
            call disp%show( covMerged )
            call disp%show("call setCovMeanMerged(cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(iweight(:ub(1))), TKC) / real(sum(iweight), TKC), uppDia)")
                            call setCovMeanMerged(cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(iweight(:ub(1))), TKC) / real(sum(iweight), TKC), uppDia)
            call disp%show("cov(:,:,2)")
            call disp%show( cov(:,:,2) )
            call disp%show("cov(:,:,0) ! reference")
            call disp%show( cov(:,:,0) )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased merged covariance of a reliability weighted multivariate sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKC), allocatable :: mean(:,:), cov(:,:,:), meanMerged(:), covMerged(:,:)
        real(TKC), allocatable :: sample(:,:)
        real(TKC), allocatable :: rweight(:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("dim = 2; lb(1) = 1; ub(1) = getUnifRand(2, 7)")
                            dim = 2; lb(1) = 1; ub(1) = getUnifRand(2, 7)
            call disp%show("do isam = 2, nsam")
            call disp%show("    lb(isam) = ub(isam - 1) + 1")
            call disp%show("    ub(isam) = ub(isam - 1) + getUnifRand(2, 7)")
            call disp%show("end do")
                            do isam = 2, nsam
                                lb(isam) = ub(isam - 1) + 1
                                ub(isam) = ub(isam - 1) + getUnifRand(2, 7)
                            end do
            call disp%show("lb")
            call disp%show( lb )
            call disp%show("ub")
            call disp%show( ub )
            call disp%show("ndim = getUnifRand(1, minval(ub - lb + 1, 1))")
                            ndim = getUnifRand(1, minval(ub - lb + 1, 1))
            call disp%show("call setRebound(cov, [1_IK, 1_IK, 0_IK], [ndim, ndim, nsam])")
                            call setRebound(cov, [1_IK, 1_IK, 0_IK], [ndim, ndim, nsam])
            call disp%show("call setRebound(mean, [1_IK, 0_IK], [ndim, nsam])")
                            call setRebound(mean, [1_IK, 0_IK], [ndim, nsam])
            call disp%show("call setResized(covMerged, [ndim, ndim])")
                            call setResized(covMerged, [ndim, ndim])
            call disp%show("call setResized(meanMerged, ndim)")
                            call setResized(meanMerged, ndim)
            call disp%show("sample = getUnifRand(-1., +1., ndim, ub(nsam))")
                            sample = getUnifRand(-1., +1., ndim, ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("rweight = getUnifRand(1., 2., size(sample, dim, IK))")
                            rweight = getUnifRand(1., 2., size(sample, dim, IK))
            call disp%show("rweight")
            call disp%show( rweight )
            call disp%show("cov(:,:,0) = getCov(sample, 2_IK, rweight)")
                            cov(:,:,0) = getCov(sample, 2_IK, rweight)
            call disp%show("cov(:,:,0) ! reference")
            call disp%show( cov(:,:,0) )
            call disp%show("mean(:,0) = getMean(sample, dim, rweight)")
                            mean(:,0) = getMean(sample, dim, rweight)
            call disp%show("mean(:,0) ! reference")
            call disp%show( mean(:,0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), 2_IK, rweight(lb(isam):ub(isam)))")
            call disp%show("    mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), 2_IK, rweight(lb(isam):ub(isam)))")
            call disp%show("end do")
                            do isam = 1, nsam
                                cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), 2_IK, rweight(lb(isam):ub(isam)))
                                mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), 2_IK, rweight(lb(isam):ub(isam)))
                            end do
            call disp%show("call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(rweight(:ub(1))), TKC) / real(sum(rweight), TKC), uppDia)")
                            call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(rweight(:ub(1))), TKC) / real(sum(rweight), TKC), uppDia)
            call disp%show("covMerged")
            call disp%show( covMerged )
            call disp%show("call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(rweight(:ub(1))), TKC) / real(sum(rweight), TKC), lowDia)")
                            call setCovMeanMerged(covMerged, meanMerged, cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(rweight(:ub(1))), TKC) / real(sum(rweight), TKC), lowDia)
            call disp%show("covMerged")
            call disp%show( covMerged )
            call disp%show("call setCovMeanMerged(cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(rweight(:ub(1))), TKC) / real(sum(rweight), TKC), uppDia)")
                            call setCovMeanMerged(cov(:,:,2), mean(:,2), cov(:,:,1), mean(:,1), real(sum(rweight(:ub(1))), TKC) / real(sum(rweight), TKC), uppDia)
            call disp%show("cov(:,:,2)")
            call disp%show( cov(:,:,2) )
            call disp%show("cov(:,:,0) ! reference")
            call disp%show( cov(:,:,0) )
            call disp%skip()
        end do
    end block

end program example