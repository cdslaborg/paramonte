program example

    use pm_kind, only: SK, IK
    use pm_kind, only: TKG => RKS ! All other real types are also supported.
    use pm_sampleCov, only: getCov
    use pm_sampleMean, only: getMean
    use pm_sampleCov, only: uppDia, lowDia
    use pm_sampleCov, only: setCovMerged
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
        real(TKG), allocatable :: mean(:,:), cov(:,:,:), covMerged(:,:)
        real(TKG), allocatable :: sample(:,:)
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
            call disp%show("call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])")
                            call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])
            call disp%show("call setResized(covMerged, [ndim, ndim])")
                            call setResized(covMerged, [ndim, ndim])
            call disp%show("sample = getUnifRand(-1., +1., ndim, ub(nsam))")
                            sample = getUnifRand(-1., +1., ndim, ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("cov(:,:,0) = getCov(sample, dim)")
                            cov(:,:,0) = getCov(sample, dim)
            call disp%show("cov(:,:,0) ! reference")
            call disp%show( cov(:,:,0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), dim)")
            call disp%show("    mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim)")
            call disp%show("end do")
                            do isam = 1, nsam
                                cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), dim)
                                mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim)
                            end do
            call disp%show("call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(ub(1), TKG) / real(ub(2), TKG), uppDia)")
                            call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(ub(1), TKG) / real(ub(2), TKG), uppDia)
            call disp%show("covMerged")
            call disp%show( covMerged )
            call disp%show("call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(ub(1), TKG) / real(ub(2), TKG), lowDia)")
                            call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(ub(1), TKG) / real(ub(2), TKG), lowDia)
            call disp%show("covMerged")
            call disp%show( covMerged )
            call disp%show("call setCovMerged(cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(ub(1), TKG) / real(ub(2), TKG), uppDia)")
                            call setCovMerged(cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(ub(1), TKG) / real(ub(2), TKG), uppDia)
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
        real(TKG), allocatable :: mean(:,:), cov(:,:,:), covMerged(:,:)
        real(TKG), allocatable :: sample(:,:)
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
            call disp%show("call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])")
                            call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])
            call disp%show("call setResized(covMerged, [ndim, ndim])")
                            call setResized(covMerged, [ndim, ndim])
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
            call disp%show("do isam = 1, nsam")
            call disp%show("    cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), dim, iweight(lb(isam):ub(isam)))")
            call disp%show("    mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim, iweight(lb(isam):ub(isam)))")
            call disp%show("end do")
                            do isam = 1, nsam
                                cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), dim, iweight(lb(isam):ub(isam)))
                                mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim, iweight(lb(isam):ub(isam)))
                            end do
            call disp%show("call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(iweight(:ub(1))), TKG) / real(sum(iweight), TKG), uppDia)")
                            call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(iweight(:ub(1))), TKG) / real(sum(iweight), TKG), uppDia)
            call disp%show("covMerged")                                                                                                                          
            call disp%show( covMerged )                                                                                                                          
            call disp%show("call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(iweight(:ub(1))), TKG) / real(sum(iweight), TKG), lowDia)")
                            call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(iweight(:ub(1))), TKG) / real(sum(iweight), TKG), lowDia)
            call disp%show("covMerged")
            call disp%show( covMerged )
            call disp%show("call setCovMerged(cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(iweight(:ub(1))), TKG) / real(sum(iweight), TKG), uppDia)")
                            call setCovMerged(cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(iweight(:ub(1))), TKG) / real(sum(iweight), TKG), uppDia)
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
        real(TKG), allocatable :: mean(:,:), cov(:,:,:), covMerged(:,:)
        real(TKG), allocatable :: sample(:,:)
        real(TKG), allocatable :: rweight(:)
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
            call disp%show("call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])")
                            call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])
            call disp%show("call setResized(covMerged, [ndim, ndim])")
                            call setResized(covMerged, [ndim, ndim])
            call disp%show("sample = getUnifRand(-1., +1., ndim, ub(nsam))")
                            sample = getUnifRand(-1., +1., ndim, ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("rweight = getUnifRand(1., 2., size(sample, dim, IK))")
                            rweight = getUnifRand(1., 2., size(sample, dim, IK))
            call disp%show("rweight")
            call disp%show( rweight )
            call disp%show("cov(:,:,0) = getCov(sample, dim, rweight)")
                            cov(:,:,0) = getCov(sample, dim, rweight)
            call disp%show("cov(:,:,0) ! reference")
            call disp%show( cov(:,:,0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), dim, rweight(lb(isam):ub(isam)))")
            call disp%show("    mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim, rweight(lb(isam):ub(isam)))")
            call disp%show("end do")
                            do isam = 1, nsam
                                cov(:,:,isam) = getCov(sample(:,lb(isam):ub(isam)), dim, rweight(lb(isam):ub(isam)))
                                mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim, rweight(lb(isam):ub(isam)))
                            end do
            call disp%show("call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(rweight(:ub(1))), TKG) / real(sum(rweight), TKG), uppDia)")
                            call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(rweight(:ub(1))), TKG) / real(sum(rweight), TKG), uppDia)
            call disp%show("covMerged")                                                                                                                          
            call disp%show( covMerged )                                                                                                                          
            call disp%show("call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(rweight(:ub(1))), TKG) / real(sum(rweight), TKG), lowDia)")
                            call setCovMerged(covMerged, cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(rweight(:ub(1))), TKG) / real(sum(rweight), TKG), lowDia)
            call disp%show("covMerged")
            call disp%show( covMerged )
            call disp%show("call setCovMerged(cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(rweight(:ub(1))), TKG) / real(sum(rweight), TKG), lowDia)")
                            call setCovMerged(cov(:,:,2), cov(:,:,1), mean(:,1) - mean(:,2), real(sum(rweight(:ub(1))), TKG) / real(sum(rweight), TKG), lowDia)
            call disp%show("cov(:,:,2)")
            call disp%show( cov(:,:,2) )
            call disp%show("cov(:,:,0) ! reference")
            call disp%show( cov(:,:,0) )
            call disp%skip()
        end do
    end block

end program example