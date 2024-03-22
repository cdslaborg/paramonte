program example

    use pm_kind, only: SK, IK, RK
    use pm_kind, only: TKC => RKS ! All other real types are also supported.
    use pm_sampleVar, only: getVar
    use pm_sampleMean, only: getMean
    use pm_sampleVar, only: getVarMerged
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
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased merged variance of a univariate sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKC) :: mean(1:nsam), var(0:nsam), varMerged
        real(TKC), allocatable :: sample(:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("lb(1) = 1; ub(1) = getUnifRand(2, 7)")
                            lb(1) = 1; ub(1) = getUnifRand(2, 7)
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
            call disp%show("sample = getUnifRand(0., 1., ub(nsam))")
                            sample = getUnifRand(0., 1., ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("var(0) = getVar(sample)")
                            var(0) = getVar(sample)
            call disp%show("var(0) ! reference")
            call disp%show( var(0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    var(isam) = getVar(sample(lb(isam):ub(isam)))")
            call disp%show("    mean(isam) = getMean(sample(lb(isam):ub(isam)))")
            call disp%show("end do")
                            do isam = 1, nsam
                                var(isam) = getVar(sample(lb(isam):ub(isam)))
                                mean(isam) = getMean(sample(lb(isam):ub(isam)))
                            end do
            call disp%show("varMerged = getVarMerged(var(2), var(1), mean(1) - mean(2), ub(1) / real(ub(2), TKC))")
                            varMerged = getVarMerged(var(2), var(1), mean(1) - mean(2), ub(1) / real(ub(2), TKC))
            call disp%show("varMerged")
            call disp%show( varMerged )
            call disp%show("var(0) ! reference")
            call disp%show( var(0) )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased merged variance of a frequency weighted univariate sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKC) :: mean(1:nsam), var(0:nsam), varMerged
        real(TKC), allocatable :: sample(:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("lb(1) = 1; ub(1) = getUnifRand(2, 7)")
                            lb(1) = 1; ub(1) = getUnifRand(2, 7)
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
            call disp%show("sample = getUnifRand(0., 1., ub(nsam))")
                            sample = getUnifRand(0., 1., ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("iweight = getUnifRand(1, 10, size(sample, 1, IK))")
                            iweight = getUnifRand(1, 10, size(sample, 1, IK))
            call disp%show("iweight")
            call disp%show( iweight )
            call disp%show("var(0) = getVar(sample, iweight)")
                            var(0) = getVar(sample, iweight)
            call disp%show("var(0) ! reference")
            call disp%show( var(0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    var(isam) = getVar(sample(lb(isam):ub(isam)), iweight(lb(isam):ub(isam)))")
            call disp%show("    mean(isam) = getMean(sample(lb(isam):ub(isam)), iweight(lb(isam):ub(isam)))")
            call disp%show("end do")
                            do isam = 1, nsam
                                var(isam) = getVar(sample(lb(isam):ub(isam)), iweight(lb(isam):ub(isam)))
                                mean(isam) = getMean(sample(lb(isam):ub(isam)), iweight(lb(isam):ub(isam)))
                            end do
            call disp%show("varMerged = getVarMerged(var(2), var(1), mean(1) - mean(2), real(sum(iweight(:ub(1))), TKC) / real(sum(iweight), TKC))")
                            varMerged = getVarMerged(var(2), var(1), mean(1) - mean(2), real(sum(iweight(:ub(1))), TKC) / real(sum(iweight), TKC))
            call disp%show("varMerged")
            call disp%show( varMerged )
            call disp%show("var(0) ! reference")
            call disp%show( var(0) )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased merged variance of a reliability weighted univariate sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKC) :: mean(1:nsam), var(0:nsam), varMerged
        real(TKC), allocatable :: sample(:), rweight(:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("lb(1) = 1; ub(1) = getUnifRand(2, 7)")
                            lb(1) = 1; ub(1) = getUnifRand(2, 7)
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
            call disp%show("sample = getUnifRand(0., 1., ub(nsam))")
                            sample = getUnifRand(0., 1., ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("rweight = getUnifRand(1., 2., size(sample, 1, IK))")
                            rweight = getUnifRand(1., 2., size(sample, 1, IK))
            call disp%show("rweight")
            call disp%show( rweight )
            call disp%show("var(0) = getVar(sample, rweight)")
                            var(0) = getVar(sample, rweight)
            call disp%show("var(0) ! reference")
            call disp%show( var(0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    var(isam) = getVar(sample(lb(isam):ub(isam)), rweight(lb(isam):ub(isam)))")
            call disp%show("    mean(isam) = getMean(sample(lb(isam):ub(isam)), rweight(lb(isam):ub(isam)))")
            call disp%show("end do")
                            do isam = 1, nsam
                                var(isam) = getVar(sample(lb(isam):ub(isam)), rweight(lb(isam):ub(isam)))
                                mean(isam) = getMean(sample(lb(isam):ub(isam)), rweight(lb(isam):ub(isam)))
                            end do
            call disp%show("varMerged = getVarMerged(var(2), var(1), mean(1) - mean(2), real(sum(rweight(:ub(1))), TKC) / real(sum(rweight), TKC))")
                            varMerged = getVarMerged(var(2), var(1), mean(1) - mean(2), real(sum(rweight(:ub(1))), TKC) / real(sum(rweight), TKC))
            call disp%show("varMerged")
            call disp%show( varMerged )
            call disp%show("var(0) ! reference")
            call disp%show( var(0) )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased merged variance of a multivariate sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKC), allocatable :: mean(:,:), var(:,:), varMerged(:)
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
            call disp%show("call setRebound(var, [1_IK, 0_IK], [ndim, nsam])")
                            call setRebound(var, [1_IK, 0_IK], [ndim, nsam])
            call disp%show("call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])")
                            call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])
            call disp%show("call setResized(varMerged, ndim)")
                            call setResized(varMerged, ndim)
            call disp%show("sample = getUnifRand(-1., +1., ndim, ub(nsam))")
                            sample = getUnifRand(-1., +1., ndim, ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("var(:,0) = getVar(sample, dim)")
                            var(:,0) = getVar(sample, dim)
            call disp%show("var(:,0) ! reference")
            call disp%show( var(:,0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    var(:,isam) = getVar(sample(:,lb(isam):ub(isam)), dim)")
            call disp%show("    mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim)")
            call disp%show("end do")
                            do isam = 1, nsam
                                var(:,isam) = getVar(sample(:,lb(isam):ub(isam)), dim)
                                mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), dim)
                            end do
            call disp%show("varMerged = getVarMerged(var(:,2), var(:,1), mean(:,1) - mean(:,2), real(ub(1), TKC) / real(ub(2), TKC))")
                            varMerged = getVarMerged(var(:,2), var(:,1), mean(:,1) - mean(:,2), real(ub(1), TKC) / real(ub(2), TKC))
            call disp%show("varMerged")
            call disp%show( varMerged )
            call disp%show("var(:,0) ! reference")
            call disp%show( var(:,0) )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased merged variance of a frequency weighted multivariate sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKC), allocatable :: mean(:,:), var(:,:), varMerged(:)
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
            call disp%show("call setRebound(var, [1_IK, 0_IK], [ndim, nsam])")
                            call setRebound(var, [1_IK, 0_IK], [ndim, nsam])
            call disp%show("call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])")
                            call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])
            call disp%show("call setResized(varMerged, ndim)")
                            call setResized(varMerged, ndim)
            call disp%show("sample = getUnifRand(-1., +1., ndim, ub(nsam))")
                            sample = getUnifRand(-1., +1., ndim, ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("iweight = getUnifRand(1, 10, size(sample, dim, IK))")
                            iweight = getUnifRand(1, 10, size(sample, dim, IK))
            call disp%show("iweight")
            call disp%show( iweight )
            call disp%show("var(:,0) = getVar(sample, 2_IK, iweight)")
                            var(:,0) = getVar(sample, 2_IK, iweight)
            call disp%show("var(:,0) ! reference")
            call disp%show( var(:,0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    var(:,isam) = getVar(sample(:,lb(isam):ub(isam)), 2_IK, iweight(lb(isam):ub(isam)))")
            call disp%show("    mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), 2_IK, iweight(lb(isam):ub(isam)))")
            call disp%show("end do")
                            do isam = 1, nsam
                                var(:,isam) = getVar(sample(:,lb(isam):ub(isam)), 2_IK, iweight(lb(isam):ub(isam)))
                                mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), 2_IK, iweight(lb(isam):ub(isam)))
                            end do
            call disp%show("varMerged = getVarMerged(var(:,2), var(:,1), mean(:,1) - mean(:,2), real(sum(iweight(:ub(1))), TKC) / real(sum(iweight), TKC))")
                            varMerged = getVarMerged(var(:,2), var(:,1), mean(:,1) - mean(:,2), real(sum(iweight(:ub(1))), TKC) / real(sum(iweight), TKC))
            call disp%show("varMerged")
            call disp%show( varMerged )
            call disp%show("var(:,0) ! reference")
            call disp%show( var(:,0) )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the biased merged variance of a reliability weighted multivariate sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKC), allocatable :: sample(:,:), mean(:,:)
        real(TKC), allocatable :: rweight(:), var(:,:), varMerged(:)
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
            call disp%show("call setRebound(var, [1_IK, 0_IK], [ndim, nsam])")
                            call setRebound(var, [1_IK, 0_IK], [ndim, nsam])
            call disp%show("call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])")
                            call setRebound(mean, [1_IK, 1_IK], [ndim, nsam])
            call disp%show("call setResized(varMerged, ndim)")
                            call setResized(varMerged, ndim)
            call disp%show("sample = getUnifRand(-1., +1., ndim, ub(nsam))")
                            sample = getUnifRand(-1., +1., ndim, ub(nsam))
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("rweight = getUnifRand(1., 2., size(sample, dim, IK))")
                            rweight = getUnifRand(1., 2., size(sample, dim, IK))
            call disp%show("rweight")
            call disp%show( rweight )
            call disp%show("var(:,0) = getVar(sample, 2_IK, rweight)")
                            var(:,0) = getVar(sample, 2_IK, rweight)
            call disp%show("var(:,0) ! reference")
            call disp%show( var(:,0) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    var(:,isam) = getVar(sample(:,lb(isam):ub(isam)), 2_IK, rweight(lb(isam):ub(isam)))")
            call disp%show("    mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), 2_IK, rweight(lb(isam):ub(isam)))")
            call disp%show("end do")
                            do isam = 1, nsam
                                var(:,isam) = getVar(sample(:,lb(isam):ub(isam)), 2_IK, rweight(lb(isam):ub(isam)))
                                mean(:,isam) = getMean(sample(:,lb(isam):ub(isam)), 2_IK, rweight(lb(isam):ub(isam)))
                            end do
            call disp%show("varMerged = getVarMerged(var(:,2), var(:,1), mean(:,1) - mean(:,2), real(sum(rweight(:ub(1))), TKC) / real(sum(rweight), TKC))")
                            varMerged = getVarMerged(var(:,2), var(:,1), mean(:,1) - mean(:,2), real(sum(rweight(:ub(1))), TKC) / real(sum(rweight), TKC))
            call disp%show("varMerged")
            call disp%show( varMerged )
            call disp%show("var(:,0) ! reference")
            call disp%show( var(:,0) )
            call disp%skip()
        end do
    end block

end program example