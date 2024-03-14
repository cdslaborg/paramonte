program example

    use pm_kind, only: SK, IK
    use pm_distNorm, only: getNormRand
    use pm_distUnif, only: getUnifRand
    use pm_arrayResize, only: setResized
    use pm_distanceKolm, only: ascending
    use pm_distanceKolm, only: getDisKolm
    use pm_arrayVerbose, only: getVerbose
    use pm_arraySort, only: getSorted
    use pm_arrayFill, only: getFilled
    use pm_statest, only: setProbKS
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    integer(IK) :: nsam1, nsam2
    integer(IK) :: itry, ntry = 10
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the KS probability of sample originating from a Uniform distribution in range `[0, 1)`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! all other real kinds are also acceptable: RK, RK32, RK64, RK128
        integer(IK), allocatable :: iweight1(:)
        real(TKC), allocatable :: rweight1(:)
        real(TKC), allocatable :: sample1(:)
        real(TKC) :: probKS, quanKS, statKS
        do itry = 1, ntry
            call disp%show("nsam1 = getUnifRand(1, 10)")
                            nsam1 = getUnifRand(1, 10)
            call disp%show("sample1 = getUnifRand(0., 1., nsam1)")
                            sample1 = getUnifRand(0., 1., nsam1)
            call disp%show("sample1")
            call disp%show( sample1 )
            call disp%show("statKS = getDisKolm(sample1) ! assuming unweighted samples.")
                            statKS = getDisKolm(sample1) ! assuming unweighted samples.
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, nsam1) ! assuming unweighted samples.")
                            call setProbKS(probKS, quanKS, statKS, nsam1) ! assuming unweighted samples.
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%skip()
            call disp%show("iweight1 = getUnifRand(1, 9, nsam1)")
                            iweight1 = getUnifRand(1, 9, nsam1)
            call disp%show("iweight1")
            call disp%show( iweight1 )
            call disp%show("rweight1 = iweight1")
                            rweight1 = iweight1
            call disp%show("statKS = getDisKolm(sample1, iweight1, sum(iweight1))")
                            statKS = getDisKolm(sample1, iweight1, sum(iweight1))
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, sum(iweight1))")
                            call setProbKS(probKS, quanKS, statKS, sum(iweight1))
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%show("statKS = getDisKolm(sample1, rweight1, sum(rweight1))")
                            statKS = getDisKolm(sample1, rweight1, sum(rweight1))
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, sum(rweight1), sum(rweight1**2))")
                            call setProbKS(probKS, quanKS, statKS, sum(rweight1), sum(rweight1**2))
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the KS probability of a sample against a Normal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! all other real kinds are also acceptable: RK, RK32, RK64, RK128
        integer(IK), allocatable :: iweight1(:)
        real(TKC), allocatable :: rweight1(:)
        real(TKC), allocatable :: sample1(:)
        real(TKC) :: probKS, quanKS, statKS
        do itry = 1, ntry
            call disp%show("nsam1 = getUnifRand(5, 10)")
                            nsam1 = getUnifRand(5, 10)
            call disp%show("sample1 = getNormRand(mean = getFilled(0., nsam1))")
                            sample1 = getNormRand(mean = getFilled(0., nsam1))
            call disp%show("sample1")
            call disp%show( sample1 )
            call disp%show("statKS = getDisKolm(sample1, getNormCDF_RKS) ! assuming unweighted samples.")
                            statKS = getDisKolm(sample1, getNormCDF_RKS) ! assuming unweighted samples.
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, nsam1) ! assuming unweighted samples.")
                            call setProbKS(probKS, quanKS, statKS, nsam1) ! assuming unweighted samples.
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%skip()
            call disp%show("iweight1 = getUnifRand(1, 9, nsam1)")
                            iweight1 = getUnifRand(1, 9, nsam1)
            call disp%show("iweight1")
            call disp%show( iweight1 )
            call disp%show("rweight1 = iweight1")
                            rweight1 = iweight1
            call disp%show("statKS = getDisKolm(sample1, iweight1, sum(iweight1), getNormCDF_RKS)")
                            statKS = getDisKolm(sample1, iweight1, sum(iweight1), getNormCDF_RKS)
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, sum(iweight1))")
                            call setProbKS(probKS, quanKS, statKS, sum(iweight1))
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%show("statKS = getDisKolm(sample1, rweight1, sum(rweight1), getNormCDF_RKS)")
                            statKS = getDisKolm(sample1, rweight1, sum(rweight1), getNormCDF_RKS)
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, sum(rweight1), sum(rweight1**2))")
                            call setProbKS(probKS, quanKS, statKS, sum(rweight1), sum(rweight1**2))
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the two sample KS probability.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! all other real kinds are also acceptable: RK, RK32, RK64, RK128
        integer(IK), allocatable :: iweight1(:), iweight2(:)
        real(TKC), allocatable :: rweight1(:), rweight2(:)
        real(TKC), allocatable :: sample1(:), sample2(:)
        real(TKC) :: probKS, quanKS, statKS
        do itry = 1, ntry
            call disp%show("nsam1 = getUnifRand(1, 10); nsam2 = getUnifRand(1, 10)")
                            nsam1 = getUnifRand(1, 10); nsam2 = getUnifRand(1, 10)
            call disp%show("sample1 = getUnifRand(0., 1., nsam1)")
                            sample1 = getUnifRand(0., 1., nsam1)
            call disp%show("sample1")
            call disp%show( sample1 )
            call disp%show("sample2 = getUnifRand(0., 1., nsam2)")
                            sample2 = getUnifRand(0., 1., nsam2)
            call disp%show("sample2")
            call disp%show( sample2 )
            call disp%show("statKS = getDisKolm(sample1, sample2) ! assuming unweighted samples.")
                            statKS = getDisKolm(sample1, sample2) ! assuming unweighted samples.
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, nsam1, nsam2) ! assuming unweighted samples.")
                            call setProbKS(probKS, quanKS, statKS, nsam1, nsam2) ! assuming unweighted samples.
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%skip()
            call disp%show("iweight1 = getUnifRand(1, 9, nsam1)")
                            iweight1 = getUnifRand(1, 9, nsam1)
            call disp%show("iweight1")
            call disp%show( iweight1 )
            call disp%show("rweight1 = iweight1")
                            rweight1 = iweight1
            call disp%show("statKS = getDisKolm(sample1, iweight1, sum(iweight1), sample2)")
                            statKS = getDisKolm(sample1, iweight1, sum(iweight1), sample2)
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, sum(iweight1), nsam2)")
                            call setProbKS(probKS, quanKS, statKS, sum(iweight1), nsam2)
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%show("statKS = getDisKolm(sample1, rweight1, sum(rweight1), sample2)")
                            statKS = getDisKolm(sample1, rweight1, sum(rweight1), sample2)
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, sum(rweight1), nsam2, sum(rweight1**2))")
                            call setProbKS(probKS, quanKS, statKS, sum(rweight1), nsam2, sum(rweight1**2))
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%skip()
            call disp%show("iweight2 = getUnifRand(1, 9, nsam2)")
                            iweight2 = getUnifRand(1, 9, nsam2)
            call disp%show("iweight2")
            call disp%show( iweight2 )
            call disp%show("rweight2 = iweight2")
                            rweight2 = iweight2
            call disp%show("statKS = getDisKolm(sample1, iweight1, sum(iweight1), sample2, iweight2, sum(iweight2))")
                            statKS = getDisKolm(sample1, iweight1, sum(iweight1), sample2, iweight2, sum(iweight2))
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, sum(iweight1), sum(iweight2))")
                            call setProbKS(probKS, quanKS, statKS, sum(iweight1), sum(iweight2))
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%show("call setProbKS(probKS, quanKS, statKS, sum(rweight1), sum(iweight2), sum(rweight1**2))")
                            call setProbKS(probKS, quanKS, statKS, sum(rweight1), sum(iweight2), sum(rweight1**2))
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%show("statKS = getDisKolm(sample1, rweight1, sum(rweight1), sample2, rweight2, sum(rweight2))")
                            statKS = getDisKolm(sample1, rweight1, sum(rweight1), sample2, rweight2, sum(rweight2))
            call disp%show("statKS")
            call disp%show( statKS )
            call disp%show("call setProbKS(probKS, quanKS, statKS, sum(rweight1), sum(rweight2), sum(rweight1**2), sum(rweight2**2))")
                            call setProbKS(probKS, quanKS, statKS, sum(rweight1), sum(rweight2), sum(rweight1**2), sum(rweight2**2))
            call disp%show("[probKS, quanKS]")
            call disp%show( [probKS, quanKS] )
            call disp%skip()
        end do
    end block

contains

    function getNormCDF_RKS(x) result(cdf)
        use pm_distNorm, only: getNormCDF
        use pm_kind, only: RKC => RKS
        real(RKC), intent(in) :: x
        real(RKC) :: cdf
        cdf = getNormCDF(x)
    end function

    function getUnifCDF_RKS(x) result(cdf)
        use pm_distUnif, only: getUnifCDF
        use pm_kind, only: RKC => RKS
        real(RKC), intent(in) :: x
        real(RKC) :: cdf
        cdf = getUnifCDF(x)
    end function

end program example