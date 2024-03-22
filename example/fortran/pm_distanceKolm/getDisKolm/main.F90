program example

    use pm_kind, only: SK, IK
    use pm_distUnif, only: getUnifRand
    use pm_arrayResize, only: setResized
    use pm_distanceKolm, only: ascending
    use pm_distanceKolm, only: getDisKolm
    use pm_arrayVerbose, only: getVerbose
    use pm_arraySort, only: getSorted
    use pm_arrayFill, only: getFilled
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    integer(IK) :: nsam1, nsam2
    integer(IK) :: itry, ntry = 10
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Kolmogorov distance of one sample w.r.t. Uniform or arbitrary distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! all other real kinds are also supported.
        integer(IK), allocatable :: iweight1(:)
        real(TKC), allocatable :: rweight1(:)
        real(TKC), allocatable :: sample1(:)
        real(TKC) :: disKolm
        do itry = 1, ntry
            call disp%show("nsam1 = getUnifRand(0, 10)")
                            nsam1 = getUnifRand(0, 10)
            call disp%show("sample1 = getUnifRand(0., 1., nsam1)")
                            sample1 = getUnifRand(0., 1., nsam1)
            call disp%show("sample1")
            call disp%show( sample1 )
            call disp%show("disKolm = getDisKolm(sample1) ! assuming unweighted samples.")
                            disKolm = getDisKolm(sample1) ! assuming unweighted samples.
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(sample1, getUnifCDF_RKS) ! assuming unweighted samples.")
                            disKolm = getDisKolm(sample1, getUnifCDF_RKS) ! assuming unweighted samples.
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(getSorted(sample1), ascending) ! assuming unweighted samples.")
                            disKolm = getDisKolm(getSorted(sample1), ascending) ! assuming unweighted samples.
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(getSorted(sample1), getUnifCDF_RKS, ascending) ! assuming unweighted samples.")
                            disKolm = getDisKolm(getSorted(sample1), getUnifCDF_RKS, ascending) ! assuming unweighted samples.
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%skip()
            call disp%show("iweight1 = getUnifRand(0, 9, nsam1)")
                            iweight1 = getUnifRand(0, 9, nsam1)
            call disp%show("iweight1")
            call disp%show( iweight1 )
            call disp%show("rweight1 = iweight1")
                            rweight1 = iweight1
            call disp%show("disKolm = getDisKolm(getVerbose(sample1, iweight1, sum(iweight1)))")
                            disKolm = getDisKolm(getVerbose(sample1, iweight1, sum(iweight1)))
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(getVerbose(sample1, iweight1, sum(iweight1)), getUnifCDF_RKS)")
                            disKolm = getDisKolm(getVerbose(sample1, iweight1, sum(iweight1)), getUnifCDF_RKS)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(sample1, iweight1, sum(iweight1))")
                            disKolm = getDisKolm(sample1, iweight1, sum(iweight1))
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(sample1, iweight1, sum(iweight1), getUnifCDF_RKS)")
                            disKolm = getDisKolm(sample1, iweight1, sum(iweight1), getUnifCDF_RKS)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(sample1, rweight1, sum(rweight1))")
                            disKolm = getDisKolm(sample1, rweight1, sum(rweight1))
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(sample1, rweight1, sum(rweight1), getUnifCDF_RKS)")
                            disKolm = getDisKolm(sample1, rweight1, sum(rweight1), getUnifCDF_RKS)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%skip()
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Kolmogorov distance of two samples.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! all other real kinds are also supported.
        integer(IK), allocatable :: iweight1(:), iweight2(:)
        real(TKC), allocatable :: rweight1(:), rweight2(:)
        real(TKC), allocatable :: sample1(:), sample2(:)
        real(TKC) :: disKolm
        do itry = 1, ntry
            call disp%show("nsam1 = getUnifRand(0, 10); nsam2 = getUnifRand(0, 10)")
                            nsam1 = getUnifRand(0, 10); nsam2 = getUnifRand(0, 10)
            call disp%show("sample1 = getUnifRand(0., 1., nsam1)")
                            sample1 = getUnifRand(0., 1., nsam1)
            call disp%show("sample1")
            call disp%show( sample1 )
            call disp%show("sample2 = getUnifRand(0., 1., nsam2)")
                            sample2 = getUnifRand(0., 1., nsam2)
            call disp%show("sample2")
            call disp%show( sample2 )
            call disp%show("disKolm = getDisKolm(sample1, sample2) ! assuming unweighted samples.")
                            disKolm = getDisKolm(sample1, sample2) ! assuming unweighted samples.
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(getSorted(sample1), getSorted(sample2), ascending) ! assuming unweighted samples.")
                            disKolm = getDisKolm(getSorted(sample1), getSorted(sample2), ascending) ! assuming unweighted samples.
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%skip()
            call disp%show("iweight1 = getUnifRand(0, 9, nsam1)")
                            iweight1 = getUnifRand(0, 9, nsam1)
            call disp%show("iweight1")
            call disp%show( iweight1 )
            call disp%show("rweight1 = iweight1")
                            rweight1 = iweight1
            call disp%show("disKolm = getDisKolm(getVerbose(sample1, iweight1, sum(iweight1)), sample2)")
                            disKolm = getDisKolm(getVerbose(sample1, iweight1, sum(iweight1)), sample2)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(sample1, iweight1, sum(iweight1), sample2)")
                            disKolm = getDisKolm(sample1, iweight1, sum(iweight1), sample2)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(sample1, rweight1, sum(rweight1), sample2)")
                            disKolm = getDisKolm(sample1, rweight1, sum(rweight1), sample2)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%skip()
            call disp%show("iweight2 = getUnifRand(0, 9, nsam2)")
                            iweight2 = getUnifRand(0, 9, nsam2)
            call disp%show("iweight2")
            call disp%show( iweight2 )
            call disp%show("rweight2 = iweight2")
                            rweight2 = iweight2
            call disp%show("disKolm = getDisKolm(getVerbose(sample1, iweight1, sum(iweight1)), getVerbose(sample2, iweight2, sum(iweight2)))")
                            disKolm = getDisKolm(getVerbose(sample1, iweight1, sum(iweight1)), getVerbose(sample2, iweight2, sum(iweight2)))
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(sample1, iweight1, sum(iweight1), sample2, iweight2, sum(iweight2))")
                            disKolm = getDisKolm(sample1, iweight1, sum(iweight1), sample2, iweight2, sum(iweight2))
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("disKolm = getDisKolm(sample1, rweight1, sum(rweight1), sample2, rweight2, sum(rweight2))")
                            disKolm = getDisKolm(sample1, rweight1, sum(rweight1), sample2, rweight2, sum(rweight2))
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%skip()
        end do
    end block

contains

    function getUnifCDF_RKS(x) result(cdf)
        use pm_distUnif, only: getUnifCDF
        use pm_kind, only: RKC => RKS
        real(RKC), intent(in) :: x
        real(RKC) :: cdf
        cdf = getUnifCDF(x)
    end function

end program example