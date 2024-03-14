program example

    use pm_kind, only: SK, IK
    use pm_distUnif, only: getUnifRand
    use pm_arrayResize, only: setResized
    use pm_distanceKolm, only: ascending
    use pm_distanceKolm, only: setDisKolm
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
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Kolmogorov distance of two unweighted samples.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! all other real kinds are also acceptable: RK, RK32, RK64, RK128
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
            call disp%show("call setDisKolm(disKolm, sample1, sample2)")
                            call setDisKolm(disKolm, sample1, sample2)
            call disp%show("sample1")
            call disp%show( sample1 )
            call disp%show("sample2")
            call disp%show( sample2 )
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("call setDisKolm(disKolm, getSorted(sample1), getSorted(sample2), ascending)")
                            call setDisKolm(disKolm, getSorted(sample1), getSorted(sample2), ascending)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%skip()
            call disp%show("iweight1 = getUnifRand(0, 9, nsam1)")
                            iweight1 = getUnifRand(0, 9, nsam1)
            call disp%show("iweight1")
            call disp%show( iweight1 )
            call disp%show("rweight1 = iweight1")
                            rweight1 = iweight1
            call disp%show("call setDisKolm(disKolm, getSorted(getVerbose(sample1, iweight1, sum(iweight1))), getSorted(sample2), ascending)")
                            call setDisKolm(disKolm, getSorted(getVerbose(sample1, iweight1, sum(iweight1))), getSorted(sample2), ascending)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("call setDisKolm(disKolm, sample1, iweight1, sum(iweight1), sample2)")
                            call setDisKolm(disKolm, sample1, iweight1, sum(iweight1), sample2)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("call setDisKolm(disKolm, sample1, rweight1, sum(rweight1), sample2)")
                            call setDisKolm(disKolm, sample1, rweight1, sum(rweight1), sample2)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%skip()
            call disp%show("iweight1 = getUnifRand(0, 9, nsam1)")
                            iweight1 = getUnifRand(0, 9, nsam1)
            call disp%show("iweight1")
            call disp%show( iweight1 )
            call disp%show("iweight2 = getUnifRand(0, 9, nsam2)")
                            iweight2 = getUnifRand(0, 9, nsam2)
            call disp%show("iweight2")
            call disp%show( iweight2 )
            call disp%show("rweight1 = iweight1; rweight2 = iweight2")
                            rweight1 = iweight1; rweight2 = iweight2
            call disp%show("call setDisKolm(disKolm, getSorted(getVerbose(sample1, iweight1, sum(iweight1))), getSorted(getVerbose(sample2, iweight2, sum(iweight2))), ascending)")
                            call setDisKolm(disKolm, getSorted(getVerbose(sample1, iweight1, sum(iweight1))), getSorted(getVerbose(sample2, iweight2, sum(iweight2))), ascending)
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("call setDisKolm(disKolm, sample1, iweight1, sum(iweight1), sample2, iweight2, sum(iweight2))")
                            call setDisKolm(disKolm, sample1, iweight1, sum(iweight1), sample2, iweight2, sum(iweight2))
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%show("call setDisKolm(disKolm, sample1, rweight1, sum(rweight1), sample2, rweight2, sum(rweight2))")
                            call setDisKolm(disKolm, sample1, rweight1, sum(rweight1), sample2, rweight2, sum(rweight2))
            call disp%show("disKolm")
            call disp%show( disKolm )
            call disp%skip()
        end do
    end block

end program example