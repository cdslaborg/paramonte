program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: IKG => IKS ! all processor kinds are supported.
    use pm_sampleECDF, only: setECDF
    use pm_arrayRange, only: getRange
    use pm_arrayResize, only: setResized
    use pm_arrayShuffle, only: getShuffled
    use pm_distUnif, only: getUnifRand
    use pm_arraySort, only: setSorted
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    integer(IK) :: itry, ntry = 10
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: RKG => RKS ! all processor kinds are supported.
        real(RKG), allocatable :: ecdf(:), lcdf(:), ucdf(:)
        integer(IK), allocatable :: iweight(:)
        real(RKG), allocatable :: rweight(:)
        do itry = 1, ntry
            call disp%show("call setResized(ecdf, getUnifRand(1_IK, 10_IK))")
                            call setResized(ecdf, getUnifRand(1_IK, 10_IK))
            call disp%show("call setResized(lcdf, size(ecdf, 1, IK))")
                            call setResized(lcdf, size(ecdf, 1, IK))
            call disp%show("call setResized(ucdf, size(ecdf, 1, IK))")
                            call setResized(ucdf, size(ecdf, 1, IK))
            call disp%show("call setECDF(ecdf, lcdf, ucdf)")
                            call setECDF(ecdf, lcdf, ucdf)
            call disp%show("lcdf")
            call disp%show( lcdf )
            call disp%show("ecdf")
            call disp%show( ecdf )
            call disp%show("ucdf")
            call disp%show( ucdf )
            call disp%show("call setECDF(ecdf)")
                            call setECDF(ecdf)
            call disp%show("ecdf")
            call disp%show( ecdf )
            call disp%skip()
            call disp%show("iweight = getUnifRand(1, 10, size(ecdf, 1, IK))")
                            iweight = getUnifRand(1, 10, size(ecdf, 1, IK))
            call disp%show("iweight")
            call disp%show( iweight )
            call disp%show("call setECDF(ecdf, iweight, sum(iweight), lcdf, ucdf)")
                            call setECDF(ecdf, iweight, sum(iweight), lcdf, ucdf)
            call disp%show("lcdf")
            call disp%show( lcdf )
            call disp%show("ecdf")
            call disp%show( ecdf )
            call disp%show("ucdf")
            call disp%show( ucdf )
            call disp%skip()
            call disp%show("rweight = iweight")
                            rweight = iweight
            call disp%show("rweight")
            call disp%show( rweight )
            call disp%show("call setECDF(ecdf, rweight, sum(rweight))")
                            call setECDF(ecdf, rweight, sum(rweight))
            call disp%show("ecdf")
            call disp%show( ecdf )
            call disp%skip()
        end do
    end block

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Compute the ecdf for varying sizes of a Normal sample.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_val2str, only: getStr
        use pm_kind, only: RKG => RKS
        use pm_distNorm, only: setNormRand
        character(:, SK), allocatable :: fileName
        real(RKG), allocatable :: sample(:), ecdf(:), lcdf(:), ucdf(:)
        integer(IK) :: i, isam, nsam
        integer(IK) :: fileUnit
        do i = 1, 4
            nsam = 10 ** i
            call setResized(lcdf, nsam)
            call setResized(ecdf, nsam)
            call setResized(ucdf, nsam)
            call setResized(sample, nsam)
            call setNormRand(sample)
            call setSorted(sample)
            call setECDF(ecdf, lcdf, ucdf, alpha = 0.98_RKG)
            fileName = "main.norm." // getStr(nsam) // ".out"
            open(newunit = fileUnit, file = fileName, status = "replace")
                write(fileUnit, "(*(g0,:,','))") "sample,ecdf,lcdf,ucdf"
                do isam = 1, nsam
                    write(fileUnit,"(*(g0,:,','))") sample(isam), ecdf(isam), lcdf(isam), ucdf(isam)
                end do
            close(fileUnit)
        end do
    end block

end program example