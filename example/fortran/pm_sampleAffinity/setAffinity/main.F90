program example

    use pm_kind, only: SK, IK, LK
    use pm_matrixChol, only: getMatChol, uppDia, lowDia
    use pm_sampleAffinity, only: upperDiag, upperUnit
    use pm_sampleAffinity, only: lowerDiag, lowerUnit
    use pm_matrixInit, only: setMatInit, lowDia
    use pm_sampleAffinity, only: setAffinity
    use pm_sampleAffinity, only: genrecmat
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_arrayFill, only: getFilled
    use pm_matrixInv, only: getMatInv
    use pm_distCov, only: getCovRand
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    integer(IK) :: itry, ntry = 5
    integer(IK) :: dim, isam, ndim, nsam
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Transform 2D sample along the second dimension.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKG => RKS
        real(RKG), allocatable :: tlate(:), tform(:,:), sample(:,:), affinInv(:,:), affinity(:,:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("dim = 2; ndim = getUnifRand(1, 3); nsam = getUnifRand(4, 5)")
                            dim = 2; ndim = getUnifRand(1, 3); nsam = getUnifRand(4, 5)
            call disp%show("[dim, ndim, nsam]")
            call disp%show( [dim, ndim, nsam] )
            call disp%show("tlate = getUnifRand(-5, +5, ndim)")
                            tlate = getUnifRand(-5, +5, ndim)
            call disp%show("tlate")
            call disp%show( tlate )
            call disp%show("sample = getUnifRand(-9, +9, ndim, nsam)")
                            sample = getUnifRand(-9, +9, ndim, nsam)
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("call setResized(affinity, shape(sample, IK))")
                            call setResized(affinity, shape(sample, IK))
            call disp%show("call setResized(affinInv, shape(sample, IK))")
                            call setResized(affinInv, shape(sample, IK))
            call disp%show("tform = getCovRand(1., ndim)")
                            tform = getCovRand(1., ndim)
            call disp%show("tform")
            call disp%show( tform )
            call disp%show("call setAffinity(affinity, sample, dim, tform, genrecmat) ! whole sample transformation.")
                            call setAffinity(affinity, sample, dim, tform, genrecmat) ! whole sample transformation.
            call disp%show("transpose(affinity)")
            call disp%show( transpose(affinity) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    call setAffinity(affinity(:,isam), sample(:,isam), tform, genrecmat) ! point-wise transformation")
            call disp%show("end do")
                            do isam = 1, nsam
                                call setAffinity(affinity(:,isam), sample(:,isam), tform, genrecmat) ! point-wise transformation
                            end do
            call disp%show("transpose(affinity)")
            call disp%show( transpose(affinity) )
            call disp%show("call setAffinity(affinInv, affinity, dim, getMatInv(tform), genrecmat) ! inverse transformation.")
                            call setAffinity(affinInv, affinity, dim, getMatInv(tform), genrecmat) ! inverse transformation.
            call disp%show("transpose(affinInv)")
            call disp%show( transpose(affinInv) )
            call disp%show("transpose(sample) ! for comparison with affinInv.")
            call disp%show( transpose(sample) )
            call disp%skip()
            call disp%show("call setMatInit(tform(2:,1:ndim-1), lowDia, vlow = 0._RKG, vdia = 0._RKG) ! make tform an upper-diagonal matrix.")
                            call setMatInit(tform(2:,1:ndim-1), lowDia, vlow = 0._RKG, vdia = 0._RKG) ! make tform an upper-diagonal matrix.
            call disp%show("tform")
            call disp%show( tform )
            call disp%show("call setAffinity(affinity, sample, dim, tform, class = upperDiag) ! whole sample transformation.")
                            call setAffinity(affinity, sample, dim, tform, class = upperDiag) ! whole sample transformation.
            call disp%show("transpose(affinity)")
            call disp%show( transpose(affinity) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    call setAffinity(affinity(:,isam), sample(:,isam), tform, upperDiag) ! point-wise transformation")
            call disp%show("end do")
                            do isam = 1, nsam
                                call setAffinity(affinity(:,isam), sample(:,isam), tform, upperDiag) ! point-wise transformation
                            end do
            call disp%show("transpose(affinity)")
            call disp%show( transpose(affinity) )
            call disp%show("call setAffinity(affinInv, affinity, dim, getMatInv(tform, upperDiag), class = upperDiag) ! inverse transformation.")
                            call setAffinity(affinInv, affinity, dim, getMatInv(tform, upperDiag), class = upperDiag) ! inverse transformation.
            call disp%show("transpose(affinInv)")
            call disp%show( transpose(affinInv) )
            call disp%show("transpose(sample) ! for comparison with affinInv.")
            call disp%show( transpose(sample) )
            call disp%skip()
            call disp%show("call setAffinity(affinity, sample, dim, tform, class = upperUnit) ! whole sample transformation.")
                            call setAffinity(affinity, sample, dim, tform, class = upperUnit) ! whole sample transformation.
            call disp%show("transpose(affinity)")
            call disp%show( transpose(affinity) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    call setAffinity(affinity(:,isam), sample(:,isam), tform, upperUnit) ! point-wise transformation")
            call disp%show("end do")
                            do isam = 1, nsam
                                call setAffinity(affinity(:,isam), sample(:,isam), tform, upperUnit) ! point-wise transformation
                            end do
            call disp%show("transpose(affinity)")
            call disp%show( transpose(affinity) )
            call disp%show("call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = upperUnit) ! inverse transformation.")
                            call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = upperUnit) ! inverse transformation.
            call disp%show("transpose(affinInv)")
            call disp%show( transpose(affinInv) )
            call disp%show("transpose(sample) ! for comparison with affinInv.")
            call disp%show( transpose(sample) )
            call disp%skip()
            call disp%show("call setAffinity(affinity, sample, dim, tform, class = lowerDiag) ! whole sample transformation.")
                            call setAffinity(affinity, sample, dim, tform, class = lowerDiag) ! whole sample transformation.
            call disp%show("transpose(affinity)")
            call disp%show( transpose(affinity) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    call setAffinity(affinity(:,isam), sample(:,isam), tform, lowerDiag) ! point-wise transformation")
            call disp%show("end do")
                            do isam = 1, nsam
                                call setAffinity(affinity(:,isam), sample(:,isam), tform, lowerDiag) ! point-wise transformation
                            end do
            call disp%show("transpose(affinity)")
            call disp%show( transpose(affinity) )
            call disp%show("call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = lowerDiag) ! inverse transformation.")
                            call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = lowerDiag) ! inverse transformation.
            call disp%show("transpose(affinInv)")
            call disp%show( transpose(affinInv) )
            call disp%show("transpose(sample) ! for comparison with affinInv.")
            call disp%show( transpose(sample) )
            call disp%skip()
            call disp%show("call setAffinity(affinity, sample, dim, tform, class = lowerUnit) ! whole sample transformation.")
                            call setAffinity(affinity, sample, dim, tform, class = lowerUnit) ! whole sample transformation.
            call disp%show("transpose(affinity)")
            call disp%show( transpose(affinity) )
            call disp%show("do isam = 1, nsam")
            call disp%show("    call setAffinity(affinity(:,isam), sample(:,isam), tform, lowerUnit) ! point-wise transformation")
            call disp%show("end do")
                            do isam = 1, nsam
                                call setAffinity(affinity(:,isam), sample(:,isam), tform, lowerUnit) ! point-wise transformation
                            end do
            call disp%show("transpose(affinity)")
            call disp%show( transpose(affinity) )
            call disp%show("call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = lowerUnit) ! inverse transformation.")
                            call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = lowerUnit) ! inverse transformation.
            call disp%show("transpose(affinInv)")
            call disp%show( transpose(affinInv) )
            call disp%show("transpose(sample) ! for comparison with affinInv.")
            call disp%show( transpose(sample) )
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Transform 2D sample along the first dimension.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKG => RKS
        real(RKG), allocatable :: tlate(:), tform(:,:), sample(:,:), affinInv(:,:), affinity(:,:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("dim = 1; ndim = getUnifRand(1, 3); nsam = getUnifRand(4, 5)")
                            dim = 1; ndim = getUnifRand(1, 3); nsam = getUnifRand(4, 5)
            call disp%show("[dim, ndim, nsam]")
            call disp%show( [dim, ndim, nsam] )
            call disp%show("tlate = getUnifRand(-5, +5, ndim)")
                            tlate = getUnifRand(-5, +5, ndim)
            call disp%show("tlate")
            call disp%show( tlate )
            call disp%show("sample = getUnifRand(-9, +9, nsam, ndim)")
                            sample = getUnifRand(-9, +9, nsam, ndim)
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("call setResized(affinity, shape(sample, IK))")
                            call setResized(affinity, shape(sample, IK))
            call disp%show("call setResized(affinInv, shape(sample, IK))")
                            call setResized(affinInv, shape(sample, IK))
            call disp%show("tform = getCovRand(1., ndim)")
                            tform = getCovRand(1., ndim)
            call disp%show("tform")
            call disp%show( tform )
            call disp%show("call setAffinity(affinity, sample, dim, tform) ! whole sample transformation.")
                            call setAffinity(affinity, sample, dim, tform, genrecmat) ! whole sample transformation.
            call disp%show("affinity")
            call disp%show( affinity )
            call disp%show("do isam = 1, nsam")
            call disp%show("    call setAffinity(affinity(isam,:), sample(isam,:), tform, genrecmat) ! point-wise transformation")
            call disp%show("end do")
                            do isam = 1, nsam
                                call setAffinity(affinity(isam,:), sample(isam,:), tform, genrecmat) ! point-wise transformation
                            end do
            call disp%show("affinity")
            call disp%show( affinity )
            call disp%show("call setAffinity(affinInv, affinity, dim, getMatInv(tform), genrecmat) ! inverse transformation.")
                            call setAffinity(affinInv, affinity, dim, getMatInv(tform), genrecmat) ! inverse transformation.
            call disp%show("affinInv")
            call disp%show( affinInv )
            call disp%show("sample ! for comparison with affinInv.")
            call disp%show( sample )
            call disp%skip()
            call disp%show("call setMatInit(tform(2:,1:ndim-1), lowDia, vlow = 0._RKG, vdia = 0._RKG) ! make tform an upper-diagonal matrix.")
                            call setMatInit(tform(2:,1:ndim-1), lowDia, vlow = 0._RKG, vdia = 0._RKG) ! make tform an upper-diagonal matrix.
            call disp%show("tform")
            call disp%show( tform )
            call disp%show("call setAffinity(affinity, sample, dim, tform, class = upperDiag) ! whole sample transformation.")
                            call setAffinity(affinity, sample, dim, tform, class = upperDiag) ! whole sample transformation.
            call disp%show("affinity")
            call disp%show( affinity )
            call disp%show("do isam = 1, nsam")
            call disp%show("    call setAffinity(affinity(isam,:), sample(isam,:), tform, upperDiag) ! point-wise transformation")
            call disp%show("end do")
                            do isam = 1, nsam
                                call setAffinity(affinity(isam,:), sample(isam,:), tform, upperDiag) ! point-wise transformation
                            end do
            call disp%show("affinity")
            call disp%show( affinity )
            call disp%show("call setAffinity(affinInv, affinity, dim, getMatInv(tform, upperDiag), class = upperDiag) ! inverse transformation.")
                            call setAffinity(affinInv, affinity, dim, getMatInv(tform, upperDiag), class = upperDiag) ! inverse transformation.
            call disp%show("affinInv")
            call disp%show( affinInv )
            call disp%show("sample ! for comparison with affinInv.")
            call disp%show( sample )
            call disp%skip()
            call disp%show("call setAffinity(affinity, sample, dim, tform, class = upperUnit) ! whole sample transformation.")
                            call setAffinity(affinity, sample, dim, tform, class = upperUnit) ! whole sample transformation.
            call disp%show("affinity")
            call disp%show( affinity )
            call disp%show("do isam = 1, nsam")
            call disp%show("    call setAffinity(affinity(isam,:), sample(isam,:), tform, upperUnit) ! point-wise transformation")
            call disp%show("end do")
                            do isam = 1, nsam
                                call setAffinity(affinity(isam,:), sample(isam,:), tform, upperUnit) ! point-wise transformation
                            end do
            call disp%show("affinity")
            call disp%show( affinity )
            call disp%show("call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = upperUnit) ! inverse transformation.")
                            call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = upperUnit) ! inverse transformation.
            call disp%show("affinInv")
            call disp%show( affinInv )
            call disp%show("sample ! for comparison with affinInv.")
            call disp%show( sample )
            call disp%skip()
            call disp%show("call setAffinity(affinity, sample, dim, tform, class = lowerDiag) ! whole sample transformation.")
                            call setAffinity(affinity, sample, dim, tform, class = lowerDiag) ! whole sample transformation.
            call disp%show("affinity")
            call disp%show( affinity )
            call disp%show("do isam = 1, nsam")
            call disp%show("    call setAffinity(affinity(isam,:), sample(isam,:), tform, lowerDiag) ! point-wise transformation")
            call disp%show("end do")
                            do isam = 1, nsam
                                call setAffinity(affinity(isam,:), sample(isam,:), tform, lowerDiag) ! point-wise transformation
                            end do
            call disp%show("affinity")
            call disp%show( affinity )
            call disp%show("call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = lowerDiag) ! inverse transformation.")
                            call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = lowerDiag) ! inverse transformation.
            call disp%show("affinInv")
            call disp%show( affinInv )
            call disp%show("sample ! for comparison with affinInv.")
            call disp%show( sample )
            call disp%skip()
            call disp%show("call setAffinity(affinity, sample, dim, tform, class = lowerUnit) ! whole sample transformation.")
                            call setAffinity(affinity, sample, dim, tform, class = lowerUnit) ! whole sample transformation.
            call disp%show("affinity")
            call disp%show( affinity )
            call disp%show("do isam = 1, nsam")
            call disp%show("    call setAffinity(affinity(isam,:), sample(isam,:), tform, lowerUnit) ! point-wise transformation")
            call disp%show("end do")
                            do isam = 1, nsam
                                call setAffinity(affinity(isam,:), sample(isam,:), tform, lowerUnit) ! point-wise transformation
                            end do
            call disp%show("affinity")
            call disp%show( affinity )
            call disp%show("call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = lowerUnit) ! inverse transformation.")
                            call setAffinity(affinInv, affinity, dim, getMatInv(tform), class = lowerUnit) ! inverse transformation.
            call disp%show("affinInv")
            call disp%show( affinInv )
            call disp%show("sample ! for comparison with affinInv.")
            call disp%show( sample )
        end do
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_mathConst, only: PI
        use pm_container, only: css_type
        use pm_matrixInit, only: uppLowDia
        use pm_matrixInit, only: getMatInit
        use pm_arraySpace, only: getLinSpace
        use pm_sampleAffinity, only: genrecmat
        use pm_io, only: getErrTableWrite, trans
        use pm_distUnifEll, only: getUnifEllRand, uppDia
        use pm_sampleShift, only: getShifted
        use pm_matrixInv, only: getMatInv
        use pm_kind, only: RKG => RKS
        type(css_type), allocatable :: shapes(:), tforms(:)
        integer(IK) :: ishape, itform
        integer(IK), parameter :: nsam = 1000, ndim = 2, dim = 2
        real(RKG), allocatable :: angle, tlate(:), cov(:,:), tform(:,:), sample(:,:), affinity(:,:), affinInv(:,:)
        tforms = [css_type(SK_"warp"), css_type(SK_"rotation")]
        shapes = [css_type(SK_"circle"), css_type(SK_"square")]
        tlate = getUnifRand(0, 5, ndim)
        do ishape = 1, size(shapes)
            if (shapes(ishape)%val == "circle") then
                sample = getUnifEllRand(getMatInit([ndim, ndim], uppLowDia, 0., 0., getLinSpace(1., real(ndim), ndim)), uppDia, nsam)
            elseif (shapes(ishape)%val == "square") then
                sample = getUnifRand(0., 1., ndim, nsam)
            else
                error stop "Unrecognized shape."
            end if
            call setResized(affinity, shape(sample, IK))
            call setResized(affinInv, shape(sample, IK))
            do itform = 1, size(tforms)
                if (tforms(itform)%val == "rotation") then
                    angle = PI / 4
                    tform = reshape([cos(angle), -sin(angle), sin(angle), cos(angle)], [ndim, ndim])
                    call setAffinity(affinity, sample, dim, tform, genrecmat, tlate)
                elseif (tforms(itform)%val == "warp") then
                    cov = reshape([1., .5, .5, 1.], shape = [ndim, ndim]) ! getCovRand(1., ndim)
                    tform = getMatChol(cov, lowDia)
                    call setAffinity(affinity, sample, dim, tform, lowerDiag, tlate)
                else
                    error stop "Unrecognized tform."
                end if
                call setAffinity(affinInv, getShifted(affinity, dim, -tlate), dim, getMatInv(tform), genrecmat)
                if (0 /= getErrTableWrite(SK_"setAffinity."//tforms(itform)%val//SK_"."//shapes(ishape)%val//SK_".sample.txt", sample, trans)) error stop "Failed table-write."
                if (0 /= getErrTableWrite(SK_"setAffinity."//tforms(itform)%val//SK_"."//shapes(ishape)%val//SK_".affinity.txt", affinity, trans)) error stop "Failed table-write."
                if (0 /= getErrTableWrite(SK_"setAffinity."//tforms(itform)%val//SK_"."//shapes(ishape)%val//SK_".affinInv.txt", affinInv, trans)) error stop "Failed table-write."
            end do
        end do
    end block

end program example