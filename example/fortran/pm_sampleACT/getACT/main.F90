program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: getErrTableRead, TAB
    use pm_sampleACT, only: batchMeansMax, batchMeansMax_type
    use pm_sampleACT, only: batchMeans, batchMeans_type
    use pm_sampleACT, only: cumSumMax, cumSumMax_type
    use pm_sampleACT, only: cumSum, cumSum_type
    use pm_sampleACT, only: getACT
    use pm_io, only: display_type

    implicit none

    integer(IK) :: skip
    type(display_type) :: disp
    character(:), allocatable :: format
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the cross-correlation of two samples.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported.
        real(TKC), allocatable :: table(:,:)
        real(TKC) :: act
        call disp%skip()
        call disp%show("if (0 /= getErrTableRead(SK_'duluth.txt', table, sep = TAB, roff = 1_IK)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(SK_'duluth.txt', table, sep = TAB, roff = 1_IK)) error stop 'table read failed.'
        call disp%show("act = getACT(table(:,2)) ! temperature see visualization below.")
                        act = getACT(table(:,2)) ! temperature see visualization below.
        call disp%show("act")
        call disp%show( act )
        call disp%show("act = getACT(table(:,2), method = batchMeans)")
                        act = getACT(table(:,2), method = batchMeans)
        call disp%show("act")
        call disp%show( act )
        call disp%show("act = getACT(table(:,2), method = batchMeans_type(size = 365_IK))")
                        act = getACT(table(:,2), method = batchMeans_type(size = 365_IK))
        call disp%show("act")
        call disp%show( act )
        call disp%show("act = getACT(table(:,2), method = batchMeansMax_type())")
                        act = getACT(table(:,2), method = batchMeansMax_type())
                        skip = int(act, IK)
        call disp%show("act")
        call disp%show( act )
        call disp%show("act = getACT(table(:,2), method = batchMeansMax_type(step = 10_IK))")
                        act = getACT(table(:,2), method = batchMeansMax_type(step = 10_IK))
        call disp%show("act")
        call disp%show( act )
        call disp%show("act = getACT(table(:,2), method = cumSum)")
                        act = getACT(table(:,2), method = cumSum)
        call disp%show("act")
        call disp%show( act )
        call disp%show("act = getACT(table(:,2), method = cumSum_type(signif = 5_IK))")
                        act = getACT(table(:,2), method = cumSum_type(signif = 5_IK))
        call disp%show("act")
        call disp%show( act )
        call disp%show("act = getACT(table(:,2), method = cumSum_type(signif = 0_IK))")
                        act = getACT(table(:,2), method = cumSum_type(signif = 0_IK))
        call disp%show("act")
        call disp%show( act )
        call disp%show("act = getACT(table(:,2), method = cumSumMax)")
                        act = getACT(table(:,2), method = cumSumMax)
        call disp%show("act")
        call disp%show( act )
        call disp%skip()
    end block

    block
        use pm_arrayResize, only: setResized
        use pm_io, only: getErrTableWrite, trans
        use pm_arrayRange, only: getRange
        use pm_kind, only: TKC => RKS ! All other real types are also supported.
        character(:, SK), allocatable :: header
        integer(IK), allocatable :: sbatch(:)
        real(TKC), allocatable :: output(:,:)
        real(TKC), allocatable :: table(:,:)
        integer(IK) :: ibatch
        if (0 /= getErrTableRead(SK_'duluth.txt', table, header, sep = TAB)) error stop 'table read failed.'
        sbatch = getRange(2_IK, size(table, 1, IK) / 2_IK)
        call setResized(output, [size(sbatch, 1, IK), 2_IK])
        do ibatch = 1, size(sbatch, 1, IK)
            output(ibatch, 1) = sbatch(ibatch)
            output(ibatch, 2) = getACT(table(:, 2), batchMeans_type(sbatch(ibatch)))
        end do
        if (0 /= getErrTableWrite(SK_'getACT.duluth.batchMeans.txt', output, header = SK_"batchSize"//TAB//SK_"ACT", sep = TAB)) error stop 'table read failed.'
        ! Write the thinned sample.
        table = table(1 : size(table, 1, IK) : skip, :)
        if (0 /= getErrTableWrite(SK_'getACT.duluth.thinned.txt', table, header//SK_" (BMM-thinned)", sep = TAB)) error stop 'table read failed.'
    end block

end program example