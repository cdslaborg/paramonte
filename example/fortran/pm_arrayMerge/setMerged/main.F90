program example

    use pm_kind, only: SK, IK
    use pm_arraySort, only: getSorted
    use pm_distUnif, only: getUnifRand
    use pm_arrayMerge, only: setMerged
    use pm_arrayResize, only: setResized
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    integer(IK) :: fileUnit, i
    integer(IK) , parameter :: NP1 = 2_IK, NP2 = 3_IK

    disp = display_type(file = SK_"main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("!character merger.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => SK ! All kinds are supported.
        character(:,TKG), allocatable :: sortedArray1, sortedArray2, sortedMerger
        integer(IK) :: len1, len2
        call disp%skip()
        call disp%show("len1 = getUnifRand(0_IK, 5_IK); len2 = getUnifRand(0_IK, 5_IK)")
                        len1 = getUnifRand(0_IK, 5_IK); len2 = getUnifRand(0_IK, 5_IK)
        call disp%show("[len1, len2]")
        call disp%show( [len1, len2] )
        call disp%show("sortedArray1 = getSorted(getUnifRand(repeat('A', len1), repeat('Z', len1)))")
                        sortedArray1 = getSorted(getUnifRand(repeat('A', len1), repeat('Z', len1)))
        call disp%show("sortedArray1")
        call disp%show( sortedArray1 , deliml = SK_"""" )
        call disp%show("sortedArray2 = getSorted(getUnifRand(repeat('A', len2), repeat('Z', len2)))")
                        sortedArray2 = getSorted(getUnifRand(repeat('A', len2), repeat('Z', len2)))
        call disp%show("sortedArray2")
        call disp%show( sortedArray2 , deliml = SK_"""" )
        call disp%show("call setResized(sortedMerger, len1 + len2)")
                        call setResized(sortedMerger, len1 + len2)
        call disp%show("call setMerged(sortedMerger, sortedArray1, sortedArray2)")
                        call setMerged(sortedMerger, sortedArray1, sortedArray2)
        call disp%show("sortedMerger")
        call disp%show( sortedMerger , deliml = SK_"""" )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!integer merger.")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => IK ! All kinds are supported.
        integer(TKG), allocatable :: sortedArray1(:), sortedArray2(:), sortedMerger(:)
        integer(IK) :: len1, len2
        call disp%skip()
        call disp%show("len1 = getUnifRand(0_IK, 5_IK); len2 = getUnifRand(0_IK, 5_IK)")
                        len1 = getUnifRand(0_IK, 5_IK); len2 = getUnifRand(0_IK, 5_IK)
        call disp%show("[len1, len2]")
        call disp%show( [len1, len2] )
        call disp%show("sortedArray1 = getSorted(getUnifRand(0, 9, len1))")
                        sortedArray1 = getSorted(getUnifRand(0, 9, len1))
        call disp%show("sortedArray1")
        call disp%show( sortedArray1 , deliml = SK_"""" )
        call disp%show("sortedArray2 = getSorted(getUnifRand(0, 9, len2))")
                        sortedArray2 = getSorted(getUnifRand(0, 9, len2))
        call disp%show("sortedArray2")
        call disp%show( sortedArray2 , deliml = SK_"""" )
        call disp%show("call setResized(sortedMerger, len1 + len2)")
                        call setResized(sortedMerger, len1 + len2)
        call disp%show("call setMerged(sortedMerger, sortedArray1, sortedArray2)")
                        call setMerged(sortedMerger, sortedArray1, sortedArray2)
        call disp%show("sortedMerger")
        call disp%show( sortedMerger , deliml = SK_"""" )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%")
    call disp%show("!real merger.")
    call disp%show("!%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => IK ! All kinds are supported.
        real(TKG), allocatable :: sortedArray1(:), sortedArray2(:), sortedMerger(:)
        integer(IK) :: len1, len2
        call disp%skip()
        call disp%show("len1 = getUnifRand(0_IK, 5_IK); len2 = getUnifRand(0_IK, 5_IK)")
                        len1 = getUnifRand(0_IK, 5_IK); len2 = getUnifRand(0_IK, 5_IK)
        call disp%show("[len1, len2]")
        call disp%show( [len1, len2] )
        call disp%show("sortedArray1 = getSorted(getUnifRand(0, 9, len1))")
                        sortedArray1 = getSorted(getUnifRand(0, 9, len1))
        call disp%show("sortedArray1")
        call disp%show( sortedArray1 , deliml = SK_"""" )
        call disp%show("sortedArray2 = getSorted(getUnifRand(0, 9, len2))")
                        sortedArray2 = getSorted(getUnifRand(0, 9, len2))
        call disp%show("sortedArray2")
        call disp%show( sortedArray2 , deliml = SK_"""" )
        call disp%show("call setResized(sortedMerger, len1 + len2)")
                        call setResized(sortedMerger, len1 + len2)
        call disp%show("call setMerged(sortedMerger, sortedArray1, sortedArray2)")
                        call setMerged(sortedMerger, sortedArray1, sortedArray2)
        call disp%show("sortedMerger")
        call disp%show( sortedMerger , deliml = SK_"""" )
        call disp%skip()
    end block

end program example