program example

    use pm_kind, only: IK, RK
    use pm_kind, only: RKS, RKD, RKH
    use pm_kind, only: IK8, IK16, IK32, IK64
    use pm_arraySpace, only: getLinSpace
    use pm_arrayMerge, only: getMerged

    implicit none

    integer(IK), parameter :: NP1 = 2_IK, NP2 = 3_IK

    ! 1-dimensional array of real values.

    real(RKS)   :: SortedArray1_RKS(NP1), SortedArray2_RKS(NP2), MergedArray_RKS(NP1+NP2)
    real(RKD)   :: SortedArray1_RKD(NP1), SortedArray2_RKD(NP2), MergedArray_RKD(NP1+NP2)
    real(RKH)   :: SortedArray1_RKH(NP1), SortedArray2_RKH(NP2), MergedArray_RKH(NP1+NP2)

    ! 1-dimensional array of integer values.

    integer(IK8)   :: SortedArray1_IK8  (NP1), SortedArray2_IK8  (NP2), MergedArray_IK8  (NP1+NP2)
    integer(IK16)   :: SortedArray1_IK16 (NP1), SortedArray2_IK16 (NP2), MergedArray_IK16 (NP1+NP2)
    integer(IK32)   :: SortedArray1_IK32 (NP1), SortedArray2_IK32 (NP2), MergedArray_IK32 (NP1+NP2)
    integer(IK64)   :: SortedArray1_IK64 (NP1), SortedArray2_IK64 (NP2), MergedArray_IK64 (NP1+NP2)

    integer(IK)     :: fileUnit, i

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Define the sorted arrays.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SortedArray1_RKS  = getLinSpace(0._RK, 4._RK, NP1)
    SortedArray1_RKD  = getLinSpace(0._RK, 4._RK, NP1)
    SortedArray1_RKH = getLinSpace(0._RK, 4._RK, NP1)

    SortedArray2_RKS  = getLinSpace(1._RK, 5._RK, NP2)
    SortedArray2_RKD  = getLinSpace(1._RK, 5._RK, NP2)
    SortedArray2_RKH = getLinSpace(1._RK, 5._RK, NP2)

    SortedArray1_IK8   = int(SortedArray1_RKS, kind = IK8)
    SortedArray1_IK16  = int(SortedArray1_RKS, kind = IK16)
    SortedArray1_IK32  = int(SortedArray1_RKS, kind = IK32)
    SortedArray1_IK64  = int(SortedArray1_RKS, kind = IK64)

    SortedArray2_IK8   = int(SortedArray2_RKS, kind = IK8)
    SortedArray2_IK16  = int(SortedArray2_RKS, kind = IK16)
    SortedArray2_IK32  = int(SortedArray2_RKS, kind = IK32)
    SortedArray2_IK64  = int(SortedArray2_RKS, kind = IK64)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Sort-merge arrays in ascending order.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    MergedArray_IK8     = getMerged(SortedArray1_IK8, SortedArray2_IK8)
    MergedArray_IK16    = getMerged(SortedArray1_IK16, SortedArray2_IK16)
    MergedArray_IK32    = getMerged(SortedArray1_IK32, SortedArray2_IK32)
    MergedArray_IK64    = getMerged(SortedArray1_IK64, SortedArray2_IK64)

    MergedArray_RKS    = getMerged(SortedArray1_RKS, SortedArray2_RKS)
    MergedArray_RKD    = getMerged(SortedArray1_RKD, SortedArray2_RKD)
    MergedArray_RKH   = getMerged(SortedArray1_RKH, SortedArray2_RKH)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Write arrays to the output file.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    open(newunit = fileUnit, file = "main.out.F90")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") "  SortedArray1_IK8",  SortedArray1_IK8
    write(fileUnit,"(*(g0,:,', '))") "  SortedArray2_IK8",  SortedArray2_IK8
    write(fileUnit,"(*(g0,:,', '))") "   MergedArray_IK8",   MergedArray_IK8
    write(fileUnit,"(*(g0,:,', '))")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") " SortedArray1_IK16",  SortedArray1_IK16
    write(fileUnit,"(*(g0,:,', '))") " SortedArray2_IK16",  SortedArray2_IK16
    write(fileUnit,"(*(g0,:,', '))") "  MergedArray_IK16",   MergedArray_IK16
    write(fileUnit,"(*(g0,:,', '))")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") " SortedArray1_IK32",  SortedArray1_IK32
    write(fileUnit,"(*(g0,:,', '))") " SortedArray2_IK32",  SortedArray2_IK32
    write(fileUnit,"(*(g0,:,', '))") "  MergedArray_IK32",   MergedArray_IK32
    write(fileUnit,"(*(g0,:,', '))")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") " SortedArray1_IK64",  SortedArray1_IK64
    write(fileUnit,"(*(g0,:,', '))") " SortedArray2_IK64",  SortedArray2_IK64
    write(fileUnit,"(*(g0,:,', '))") "  MergedArray_IK64",   MergedArray_IK64
    write(fileUnit,"(*(g0,:,', '))")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") " SortedArray1_RKS",  SortedArray1_RKS
    write(fileUnit,"(*(g0,:,', '))") " SortedArray2_RKS",  SortedArray2_RKS
    write(fileUnit,"(*(g0,:,', '))") "  MergedArray_RKS",   MergedArray_RKS
    write(fileUnit,"(*(g0,:,', '))")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") " SortedArray1_RKD",  SortedArray1_RKD
    write(fileUnit,"(*(g0,:,', '))") " SortedArray2_RKD",  SortedArray2_RKD
    write(fileUnit,"(*(g0,:,', '))") "  MergedArray_RKD",   MergedArray_RKD
    write(fileUnit,"(*(g0,:,', '))")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") "SortedArray1_RKH", SortedArray1_RKH
    write(fileUnit,"(*(g0,:,', '))") "SortedArray2_RKH", SortedArray2_RKH
    write(fileUnit,"(*(g0,:,', '))") " MergedArray_RKH",  MergedArray_RKH
    write(fileUnit,"(*(g0,:,', '))")

    close(fileUnit)

end program example