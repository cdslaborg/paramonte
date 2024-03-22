program example

    use pm_kind, only: IK, RK
    use pm_kind, only: RKS, RKD, RKH
    use pm_kind, only: IKL, IKS, IKD, IKH
    use pm_arraySpace, only: getLinSpace
    use pm_arrayMerge, only: getMerged

    implicit none

    integer(IK), parameter :: NP1 = 2_IK, NP2 = 3_IK

    ! 1-dimensional array of real values.

    real(RKS)   :: SortedArray1_RKS(NP1), SortedArray2_RKS(NP2), MergedArray_RKS(NP1+NP2)
    real(RKD)   :: SortedArray1_RKD(NP1), SortedArray2_RKD(NP2), MergedArray_RKD(NP1+NP2)
    real(RKH)   :: SortedArray1_RKH(NP1), SortedArray2_RKH(NP2), MergedArray_RKH(NP1+NP2)

    ! 1-dimensional array of integer values.

    integer(IKL)    :: SortedArray1_IKL(NP1), SortedArray2_IKL(NP2), MergedArray_IKL(NP1+NP2)
    integer(IKS)    :: SortedArray1_IKS(NP1), SortedArray2_IKS(NP2), MergedArray_IKS(NP1+NP2)
    integer(IKD)    :: SortedArray1_IKD(NP1), SortedArray2_IKD(NP2), MergedArray_IKD(NP1+NP2)
    integer(IKH)    :: SortedArray1_IKH(NP1), SortedArray2_IKH(NP2), MergedArray_IKH(NP1+NP2)

    integer(IK)     :: fileUnit, i

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Define the sorted arrays.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SortedArray1_RKS = getLinSpace(0._RK, 4._RK, NP1)
    SortedArray1_RKD = getLinSpace(0._RK, 4._RK, NP1)
    SortedArray1_RKH = getLinSpace(0._RK, 4._RK, NP1)

    SortedArray2_RKS = getLinSpace(1._RK, 5._RK, NP2)
    SortedArray2_RKD = getLinSpace(1._RK, 5._RK, NP2)
    SortedArray2_RKH = getLinSpace(1._RK, 5._RK, NP2)

    SortedArray1_IKL  = int(SortedArray1_RKS, kind = IKL)
    SortedArray1_IKS  = int(SortedArray1_RKS, kind = IKS)
    SortedArray1_IKD  = int(SortedArray1_RKS, kind = IKD)
    SortedArray1_IKH  = int(SortedArray1_RKS, kind = IKH)

    SortedArray2_IKL  = int(SortedArray2_RKS, kind = IKL)
    SortedArray2_IKS  = int(SortedArray2_RKS, kind = IKS)
    SortedArray2_IKD  = int(SortedArray2_RKS, kind = IKD)
    SortedArray2_IKH  = int(SortedArray2_RKS, kind = IKH)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Sort-merge arrays in ascending order.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    MergedArray_IKL = getMerged(SortedArray1_IKL, SortedArray2_IKL)
    MergedArray_IKS = getMerged(SortedArray1_IKS, SortedArray2_IKS)
    MergedArray_IKD = getMerged(SortedArray1_IKD, SortedArray2_IKD)
    MergedArray_IKH = getMerged(SortedArray1_IKH, SortedArray2_IKH)

    MergedArray_RKS = getMerged(SortedArray1_RKS, SortedArray2_RKS)
    MergedArray_RKD = getMerged(SortedArray1_RKD, SortedArray2_RKD)
    MergedArray_RKH = getMerged(SortedArray1_RKH, SortedArray2_RKH)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Write arrays to the output file.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    open(newunit = fileUnit, file = "main.out.F90")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") "  SortedArray1_IKL",  SortedArray1_IKL
    write(fileUnit,"(*(g0,:,', '))") "  SortedArray2_IKL",  SortedArray2_IKL
    write(fileUnit,"(*(g0,:,', '))") "   MergedArray_IKL",   MergedArray_IKL
    write(fileUnit,"(*(g0,:,', '))")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") " SortedArray1_IKS",  SortedArray1_IKS
    write(fileUnit,"(*(g0,:,', '))") " SortedArray2_IKS",  SortedArray2_IKS
    write(fileUnit,"(*(g0,:,', '))") "  MergedArray_IKS",   MergedArray_IKS
    write(fileUnit,"(*(g0,:,', '))")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") " SortedArray1_IKD",  SortedArray1_IKD
    write(fileUnit,"(*(g0,:,', '))") " SortedArray2_IKD",  SortedArray2_IKD
    write(fileUnit,"(*(g0,:,', '))") "  MergedArray_IKD",   MergedArray_IKD
    write(fileUnit,"(*(g0,:,', '))")

    write(fileUnit,"(*(g0,:,', '))")
    write(fileUnit,"(*(g0,:,', '))") " SortedArray1_IKH",  SortedArray1_IKH
    write(fileUnit,"(*(g0,:,', '))") " SortedArray2_IKH",  SortedArray2_IKH
    write(fileUnit,"(*(g0,:,', '))") "  MergedArray_IKH",   MergedArray_IKH
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