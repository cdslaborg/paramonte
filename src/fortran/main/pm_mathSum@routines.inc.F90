!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This include file contains procedure implementations of [pm_mathSum](@ref pm_mathSum).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, August 8, 2024, 10:23 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define GET_CONJG(X)conjg(X)
#define TYPE_KIND(X)complex(X)
#elif   RK_ENABLED
#define GET_CONJG(X)X
#define TYPE_KIND(X)real(X)
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%
#if     getSum_ENABLED
        !%%%%%%%%%%%%%

        TYPE_KIND(TKG), parameter :: ONE = 1._TKG, ZERO = 0._TKG
#if     Def_ENABLED || FAB_ENABLED
        integer(IK), parameter :: blockSize = 128_IK
        integer(IK) :: ipart, npart, lenx, ibeg, iend
        real(TKG), parameter :: blockSizeInv = 1._TKG / real(blockSize, TKG)
        TYPE_KIND(TKG) :: temp, ydum, error
        lenx = size(x, 1, IK)
        npart = ceiling(lenx * blockSizeInv, IK)
        sumres = ZERO
        error = ZERO
        do ipart = 1_IK, npart
            ibeg = 1_IK + (ipart - 1_IK) * blockSize
            iend = min(ipart * blockSize, lenx)
            temp = sumres
            ydum = error + sum(x(ibeg : iend))
            sumres = temp + ydum
            error = (temp - sumres) + ydum
        end do
#elif   NAB_ENABLED
        integer(IK), parameter :: blockSize = 128_IK
        real(TKG), parameter :: blockSizeInv = 1._TKG / real(blockSize, TKG)
        integer(IK) :: ipart, npart, lenx, ibeg, iend
        lenx = size(x, 1, IK)
        npart = ceiling(lenx * blockSizeInv, IK)
        sumres = ZERO
        do ipart = 1_IK, npart
            ibeg = 1_IK + (ipart - 1_IK) * blockSize
            iend = min(ipart * blockSize, lenx)
            sumres = sumres + sum(x(ibeg : iend))
        end do
#elif   KAB_ENABLED
        TYPE_KIND(TKG) :: temp, ydum, error
        integer(IK) :: i
        sumres = ZERO
        error = ZERO
        do i = 1_IK, size(x, 1, IK)
            temp = sumres
            ydum = error + x(i)
            sumres = temp + ydum
            error = (temp - sumres) + ydum
        end do
#elif   Ite_ENABLED
        integer(IK) :: i
        sumres = ZERO
        do i = 1_IK, size(x, 1, IK)
            sumres = sumres + x(i)
        end do
#elif   Rec_ENABLED
        integer(IK), parameter :: blockSize = 128_IK
        integer(IK) :: lenx, half
        lenx = size(x, 1, IK)
        if (lenx < blockSize) then
            sumres = sum(x)
        else
            half = lenx / 2_IK
            sumres = getSum(x(1 : half), recursion) + getSum(x(half + 1 : lenx), recursion)
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%
#elif   getDot_ENABLED
        !%%%%%%%%%%%%%

        TYPE_KIND(TKG), parameter :: ONE = 1._TKG, ZERO = 0._TKG
#if     Def_ENABLED || FAB_ENABLED
        integer(IK), parameter :: blockSize = 128_IK
        integer(IK) :: ipart, npart, lenx, ibeg, iend
        real(TKG), parameter :: blockSizeInv = 1._TKG / real(blockSize, TKG)
        TYPE_KIND(TKG) :: temp, ydum, error
        lenx = size(x, 1, IK)
        npart = ceiling(lenx * blockSizeInv, IK)
        dotres = ZERO
        error = ZERO
        do ipart = 1_IK, npart
            ibeg = 1_IK + (ipart - 1_IK) * blockSize
            iend = min(ipart * blockSize, lenx)
            temp = dotres
            ydum = error + dot_product(x(ibeg : iend), y(ibeg : iend))
            dotres = temp + ydum
            error = (temp - dotres) + ydum
        end do
#elif   NAB_ENABLED
        integer(IK), parameter :: blockSize = 128_IK
        real(TKG), parameter :: blockSizeInv = 1._TKG / real(blockSize, TKG)
        integer(IK) :: ipart, npart, lenx, ibeg, iend
        lenx = size(x, 1, IK)
        npart = ceiling(lenx * blockSizeInv, IK)
        dotres = ZERO
        do ipart = 1_IK, npart
            ibeg = 1_IK + (ipart - 1_IK) * blockSize
            iend = min(ipart * blockSize, lenx)
            dotres = dotres + dot_product(x(ibeg : iend), y(ibeg : iend))
        end do
#elif   KAB_ENABLED
        TYPE_KIND(TKG) :: temp, ydum, error
        integer(IK) :: i
        dotres = ZERO
        error = ZERO
        do i = 1_IK, size(x, 1, IK)
            temp = dotres
            ydum = error + GET_CONJG(x(i)) * y(i)
            dotres = temp + ydum
            error = (temp - dotres) + ydum
        end do
#elif   Ite_ENABLED
        integer(IK) :: i
        dotres = ZERO
        do i = 1_IK, size(x, 1, IK)
            dotres = dotres + GET_CONJG(x(i)) * y(i)
        end do
#elif   Rec_ENABLED
        integer(IK), parameter :: blockSize = 128_IK
        integer(IK) :: lenx, half
        lenx = size(x, 1, IK)
        if (lenx < blockSize) then
            dotres = dot_product(x, y)
        else
            half = lenx / 2_IK
            dotres = getDot(x(1 : half), y(1 : half), recursion) + getDot(x(half + 1 : lenx), y(half + 1 : lenx), recursion)
        end if
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, size(x, 1, IK) == size(y, 1, IK), SK_"@getDot(): The condition `size(x) == size(y)` must hold for corresponding input arguments. size(x), size(y) = "//getStr([size(x, 1, IK), size(y, 1, IK)]))
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  TYPE_KIND
#undef  GET_CONJG