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
!>  This include file contains procedure implementation of [pm_mathCumSum](@ref pm_mathCumSum).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 3:33 AM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%
#if     getCumSum_ENABLED
        !%%%%%%%%%%%%%%%%

        if (0_IK < size(array, 1, IK)) then
            if (present(direction) .and. present(action)) then
                if (same_type_as(direction, forward)) then
                    if (same_type_as(action, nothing)) then
                        call setCumSum(cumsum, array, forward, nothing)
                    elseif (same_type_as(action, reverse)) then
                        call setCumSum(cumsum, array, forward, reverse)
                    else
                        error stop "@getCumSum(): Unrecognized action."
                    end if
                elseif (same_type_as(direction, backward)) then
                    if (same_type_as(action, nothing)) then
                        call setCumSum(cumsum, array, backward, nothing)
                    elseif (same_type_as(action, reverse)) then
                        call setCumSum(cumsum, array, backward, reverse)
                    else
                        error stop "@getCumSum(): Unrecognized action."
                    end if
                else
                    error stop "@getCumSum(): Unrecognized direction."
                end if
            elseif (present(direction)) then
                if (same_type_as(direction, forward)) then
                    call setCumSum(cumsum, array, forward, nothing)
                elseif (same_type_as(direction, backward)) then
                    call setCumSum(cumsum, array, backward, nothing)
                else
                    error stop "@getCumSum(): Unrecognized direction."
                end if
            elseif (present(action)) then
                if (same_type_as(action, nothing)) then
                    call setCumSum(cumsum, array, forward, nothing)
                elseif (same_type_as(action, reverse)) then
                    call setCumSum(cumsum, array, forward, reverse)
                else
                    error stop "@getCumSum(): Unrecognized action."
                end if
            else
                call setCumSum(cumsum, array)
            end if
        end if

        !%%%%%%%%%%%%%%%%
#elif   setCumSum_ENABLED
        !%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define TYPE_KIND complex(TKG)
#elif   RK_ENABLED
#define TYPE_KIND real(TKG)
#elif   !IK_ENABLED
#error  "Unrecognized interface."
#endif

#if     New_ENABLED || (For_ENABLED && Non_ENABLED) || (Bac_ENABLED && Rev_ENABLED)
        integer(IK) :: i
#endif
        integer(IK) :: lenArray
        integer(IK), parameter :: BLOCKSIZE = 128_IK
        lenArray = size(array, kind = IK)
        CHECK_ASSERTION(__LINE__, 0_IK < lenArray, SK_"@setCumSum(): The condition `0 < size(array)` must hold. size(array) = "//getStr(lenArray))
#if     Old_ENABLED && For_ENABLED && Non_ENABLED && (IK_ENABLED || CK_ENABLED || RK_ENABLED)
        do i = 2, lenArray
            array(i) = array(i - 1) + array(i)
        end do
#elif   Old_ENABLED && For_ENABLED && Non_ENABLED && (CK_ENABLED || RK_ENABLED) && 0
        ! This recursive approach is about 4 and 5 times slower
        ! than the naive approach for 32 and 64 bit numbers.
        block
            integer(IK) :: lenArrayHalf
            if (lenArray < BLOCKSIZE) then
                do i = 2, lenArray
                    array(i) = array(i - 1) + array(i)
                end do
            else
                lenArrayHalf = lenArray / 2_IK
                call setCumSum(array(1 : lenArrayHalf))
                call setCumSum(array(lenArrayHalf + 1 : lenArray))
                array(lenArrayHalf + 1 : lenArray) = array(lenArrayHalf) + array(lenArrayHalf + 1 : lenArray)
            end if
        end block
        !block
        !    TYPE_KIND :: error, temp, y
        !    integer(IK) :: ipart, npart, ibeg, iend
        !    npart = (lenArray - 1) / BLOCKSIZE
        !    error = 0._TKG
        !    do ipart = 1, npart
        !        iend = min(ipart * BLOCKSIZE, lenArray)
        !        ibeg = 1 + (ipart - 1) * BLOCKSIZE
        !        temp = array(ibeg)
        !        array(ibeg + 1) = error
        !        do i = ibeg + 2, iend
        !            array(i) = array(i - 1) + array(i)
        !        end do
        !        y = array(iend)
        !        array(ibeg : iend) = temp + array(ibeg : iend)
        !        error = (temp - array(iend)) + y;
        !    end do
        !end block
#elif   Old_ENABLED && For_ENABLED && Non_ENABLED && (CK_ENABLED || RK_ENABLED) && 0
        ! This compensated summation algorithm is more accurate than naive
        ! summation for extremely large array sizes (> 10**7 elements),
        ! but is roughly 4X and 3X slower than the naive approach
        ! for 32 and 64 bit arrays. Hence, we use the naive approach.
        block
            TYPE_KIND :: error, temp, y
            error = 0._TKG
            do i = 2, lenArray
                temp = array(i - 1)
                y = error + array(i)
                array(i) = temp + y
                error = (temp - array(i)) + y
            end do
        end block
#elif   Old_ENABLED && For_ENABLED && Rev_ENABLED
        call setCumSum(array)
        call setReversed(array)
#elif   Old_ENABLED && Bac_ENABLED && Non_ENABLED
        call setReversed(array)
        call setCumSum(array)
#elif   Old_ENABLED && Bac_ENABLED && Rev_ENABLED
        do i = lenArray - 1, 1, -1
            array(i) = array(i + 1) + array(i)
        end do
#elif   New_ENABLED
        CHECK_ASSERTION(__LINE__, lenArray == size(cumsum, 1, IK), SK_"@setCumSum(): The condition `size(array, 1, IK) == size(cumsum, 1, IK)` must hold. size(array), size(cumsum) = "//getStr([lenArray, size(cumsum, 1, IK)]))
#if     For_ENABLED && Non_ENABLED
        cumsum(1) = array(1)
        do i = 2, lenArray
            cumsum(i) = cumsum(i - 1) + array(i)
        end do
#elif   For_ENABLED && Rev_ENABLED
        cumsum(lenArray) = array(1)
        do i = 2, lenArray
            cumsum(lenArray - i + 1) = cumsum(lenArray - i + 2) + array(i)
        end do
#elif   Bac_ENABLED && Non_ENABLED
        cumsum(1) = array(lenArray)
        do i = 2, lenArray
            cumsum(i) = cumsum(i - 1) + array(lenArray - i + 1)
        end do
#elif   Bac_ENABLED && Rev_ENABLED
        cumsum(lenArray) = array(lenArray)
        do i = lenArray - 1, 1, -1
            cumsum(i) = cumsum(i + 1) + array(i)
        end do
#else
#error  "Unrecognized interface."
#endif
#endif

#undef  TYPE_KIND

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif