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

        !%%%%%%%%%%%%%%%%%%%%
#if     getCumPropExp_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

#if     Def_ENABLED
        type(sequence_type), parameter :: control = sequence_type()
#elif   !(Sel_ENABLED || Seq_ENABLED)
#error  "Unrecognized interface."
#endif
        if (0_IK < size(array, 1, IK)) then
            if (present(direction) .and. present(action)) then
                if (same_type_as(direction, forward)) then
                    if (same_type_as(action, nothing)) then
                        call setCumPropExp(cumPropExp, array, maxArray, control, forward, nothing)
                    elseif (same_type_as(action, reverse)) then
                        call setCumPropExp(cumPropExp, array, maxArray, control, forward, reverse)
                    else
                        error stop "@getCumPropExp(): Unrecognized action."
                    end if
                elseif (same_type_as(direction, backward)) then
                    if (same_type_as(action, nothing)) then
                        call setCumPropExp(cumPropExp, array, maxArray, control, backward, nothing)
                    elseif (same_type_as(action, reverse)) then
                        call setCumPropExp(cumPropExp, array, maxArray, control, backward, reverse)
                    else
                        error stop "@getCumPropExp(): Unrecognized action."
                    end if
                else
                    error stop "@getCumPropExp(): Unrecognized direction."
                end if
            elseif (present(direction)) then
                if (same_type_as(direction, forward)) then
                    call setCumPropExp(cumPropExp, array, maxArray, control, forward, nothing)
                elseif (same_type_as(direction, backward)) then
                    call setCumPropExp(cumPropExp, array, maxArray, control, backward, nothing)
                else
                    error stop "@getCumPropExp(): Unrecognized direction."
                end if
            elseif (present(action)) then
                if (same_type_as(action, nothing)) then
                    call setCumPropExp(cumPropExp, array, maxArray, control, forward, nothing)
                elseif (same_type_as(action, reverse)) then
                    call setCumPropExp(cumPropExp, array, maxArray, control, forward, reverse)
                else
                    error stop "@getCumPropExp(): Unrecognized action."
                end if
            else
                call setCumPropExp(cumPropExp, array, maxArray, control)
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCumPropExp_ENABLED && Old_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKC), parameter :: LOGTINY = log(tiny(0._RKC))
        integer(IK) :: lenArray
        lenArray = size(array, kind = IK)
        CHECK_ASSERTION(__LINE__, 0_IK < lenArray, SK_"@setCumPropExp(): The condition `0 < size(array)` must hold. size(array) = "//getStr(lenArray)) ! fpp
        CHECK_ASSERTION(__LINE__, maxval(array, 1) == maxArray, SK_"@setCumPropExp(): The condition `maxval(array, 1) == maxArray` must hold. maxval(array, 1), maxArray = "//getStr([maxval(array, 1), maxArray])) ! fpp
#if     For_ENABLED && Non_ENABLED
        block
            integer(IK) :: i
            real(RKC) :: cumPropExpInv
#if         Seq_ENABLED
            array(1) = exp(array(1) - maxArray)
#elif       Sel_ENABLED
            real(RKC) :: exponent
            exponent = array(1) - maxArray
            array(1) = 0._RKC
            if (LOGTINY < exponent) array(1) = exp(exponent)
#endif
            do i = 2, lenArray
#if             Seq_ENABLED
                array(i) = array(i - 1) + exp(array(i) - maxArray)
#elif           Sel_ENABLED
                exponent = array(i) - maxArray
                array(i) = array(i - 1)
                if (LOGTINY < exponent) array(i) = array(i) + exp(exponent)
#endif
            end do
            cumPropExpInv = 1._RKC / array(lenArray)
            do i = 1, lenArray - 1
                array(i) = array(i) * cumPropExpInv
            end do
            array(lenArray) = 1._RKC
        end block
#elif   For_ENABLED && Rev_ENABLED
        call setCumPropExp(array, maxArray, control)
        call setReversed(array)
#elif   Bac_ENABLED && Non_ENABLED
        call setReversed(array)
        call setCumPropExp(array, maxArray, control)
#elif   Bac_ENABLED && Rev_ENABLED
        block
            integer(IK) :: i
            real(RKC) :: cumPropExpInv
#if         Seq_ENABLED
            array(lenArray) = exp(array(lenArray) - maxArray)
#elif       Sel_ENABLED
            real(RKC) :: exponent
            exponent = array(lenArray) - maxArray
            array(lenArray) = 0._RKC
            if (LOGTINY < exponent) array(lenArray) = exp(exponent)
#endif
            do i = lenArray - 1, 1, -1
#if             Seq_ENABLED
                array(i) = array(i + 1) + exp(array(i) - maxArray)
#elif           Sel_ENABLED
                exponent = array(i) - maxArray
                array(i) = array(i + 1)
                if (LOGTINY < exponent) array(i) = array(i) + exp(exponent)
#endif
            end do
            cumPropExpInv = 1._RKC / array(1)
            array(1) = 1._RKC
            do i = 2, lenArray
                array(i) = array(i) * cumPropExpInv
            end do
        end block
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCumPropExp_ENABLED && New_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     Sel_ENABLED
        real(RKC) :: exponent
#elif   !Seq_ENABLED
#error  "Unrecognized interface."
#endif
        integer(IK) :: i, j
        integer(IK) :: lenArray
        real(RKC), parameter :: LOGTINY = log(tiny(0._RKC))
        real(RKC) :: cumPropExpInv
        lenArray = size(array, kind = IK)
        CHECK_ASSERTION(__LINE__, size(array, 1, IK) == size(cumPropExp, 1, IK), SK_"@setCumPropExp(): The condition `size(array, 1) == size(cumPropExp, 1)` must hold. size(array), size(cumPropExp) = "//getStr([size(array, 1, IK), size(cumPropExp,1, IK)]))
        CHECK_ASSERTION(__LINE__, maxval(array, 1) == maxArray, SK_"@setCumPropExp(): The condition `maxval(array, 1) == maxArray` must hold. maxval(array, 1), maxArray = "//getStr([maxval(array, 1), maxArray]))
        CHECK_ASSERTION(__LINE__, 0_IK < lenArray, SK_"@getCumPropExp()/setCumPropExp(): The condition `0 < size(array)` must hold. size(array) = "//getStr(lenArray))

#if     For_ENABLED && Non_ENABLED

#if     Seq_ENABLED
        cumPropExp(1) = exp(array(1) - maxArray)
#elif   Sel_ENABLED
        exponent = array(1) - maxArray
        cumPropExp(1) = 0._RKC
        if (LOGTINY < exponent) cumPropExp(1) = exp(exponent)
#endif
        do i = 2, lenArray
#if         Seq_ENABLED
            cumPropExp(i) = cumPropExp(i - 1) + exp(array(i) - maxArray)
#elif       Sel_ENABLED
            exponent = array(i) - maxArray
            if (LOGTINY < exponent) then
                cumPropExp(i) = cumPropExp(i - 1) + exp(exponent)
            else
                cumPropExp(i) = cumPropExp(i - 1)
            end if
#endif
        end do
        cumPropExpInv = 1._RKC / cumPropExp(lenArray)
        cumPropExp(lenArray) = 1._RKC
        do j = lenArray - 1, 1, -1
            cumPropExp(j) = cumPropExp(j) * cumPropExpInv
        end do

#elif   For_ENABLED && Rev_ENABLED

#if     Seq_ENABLED
        cumPropExp(lenArray) = exp(array(1) - maxArray)
#elif   Sel_ENABLED
        exponent = array(1) - maxArray
        cumPropExp(lenArray) = 0._RKC
        if (LOGTINY < exponent) cumPropExp(lenArray) = exp(exponent)
#endif
        do i = 2, lenArray
            j = lenArray - i + 1_IK
#if         Seq_ENABLED
            cumPropExp(j) = cumPropExp(j + 1) + exp(array(i) - maxArray)
#elif       Sel_ENABLED
            cumPropExp(j) = cumPropExp(j + 1)
            exponent = array(i) - maxArray
            if (LOGTINY < exponent) cumPropExp(j) = cumPropExp(j) + exp(exponent)
#endif
        end do
        cumPropExpInv = 1._RKC / cumPropExp(1)
        cumPropExp(1) = 1._RKC
        do j = 2, lenArray
            cumPropExp(j) = cumPropExp(j) * cumPropExpInv
        end do

#elif   Bac_ENABLED && Non_ENABLED

#if     Seq_ENABLED
        cumPropExp(1) = exp(array(lenArray) - maxArray)
#elif   Sel_ENABLED
        cumPropExp(1) = 0._RKC
        exponent = array(lenArray) - maxArray
        if (LOGTINY < exponent) cumPropExp(1) = exp(exponent)
#endif
        do i = 2, lenArray
#if         Seq_ENABLED
            cumPropExp(i) = cumPropExp(i - 1) + exp(array(lenArray - i + 1) - maxArray)
#elif       Sel_ENABLED
            cumPropExp(i) = cumPropExp(i - 1)
            exponent = array(lenArray - i + 1) - maxArray
            if (LOGTINY < exponent) cumPropExp(i) = cumPropExp(i) + exp(exponent)
#endif
        end do
        cumPropExpInv = 1._RKC / cumPropExp(lenArray)
        cumPropExp(lenArray) = 1._RKC
        do j = lenArray - 1, 1, -1
            cumPropExp(j) = cumPropExp(j) * cumPropExpInv
        end do

#elif   Bac_ENABLED && Rev_ENABLED

#if     Seq_ENABLED
        cumPropExp(lenArray) = exp(array(lenArray) - maxArray)
#elif   Sel_ENABLED
        cumPropExp(lenArray) = 0._RKC
        exponent = array(lenArray) - maxArray
        if (LOGTINY < exponent) cumPropExp(lenArray) = exp(exponent)
#endif
        do i = lenArray - 1, 1, -1
#if         Seq_ENABLED
            cumPropExp(i) = cumPropExp(i + 1) + exp(array(i) - maxArray)
#elif       Sel_ENABLED
            cumPropExp(i) = cumPropExp(i + 1)
            exponent = array(i) - maxArray
            if (LOGTINY < exponent) cumPropExp(i) = cumPropExp(i) + exp(exponent)
#endif
        end do
        cumPropExpInv = 1._RKC / cumPropExp(1)
        cumPropExp(1) = 1._RKC
        do j = 2, lenArray
            cumPropExp(j) = cumPropExp(j) * cumPropExpInv
        end do
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif