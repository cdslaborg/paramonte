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
!>  This include file contains procedure implementation of [pm_mathNumSys](@ref pm_mathNumSys).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Sunday 3:33 AM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getDecimal_ENABLED || getNumeral_ENABLED
        character(*,SKC), parameter :: NUM_ALPHA = SKC_"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
#endif

        !%%%%%%%%%%%%%%%%%
#if     getDecimal_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IKC) :: i, numlen
        numlen = len(numeral, IKC)
        CHECK_ASSERTION(__LINE__, base > 0_IKC, SK_"@getDecimal(): The input `base` must be positive integer. base = "//getStr(base))
        CHECK_ASSERTION(__LINE__, numlen > 0_IKC, SK_"@getDecimal(): The input `numeral` must be non-empty. len(numeral) = "//getStr(numlen))
        !   \todo
        !   \pvhigh
        !   \warning
        !   The presence of `getStr()` in the followng statement `getStr(numeral), getStr(NUM_ALPHA(1:base))` is because of GNU gfortran bug as of version 11,
        !   which mistakenly converts the kind of a slices string to the default kind. It MUST be removed as soon as the gfortran bug is resolved.
        CHECK_ASSERTION(__LINE__, verify(getStr(numeral), getStr(NUM_ALPHA(1:base)), kind = IKC) == 0_IKC, \
        SK_"@getDecimal(): The input `numeral` must be only ASCII digits and the first `base` upper-case English alphabets. base, numeral, NUM_ALPHA(1:base) = " \
        //getStr(base)//SK_", "//getStr(numeral)//SK_", "//getStr(NUM_ALPHA(1:base)))
        decimal = 0_IKC
        do i = 1_IKC, numlen
            !   \todo
            !   \pvhigh
            !   \warning
            !   The following fence bypasses gfortran bug as of Version 11. It MUST be removed as soon as GFortran bug is resolved.
#if         __GFORTRAN__
            block
                character(1,SKC) :: numeral_char
                numeral_char = numeral(i:i)
            decimal = decimal + (index(NUM_ALPHA, numeral_char, kind = IKC) - 1_IKC) * base**(numlen - i)
            end block
#else
            decimal = decimal + (index(NUM_ALPHA, numeral(i:i), kind = IKC) - 1_IKC) * base**(numlen - i)
#endif
        end do

        !%%%%%%%%%%%%%%%%%
#elif   getNumeral_ENABLED
        !%%%%%%%%%%%%%%%%%

        !   \todo
        !   \warning
        !   `RK64` is good for handling integers as big as `huge(1_int1024)` that has 306 digits:
        !   7022238808055921514567598401519627865695222573993385049743362545223932648652381372371
        !   4248954065443758250044484324763030335464753443131493161268527593544579835065583369088
        !   0801860555545317367555154113605281582053784524026102900245630757473088050106395169337
        !   932361665227499793929447186391815763110662594625536
        !   As of 2022, the largest integer kind is made available by GNU Fortran compiler as `int128` whose huge has only 39 digits: 170141183460469231731687303715884105728
        integer, parameter :: RK64 = selected_real_kind(15)
        integer(IKC) :: remainderPlusOne, lamiced
        integer(IK) :: i, numlen
        CHECK_ASSERTION(__LINE__, range(0_IKC) < precision(0._RK64), \
        SK_"@getNumeral(): Hello, the distant future user. The condition `range(0_IKC) < precision(0._RK64)` must hold. This is an internal library limitation. Please report this to the developers if they are still alive.")
        CHECK_ASSERTION(__LINE__, base > 0_IKC, SK_"@getNumeral(): The input `base` must be a positive number. base = "//getStr(base))
        CHECK_ASSERTION(__LINE__, decimal > 0_IKC, SK_"@getNumeral(): The input `decimal` must be a positive number. decimal = "//getStr(decimal))
        numlen = int(log10(real(huge(1_IKC),RK64)) / log10(real(base,RK64)), IK) + 1_IK
        allocate(character(numlen,SKC) :: numeral)
        lamiced = decimal
        do i = numlen, 1_IK, -1_IK
            if(lamiced < base) then
                numeral(i:i) = NUM_ALPHA(lamiced + 1_IK : lamiced + 1_IK)
                exit
            end if
            remainderPlusOne = mod(lamiced, base) + 1_IKC
            numeral(i:i) = NUM_ALPHA(remainderPlusOne : remainderPlusOne)
            lamiced = lamiced / base
        end do
        CHECK_ASSERTION(__LINE__, i > 0_IK, SK_"@getNumeral(): Internal library error occured: i < 1. i = "//getStr(i))
        numeral = numeral(i:numlen)

        !%%%%%%%%%%%%%%%%%%%%
#elif   getCountDigit_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

#if     1
        integer(IKC) :: lav
        lav = val / 10_IKC
        count = 1_IK
        do
            if (lav == 0_IKC) exit
            count = count + 1_IK
            lav = lav / 10_IKC
        end do
#else
        !   \todo
        !   \warning
        !   `RK64` is good for handling integers as big as `huge(1_int1024)` that has 306 digits:
        !   7022238808055921514567598401519627865695222573993385049743362545223932648652381372371
        !   4248954065443758250044484324763030335464753443131493161268527593544579835065583369088
        !   0801860555545317367555154113605281582053784524026102900245630757473088050106395169337
        !   932361665227499793929447186391815763110662594625536
        !   As of 2022, the largest integer kind is made available by GNU Fortran compiler as `int128` whose huge has only 39 digits: 170141183460469231731687303715884105728
        !   \todo
        !   The performance of this algorithm can be improved by converting the integer values to appropriate reals with lower precision than RK64.
        integer, parameter :: RK64 = selected_real_kind(15)
        count = int(log10(real(abs(lav), RK64))) + 1_IK
#endif

        !%%%%%%%%%%%%%%%%%
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%