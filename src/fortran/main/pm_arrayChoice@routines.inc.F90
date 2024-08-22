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
!>  This file contains the implementation details of the routines under the generic interface [setChoice](@ref pm_arrayResize::setChoice).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%
#if     getChoice_ENABLED
        !%%%%%%%%%%%%%%%%

#if     D0_D0_ENABLED || D1_D0_ENABLED
        call setChoice(rngf, choice, array)
#elif   D0_S1_ENABLED || D1_D1_ENABLED
        call setChoice(rngf, choice, array, unique)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setChoice_ENABLED && Def_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: index, lenArray
#if     D0_D0_ENABLED || D1_D1_ENABLED
        integer(IK) :: ichoice, lenChoice
#endif
        ! Check string length compatibility.
#if     SK_ENABLED && D0_D0_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
        !check_assertion(__LINE__, 0_IK < len(choice, IK), SK_"@setChoice(): The condition `1 == len(choice)` must hold. len(choice) = "//getStr(len(choice, IK))) ! fpp
#elif   SK_ENABLED || IK_ENABLED || LK_ENABLED || CK_ENABLED || RK_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE size
#if     SK_ENABLED
        CHECK_ASSERTION(__LINE__, len(array, IK) == len(choice, IK), SK_"@setChoice(): The condition `len(array) == len(choice)` must hold. len(array), len(choice) = "//getStr([len(array, IK), len(choice, IK)])) ! fpp
#endif
#else
#error  "Unrecognized interface."
#endif
        lenArray = GET_SIZE(array, kind = IK)
        CHECK_ASSERTION(__LINE__, 0_IK < lenArray, SK_"@setChoice(): The length of the input `array` must be non-zero. lenArray = "//getStr(lenArray)) ! fpp
#if     D1_D0_ENABLED
        if (1_IK < lenArray) then
            call setUnifRand(rng, index, 1_IK, lenArray)
            choice = array(GET_INDEX(index))
        else
            choice = array(GET_INDEX(1))
        end if
#elif   D0_D0_ENABLED || D1_D1_ENABLED
        lenChoice = GET_SIZE(choice, kind = IK)
        !check_assertion(__LINE__, 0_IK < lenChoice, SK_"@setChoice(): The length of the input `choice` must be non-zero. lenChoice = "//getStr(lenChoice)) ! fpp
        if (present(unique)) then
            if (unique) then
                CHECK_ASSERTION(__LINE__, lenChoice <= lenArray, SK_"@setChoice(): The size of the input `choice` must smaller than or equal to the size of `array`. lenChoice, lenArray = "//getStr([lenChoice, lenArray])) ! fpp
                block
                    integer(IK) :: shuffle(lenArray)
                    call setRange(shuffle, 1_IK)
                    call setShuffled(rng, shuffle, lenChoice)
                    do concurrent(index = 1 : lenChoice)
                        choice(GET_INDEX(index)) = array(GET_INDEX(shuffle(index)))
                    end do
                    return
                end block
                return
            end if
        end if
        if (1_IK < lenArray) then
            ! \todo
            ! This must be improved. No need for uniform CDF in the default case. Use simple shuffling.
            block
                real(RKD) :: unifrnd, cdf(lenArray), lenArray_RKD
                lenArray_RKD = real(lenArray, RKD)
                call setLinSpace(cdf, 0._RKD, (lenArray_RKD - 1._RKD) / lenArray_RKD)
                do ichoice = 1, lenChoice
                    call setUnifRand(rng, unifrnd)
                    index = getBin(cdf, unifrnd)
                    choice(GET_INDEX(ichoice)) = array(GET_INDEX(index))
                end do
            end block
        else
            do concurrent(index = 1 : lenChoice)
                choice(GET_INDEX(index)) = array(GET_INDEX(1))
            end do
        end if
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface"
        !%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  GET_INDEX
#undef  GET_SIZE
