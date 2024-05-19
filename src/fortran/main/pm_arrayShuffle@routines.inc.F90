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
!>  This file contains the implementation details of the routines of the module [pm_arrayResize](@ref pm_arrayResize).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the indexing rule.
#if     D0_ENABLED && SK_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
#elif   D1_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE size
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%
#if     getShuffled_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, index(GET_SIZE(array, kind = IK))
        if (present(count)) then
            do concurrent(i = 1_IK : size(index, 1, IK))
                index(i) = i
            end do
            call setShuffled(index, count)
#if         D0_ENABLED && SK_ENABLED
            allocate(character(count,SKG) :: arrayShuffled)
#elif       D1_ENABLED && SK_ENABLED
            allocate(character(len(array,IK),SKG) :: arrayShuffled(count))
#elif       D1_ENABLED
            allocate(arrayShuffled(count))
#else
#error      "Unrecognized interface."
#endif
            do concurrent(i = 1_IK : count)
                arrayShuffled(GET_INDEX(i)) = array(GET_INDEX(index(i)))
            end do
        else
            arrayShuffled = array
            call setShuffled(arrayShuffled)
        end if

        !%%%%%%%%%%%%%%%%%%
#elif   setShuffled_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     D0_ENABLED && SK_ENABLED
        character(1,SKG) :: temp
#elif   D1_ENABLED && SK_ENABLED
        character(len(array,IK),SKG) :: temp
#elif   D1_ENABLED && IK_ENABLED
        integer(IKG) :: temp
#elif   D1_ENABLED && LK_ENABLED
        logical(LKG) :: temp
#elif   D1_ENABLED && CK_ENABLED
        complex(CKG) :: temp
#elif   D1_ENABLED && RK_ENABLED
        real(RKG) :: temp
#elif   D1_ENABLED && PSSK_ENABLED
        type(css_pdt(SKG)) :: temp
#elif   D1_ENABLED && BSSK_ENABLED
        type(css_type) :: temp
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: lenArray, index, randLoc
#if     RNGD_ENABLED
#define RNG
#elif   RNGF_ENABLED || RNGX_ENABLED
#define RNG rng,
#else
#error  "Unrecognized interface"
#endif
        !lenArray = GET_SIZE(array, kind = IK)
        !do index = lenArray, 2_IK, -1_IK
        !    call setUnifRand(rng, randLoc, 1_IK, index)
        !    temp = array(GET_INDEX(randLoc))
        !    array(GET_INDEX(randLoc)) = array(GET_INDEX(index))
        !    array(GET_INDEX(index)) = temp
        !end do
        lenArray = GET_SIZE(array, kind = IK)
        if (present(count)) then
            CHECK_ASSERTION(__LINE__, 0_IK <= count, SK_"@setShuffled(): The condition `0 <= count` must hold. count = "//getStr(count)) ! fpp
            CHECK_ASSERTION(__LINE__, count <= lenArray, SK_"@setShuffled(): The condition `count <= lenArray` must hold. count, lenArray = "//getStr([count, lenArray])) ! fpp
            randLoc = count
        else
            randLoc = lenArray - 1_IK
        end if
        do index = 1_IK, randLoc
            call setUnifRand(RNG randLoc, index, lenArray)
            temp = array(GET_INDEX(randLoc))
            array(GET_INDEX(randLoc)) = array(GET_INDEX(index))
            array(GET_INDEX(index)) = temp
        end do
#else
        !%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface"
        !%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_INDEX
#undef  GET_SIZE
#undef  RNG