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
!>  This include file contains procedure implementations of [pm_arrayCompact](@ref pm_arrayCompact).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Saturday 1:48 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define logical comparison.
#if     LK_ENABLED
#define IS_NEQ(a,b) a .neqv. b
#elif   SK_ENABLED || IK_ENABLED || CK_ENABLED || RK_ENABLED
#define IS_NEQ(a,b) a /= b
#else
#error  "Unrecognized interface."
#endif
        ! Define array indexing rule.
#if     SK_ENABLED && D0_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
#elif   D1_ENABLED || D2_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE size
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: ip
#if     D2_ENABLED
        integer(IK) :: nd, np
#endif
#if     getCompact_ENABLED
#define EVALUATE(THIS)
#define ARRAY compact
        integer(IK) :: csize
#if     SK_ENABLED && D0_ENABLED
        allocate(character(len(array,IK),SKG) :: compact)
#else
        allocate(compact, mold = array)
#endif
#elif   setCompact_ENABLED
#define EVALUATE(THIS) THIS
#else
#error  "Unrecognized interface."
#endif
        if (GET_SIZE(array, kind = IK) == 0_IK) then
#if         setCompact_ENABLED
            csize = 0_IK
#elif       !getCompact_ENABLED
#error      "Unrecognized interface."
#endif
            return
        end if
        !%%%%%%%%%%%%%%%%%%%%%%%
#if     D0_ENABLED || D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%
        csize = 1_IK
#if     getCompact_ENABLED
        ARRAY(GET_INDEX(1)) = array(GET_INDEX(1)) ! fpp
#elif   setCompact_ENABLED
        CHECK_ASSERTION(__LINE__, size(weight, kind = IK) == GET_SIZE(array, kind = IK), \
        SK_"The size of `weight` must equal the size of `array`. size(array), size(weight) = "\
        //getStr([GET_SIZE(array, kind = IK), size(weight, kind = IK)])) ! fpp
        weight(csize) = 1_IK
#else
#error  "Unrecognized interface."
#endif
        do ip = 2_IK, GET_SIZE(array, kind = IK) ! fpp
            if (IS_NEQ(array(GET_INDEX(ip-1_IK)), array(GET_INDEX(ip)))) then ! fpp
                csize = csize + 1_IK
                EVALUATE(weight(csize) = 1_IK) ! fpp
                EVALUATE(if (csize /= ip)) ARRAY(GET_INDEX(csize)) = array(GET_INDEX(ip)) ! fpp
            else
                EVALUATE(weight(csize) = weight(csize) + 1_IK)
            end if
        end do
#if     getCompact_ENABLED
        ARRAY = ARRAY(1:csize) ! fpp
#endif
        !%%%%%%%%%
#elif   D2_ENABLED
        !%%%%%%%%%
        CHECK_ASSERTION(__LINE__, dim == 1_IK .or. dim == 2_IK, \
        SK_"The input `dim` must be either 1 or 2. dim = "//getStr(dim)) ! fpp
#if     setCompact_ENABLED
        CHECK_ASSERTION(__LINE__, size(weight, kind = IK) == size(array, dim, IK), \
        SK_"The size of `weight` must equal the size of `array` along dimension `dim`. dim, size(array, dim), size(weight) = "\
        //getStr([dim, size(array, dim, IK), size(weight, 1, IK)])) ! fpp
#endif
        np = size(array, dim, IK)
        if (dim == 2_IK) then
            nd = size(array, 1_IK, IK)
            csize = 1_IK
#if         getCompact_ENABLED
            ARRAY(1:nd,GET_INDEX(1)) = array(1:nd,GET_INDEX(1)) ! fpp
#endif
            EVALUATE(weight(csize) = 1_IK)
            do ip = 2_IK, np
                if (any(IS_NEQ(array(1:nd,ip-1), array(1:nd,ip)))) then ! fpp
                    csize = csize + 1_IK
                    EVALUATE(weight(csize) = 1_IK) ! fpp
                    EVALUATE(if (csize /= ip)) ARRAY(1:nd,csize) = array(1:nd,ip) ! fpp
                else
                    EVALUATE(weight(csize) = weight(csize) + 1_IK) ! fpp
                end if
            end do
#if         getCompact_ENABLED
            ARRAY = ARRAY(1:nd,1:csize) ! fpp
#endif
        elseif (dim == 1_IK) then
            nd = size(array, dim = 2_IK, kind = IK)
            csize = 1_IK
#if         getCompact_ENABLED
            ARRAY(GET_INDEX(1),1:nd) = array(GET_INDEX(1),1:nd) ! fpp
#endif
            EVALUATE(weight(csize) = 1_IK) ! fpp
            do ip = 2_IK, np
                if (any(IS_NEQ(array(ip-1,1:nd), array(ip,1:nd)))) then ! fpp
                    csize = csize + 1_IK
                    EVALUATE(weight(csize) = 1_IK) ! fpp
                    EVALUATE(if (csize /= ip)) ARRAY(csize,1:nd) = array(ip,1:nd) ! fpp
                else
                    EVALUATE(weight(csize) = weight(csize) + 1_IK) ! fpp
                end if
            end do
#if         getCompact_ENABLED
            ARRAY = ARRAY(1:csize,1:nd) ! fpp
#endif
        end if
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  CHECK_SUM_WEIGHT
#undef  GET_INDEX
#undef  GET_SIZE
#undef  EVALUATE
#undef  IS_NEQ
#undef  ARRAY