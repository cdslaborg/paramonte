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
!>  This file contains the implementation details of the routines under the generic interfaces in [pm_arrayStrip](@ref pm_arrayStrip).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%
#if     getStripped_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     SB_ENABLED && CusCom_ENABLED
        arrayStripped = array(getSIL(array, pattern, iseq) : getSIR(array, pattern, iseq))
#elif   SB_ENABLED && DefCom_ENABLED
        arrayStripped = array(getSIL(array, pattern) : getSIR(array, pattern))
#elif   SL_ENABLED && CusCom_ENABLED
        arrayStripped = array(getSIL(array, pattern, iseq) : )
#elif   SL_ENABLED && DefCom_ENABLED
        arrayStripped = array(getSIL(array, pattern) : )
#elif   SR_ENABLED && CusCom_ENABLED
        arrayStripped = array(1 : getSIR(array, pattern, iseq))
#elif   SR_ENABLED && DefCom_ENABLED
        arrayStripped = array(1 : getSIR(array, pattern))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getSIL_ENABLED || getSIR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: lenArray
#if     SK_ENABLED && D0_D0_ENABLED
#define GET_SIZE len
#else
#define GET_SIZE size
#endif
        ! Define the incrementing rules.
#if     getSIL_ENABLED
#define INCREMENT(a,b) a + b
#define START_INDEX 1_IK
#define STOP_INDEX lenArray - lenPatternMinusOne
#define OUT_OF_RANGE(a,b) a > b
        !integer(IK) , parameter :: SIGN = +1_IK
#elif   getSIR_ENABLED
#define OUT_OF_RANGE(a,b) a < b
#define INCREMENT(a,b) a - b
#define START_INDEX lenArray
#define STOP_INDEX lenPattern
        !integer(IK) , parameter :: SIGN = -1_IK
#else
#error  "Unrecognized interface."
#endif
        ! Define the indexing rules.
#if     D1_D0_ENABLED
#define GET_INDEX(i) i
        integer(IK) , parameter :: lenPattern = 1_IK, lenPatternMinusOne = 0_IK
#elif   D0_D0_ENABLED || D1_D1_ENABLED
#if     getSIL_ENABLED
#define GET_INDEX(i) i : i + lenPatternMinusOne
#elif   getSIR_ENABLED
#define GET_INDEX(i) i - lenPatternMinusOne : i
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: lenPattern, lenPatternMinusOne
        lenPattern = GET_SIZE(pattern, kind = IK) ! fpp
        lenPatternMinusOne = lenPattern - 1_IK
#else
#error  "Unrecognized interface."
#endif
        ! The order of conditions should not be changed here.
#if     CusCom_ENABLED && (D0_D0_ENABLED || D1_D0_ENABLED)
#define IS_EQ(a,b,lenb) iseq(a,b)
#elif   CusCom_ENABLED && D1_D1_ENABLED
#define IS_EQ(a,b,lenb) iseq(a,b,lenb)
#elif   DefCom_ENABLED && D1_D0_ENABLED && LK_ENABLED
#define IS_EQ(a,b,lenb) a .eqv. b
#elif   DefCom_ENABLED && D1_D1_ENABLED && LK_ENABLED
#define IS_EQ(a,b,lenb) all(a .eqv. b)
#elif   DefCom_ENABLED && (D1_D0_ENABLED || D0_D0_ENABLED)
#define IS_EQ(a,b,lenb) a == b
#elif   DefCom_ENABLED && D1_D1_ENABLED
#define IS_EQ(a,b,lenb) all(a == b)
#else
#error  "Unrecognized interface."
#endif
        lenArray = GET_SIZE(array, kind = IK) ! fpp
        if (lenArray < lenPattern .or. lenArray == 0_IK .or. lenPattern == 0_IK) then
#if         getSIL_ENABLED
            index = 1_IK
#elif       getSIR_ENABLED
            index = lenArray
#endif
            return
        end if
        index = START_INDEX
        do
            if (.not. IS_EQ(array(GET_INDEX(index)), pattern, lenPattern)) return
            index = INCREMENT(index, lenPattern) ! fpp
            if (OUT_OF_RANGE(index, STOP_INDEX)) return
        end do
#undef  OUT_OF_RANGE
#undef  START_INDEX
#undef  STOP_INDEX
#undef  INCREMENT
#undef  GET_INDEX
#undef  GET_SIZE
#undef  IS_EQ

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif