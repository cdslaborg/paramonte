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
!>  This file contains procedure implementations of [pm_arrayCompareLex](@ref pm_arrayCompareLex).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define logical comparison.
#if     LK_ENABLED
#define EQUIVALENT_TO .eqv.
#elif   SK_ENABLED || IK_ENABLED || CK_ENABLED || RK_ENABLED
#define EQUIVALENT_TO ==
#else
#error  "Unrecognized interface."
#endif
        ! Define array indexing rules.
#if     SK_ENABLED && D0_D0_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
#elif   D1_D1_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE size
#else
#error  "Unrecognized interface."
#endif
        ! Define comparison rules.
#if     (isllt_ENABLED || islle_ENABLED) && D1_D1_ENABLED && SK_ENABLED
#define COMPARABLE_TO .llt.
#elif   isllt_ENABLED || islle_ENABLED
#define COMPARABLE_TO <
#elif   (islgt_ENABLED || islge_ENABLED) && D1_D1_ENABLED && SK_ENABLED
#define COMPARABLE_TO .lgt.
#elif   islgt_ENABLED || islge_ENABLED
#define COMPARABLE_TO >
#endif
        integer(IK) :: i
        do i = 1_IK, min(GET_SIZE(array1, kind = IK), GET_SIZE(array2, kind = IK))
            if (array1(GET_INDEX(i)) EQUIVALENT_TO array2(GET_INDEX(i))) cycle
            compares = logical(array1(GET_INDEX(i)) COMPARABLE_TO array2(GET_INDEX(i)), kind = LK)
            return
        end do
        compares = logical(GET_SIZE(array1, kind = IK) & ! LCOV_EXCL_LINE
#if     isllt_ENABLED
        < & ! LCOV_EXCL_LINE
#elif   islle_ENABLED
        <= & ! LCOV_EXCL_LINE
#elif   islge_ENABLED
        >= & ! LCOV_EXCL_LINE
#elif   islgt_ENABLED
        > & ! LCOV_EXCL_LINE
#endif
        GET_SIZE(array2, kind = IK), kind = LK)
#undef  COMPARABLE_TO
#undef  EQUIVALENT_TO
#undef  GET_INDEX
#undef  GET_SIZE