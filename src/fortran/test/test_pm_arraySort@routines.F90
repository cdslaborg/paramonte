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
!>  This file contains procedure implementations of [pm_arraySort](@ref pm_arraySort).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_arraySort) routines ! LCOV_EXCL_LINE

    use pm_kind, only: LK
    use pm_val2str, only: getStr
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_arrayRemap, only: getRemapped
    use pm_arrayResize, only: setResized
    use pm_arrayReverse, only: getReversed
    use pm_logicalCompare, only: operator(>)
    use pm_complexCompareLex, only: operator(>)
    use pm_container, only: operator(>)
    !use pm_container, only: css_pdt
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isAscending_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_isAscending_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_isAscending_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_isAscending_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_isAscending_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_isAscending_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_isAscending_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_isAscending_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_isAscending_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_isAscending_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_isAscending_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_isAscending_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_isAscending_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_isAscending_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_isAscending_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_isAscending_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_isAscending_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_isAscending_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_isAscending_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_isAscending_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_isAscending_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isAscending_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isAscending_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isAscending_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isAscending_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isAscending_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_isAscending_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_isAscending_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_isAscending_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_isAscending_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_isAscending_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_isAscending_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_isAscending_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_isAscending_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_isAscending_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_isAscending_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isAscending_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isDescending_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_isDescending_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_isDescending_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_isDescending_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_isDescending_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_isDescending_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_isDescending_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_isDescending_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_isDescending_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_isDescending_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_isDescending_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_isDescending_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_isDescending_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_isDescending_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_isDescending_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_isDescending_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_isDescending_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_isDescending_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_isDescending_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_isDescending_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_isDescending_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isDescending_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isDescending_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isDescending_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isDescending_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isDescending_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_isDescending_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_isDescending_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_isDescending_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_isDescending_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_isDescending_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_isDescending_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_isDescending_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_isDescending_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_isDescending_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_isDescending_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isDescending_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_isSorted_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_isSorted_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_isSorted_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_isSorted_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_isSorted_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_isSorted_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_isSorted_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_isSorted_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_isSorted_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_isSorted_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_isSorted_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_isSorted_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_isSorted_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_isSorted_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_isSorted_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_isSorted_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_isSorted_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_isSorted_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_isSorted_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_isSorted_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_isSorted_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_isSorted_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_isSorted_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_isSorted_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_isSorted_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_isSorted_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_isSorted_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_isSorted_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_isSorted_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_isSorted_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_isSorted_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_isSorted_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_isSorted_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_isSorted_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_isSorted_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt, operator(<=)
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ind_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedIndDef_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedIndDef_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedIndDef_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedIndDef_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedIndDef_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedIndDef_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedIndDef_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedIndDef_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedIndDef_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedIndDef_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedIndDef_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedIndDef_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedIndDef_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedIndDef_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedIndDef_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedIndDef_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedIndDef_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedIndDef_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedIndDef_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedIndDef_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedIndDef_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedIndDef_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedIndDef_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedIndDef_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedIndDef_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedIndDef_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedIndDef_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedIndDef_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedIndDef_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedIndDef_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedIndDef_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedIndDef_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedIndDef_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedIndDef_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedIndDef_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ind_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Qsorti_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrQsorti_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrQsorti_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrQsorti_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrQsorti_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrQsorti_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrQsorti_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrQsorti_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrQsorti_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrQsorti_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrQsorti_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrQsorti_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrQsorti_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrQsorti_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrQsorti_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrQsorti_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrQsorti_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrQsorti_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrQsorti_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrQsorti_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrQsorti_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrQsorti_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrQsorti_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrQsorti_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrQsorti_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrQsorti_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrQsorti_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrQsorti_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrQsorti_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrQsorti_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrQsorti_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrQsorti_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrQsorti_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrQsorti_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrQsorti_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrQsorti_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Qsorti_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Qsortr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrQsortr_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrQsortr_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrQsortr_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrQsortr_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrQsortr_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrQsortr_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrQsortr_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrQsortr_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrQsortr_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrQsortr_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrQsortr_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrQsortr_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrQsortr_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrQsortr_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrQsortr_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrQsortr_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrQsortr_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrQsortr_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrQsortr_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrQsortr_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrQsortr_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrQsortr_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrQsortr_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrQsortr_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrQsortr_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrQsortr_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrQsortr_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrQsortr_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrQsortr_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrQsortr_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrQsortr_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrQsortr_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrQsortr_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrQsortr_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrQsortr_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Qsortr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Qsortrdp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrQsortrdp_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrQsortrdp_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrQsortrdp_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrQsortrdp_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrQsortrdp_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrQsortrdp_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Qsortrdp_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Bubble_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrBubble_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrBubble_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrBubble_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrBubble_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrBubble_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrBubble_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrBubble_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrBubble_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrBubble_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrBubble_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrBubble_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrBubble_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrBubble_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrBubble_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrBubble_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrBubble_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrBubble_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrBubble_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrBubble_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrBubble_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrBubble_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrBubble_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrBubble_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrBubble_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrBubble_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrBubble_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrBubble_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrBubble_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrBubble_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrBubble_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrBubble_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrBubble_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrBubble_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrBubble_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrBubble_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Bubble_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Heapi_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrHeapi_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrHeapi_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrHeapi_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrHeapi_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrHeapi_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrHeapi_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrHeapi_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrHeapi_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrHeapi_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrHeapi_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrHeapi_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrHeapi_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrHeapi_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrHeapi_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrHeapi_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrHeapi_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrHeapi_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrHeapi_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrHeapi_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrHeapi_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrHeapi_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrHeapi_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrHeapi_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrHeapi_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrHeapi_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrHeapi_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrHeapi_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrHeapi_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrHeapi_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrHeapi_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrHeapi_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrHeapi_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrHeapi_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrHeapi_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrHeapi_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Heapi_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Heapr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrHeapr_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrHeapr_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrHeapr_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrHeapr_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrHeapr_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrHeapr_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrHeapr_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrHeapr_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrHeapr_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrHeapr_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrHeapr_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrHeapr_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrHeapr_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrHeapr_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrHeapr_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrHeapr_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrHeapr_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrHeapr_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrHeapr_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrHeapr_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrHeapr_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrHeapr_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrHeapr_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrHeapr_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrHeapr_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrHeapr_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrHeapr_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrHeapr_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrHeapr_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrHeapr_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrHeapr_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrHeapr_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrHeapr_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrHeapr_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrHeapr_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Heapr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Insertionl_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrInsertionl_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrInsertionl_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrInsertionl_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrInsertionl_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrInsertionl_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrInsertionl_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrInsertionl_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrInsertionl_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrInsertionl_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrInsertionl_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrInsertionl_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrInsertionl_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrInsertionl_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrInsertionl_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrInsertionl_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrInsertionl_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrInsertionl_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrInsertionl_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrInsertionl_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrInsertionl_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrInsertionl_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrInsertionl_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrInsertionl_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrInsertionl_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrInsertionl_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrInsertionl_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrInsertionl_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrInsertionl_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrInsertionl_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrInsertionl_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrInsertionl_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrInsertionl_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrInsertionl_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrInsertionl_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrInsertionl_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Insertionl_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Insertionb_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrInsertionb_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrInsertionb_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrInsertionb_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrInsertionb_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrInsertionb_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrInsertionb_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrInsertionb_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrInsertionb_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrInsertionb_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrInsertionb_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrInsertionb_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrInsertionb_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrInsertionb_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrInsertionb_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrInsertionb_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrInsertionb_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrInsertionb_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrInsertionb_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrInsertionb_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrInsertionb_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrInsertionb_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrInsertionb_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrInsertionb_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrInsertionb_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrInsertionb_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrInsertionb_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrInsertionb_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrInsertionb_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrInsertionb_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrInsertionb_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrInsertionb_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrInsertionb_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrInsertionb_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrInsertionb_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrInsertionb_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Insertionb_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Merger_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrMerger_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrMerger_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrMerger_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrMerger_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrMerger_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrMerger_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrMerger_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrMerger_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrMerger_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrMerger_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrMerger_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrMerger_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrMerger_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrMerger_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrMerger_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrMerger_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrMerger_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrMerger_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrMerger_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrMerger_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrMerger_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrMerger_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrMerger_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrMerger_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrMerger_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrMerger_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrMerger_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrMerger_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrMerger_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrMerger_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrMerger_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrMerger_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrMerger_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrMerger_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrMerger_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Merger_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Selection_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrSelection_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrSelection_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrSelection_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrSelection_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrSelection_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrSelection_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrSelection_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrSelection_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrSelection_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrSelection_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrSelection_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrSelection_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrSelection_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrSelection_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrSelection_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrSelection_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrSelection_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrSelection_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrSelection_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrSelection_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrSelection_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrSelection_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrSelection_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrSelection_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrSelection_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrSelection_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrSelection_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrSelection_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrSelection_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrSelection_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrSelection_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrSelection_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrSelection_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrSelection_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrSelection_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Selection_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setSorted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Arr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Shell_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrShell_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrShell_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrShell_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrShell_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrShell_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrShell_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrShell_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrShell_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrShell_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrShell_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setSortedArrShell_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setSortedArrShell_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setSortedArrShell_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setSortedArrShell_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setSortedArrShell_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setSortedArrShell_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setSortedArrShell_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setSortedArrShell_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setSortedArrShell_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setSortedArrShell_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setSortedArrShell_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setSortedArrShell_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setSortedArrShell_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setSortedArrShell_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setSortedArrShell_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setSortedArrShell_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setSortedArrShell_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setSortedArrShell_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setSortedArrShell_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setSortedArrShell_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED && 0
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setSortedArrShell_D1_PSSK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setSortedArrShell_D1_PSSK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setSortedArrShell_D1_PSSK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setSortedArrShell_D1_PSSK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setSortedArrShell_D1_PSSK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arraySort@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Shell_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Arr_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setSorted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines