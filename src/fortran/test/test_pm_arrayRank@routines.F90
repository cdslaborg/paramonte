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
!>  This file contains procedure implementations of [pm_arrayRank](@ref pm_arrayRank).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_arrayRank) routines ! LCOV_EXCL_LINE

    use pm_val2str, only: getStr
    use pm_container, only: operator(>)
    use pm_logicalCompare, only: operator(>)
    use pm_complexCompareLex, only: operator(>)
    use pm_arraySort, only: setSorted
    use pm_distUnif, only: getUnifRand, setUnifRand
    use pm_arrayReverse, only: setReversed
    use pm_arrayRange, only: getRange
    use pm_kind, only: LK, RKR => RK
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Dense_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRankDense_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankDense_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankDense_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankDense_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankDense_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
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
    module procedure test_getRankDense_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankDense_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankDense_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankDense_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankDense_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getRankDense_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getRankDense_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getRankDense_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getRankDense_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getRankDense_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getRankDense_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getRankDense_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getRankDense_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getRankDense_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getRankDense_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getRankDense_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getRankDense_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getRankDense_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getRankDense_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getRankDense_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRankDense_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRankDense_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRankDense_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRankDense_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRankDense_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRankDense_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankDense_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankDense_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankDense_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankDense_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Dense_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Dense_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRankDense_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankDense_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankDense_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankDense_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankDense_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
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
    module procedure test_setRankDense_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankDense_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankDense_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankDense_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankDense_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setRankDense_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setRankDense_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setRankDense_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setRankDense_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setRankDense_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setRankDense_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setRankDense_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setRankDense_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setRankDense_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setRankDense_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setRankDense_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setRankDense_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setRankDense_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setRankDense_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setRankDense_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRankDense_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRankDense_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRankDense_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRankDense_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRankDense_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRankDense_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankDense_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankDense_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankDense_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankDense_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Dense_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fractional_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRankFractional_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankFractional_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankFractional_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankFractional_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankFractional_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
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
    module procedure test_getRankFractional_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankFractional_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankFractional_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankFractional_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankFractional_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getRankFractional_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getRankFractional_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getRankFractional_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getRankFractional_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getRankFractional_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getRankFractional_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getRankFractional_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getRankFractional_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getRankFractional_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getRankFractional_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getRankFractional_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getRankFractional_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getRankFractional_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getRankFractional_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getRankFractional_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRankFractional_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRankFractional_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRankFractional_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRankFractional_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRankFractional_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRankFractional_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankFractional_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankFractional_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankFractional_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankFractional_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fractional_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fractional_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRankFractional_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankFractional_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankFractional_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankFractional_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankFractional_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
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
    module procedure test_setRankFractional_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankFractional_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankFractional_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankFractional_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankFractional_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setRankFractional_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setRankFractional_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setRankFractional_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setRankFractional_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setRankFractional_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setRankFractional_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setRankFractional_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setRankFractional_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setRankFractional_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setRankFractional_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setRankFractional_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setRankFractional_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setRankFractional_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setRankFractional_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setRankFractional_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRankFractional_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRankFractional_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRankFractional_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRankFractional_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRankFractional_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRankFractional_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankFractional_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankFractional_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankFractional_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankFractional_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fractional_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Modified_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRankModified_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankModified_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankModified_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankModified_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankModified_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
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
    module procedure test_getRankModified_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankModified_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankModified_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankModified_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankModified_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getRankModified_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getRankModified_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getRankModified_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getRankModified_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getRankModified_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getRankModified_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getRankModified_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getRankModified_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getRankModified_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getRankModified_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getRankModified_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getRankModified_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getRankModified_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getRankModified_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getRankModified_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRankModified_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRankModified_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRankModified_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRankModified_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRankModified_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRankModified_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankModified_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankModified_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankModified_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankModified_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Modified_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Modified_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRankModified_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankModified_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankModified_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankModified_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankModified_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
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
    module procedure test_setRankModified_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankModified_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankModified_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankModified_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankModified_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setRankModified_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setRankModified_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setRankModified_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setRankModified_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setRankModified_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setRankModified_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setRankModified_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setRankModified_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setRankModified_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setRankModified_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setRankModified_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setRankModified_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setRankModified_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setRankModified_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setRankModified_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRankModified_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRankModified_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRankModified_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRankModified_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRankModified_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRankModified_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankModified_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankModified_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankModified_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankModified_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Modified_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ordinal_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRankOrdinal_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankOrdinal_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankOrdinal_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankOrdinal_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankOrdinal_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
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
    module procedure test_getRankOrdinal_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankOrdinal_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankOrdinal_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankOrdinal_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankOrdinal_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getRankOrdinal_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getRankOrdinal_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getRankOrdinal_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getRankOrdinal_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getRankOrdinal_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getRankOrdinal_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getRankOrdinal_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getRankOrdinal_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getRankOrdinal_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getRankOrdinal_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getRankOrdinal_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getRankOrdinal_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getRankOrdinal_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getRankOrdinal_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getRankOrdinal_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRankOrdinal_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRankOrdinal_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRankOrdinal_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRankOrdinal_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRankOrdinal_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRankOrdinal_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankOrdinal_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankOrdinal_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankOrdinal_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankOrdinal_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ordinal_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ordinal_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRankOrdinal_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankOrdinal_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankOrdinal_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankOrdinal_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankOrdinal_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
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
    module procedure test_setRankOrdinal_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankOrdinal_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankOrdinal_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankOrdinal_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankOrdinal_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setRankOrdinal_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setRankOrdinal_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setRankOrdinal_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setRankOrdinal_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setRankOrdinal_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setRankOrdinal_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setRankOrdinal_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setRankOrdinal_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setRankOrdinal_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setRankOrdinal_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setRankOrdinal_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setRankOrdinal_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setRankOrdinal_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setRankOrdinal_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setRankOrdinal_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRankOrdinal_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRankOrdinal_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRankOrdinal_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRankOrdinal_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRankOrdinal_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRankOrdinal_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankOrdinal_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankOrdinal_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankOrdinal_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankOrdinal_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ordinal_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Standard_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRankStandard_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankStandard_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankStandard_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankStandard_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankStandard_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
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
    module procedure test_getRankStandard_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankStandard_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankStandard_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankStandard_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankStandard_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getRankStandard_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getRankStandard_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getRankStandard_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getRankStandard_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getRankStandard_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getRankStandard_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getRankStandard_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getRankStandard_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getRankStandard_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getRankStandard_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getRankStandard_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getRankStandard_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getRankStandard_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getRankStandard_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getRankStandard_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRankStandard_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRankStandard_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRankStandard_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRankStandard_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRankStandard_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRankStandard_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRankStandard_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRankStandard_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRankStandard_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRankStandard_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Standard_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRank_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Standard_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRankStandard_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankStandard_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankStandard_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankStandard_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankStandard_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
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
    module procedure test_setRankStandard_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankStandard_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankStandard_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankStandard_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankStandard_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setRankStandard_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setRankStandard_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setRankStandard_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setRankStandard_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setRankStandard_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setRankStandard_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setRankStandard_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setRankStandard_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setRankStandard_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setRankStandard_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setRankStandard_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setRankStandard_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setRankStandard_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setRankStandard_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setRankStandard_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRankStandard_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRankStandard_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRankStandard_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRankStandard_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRankStandard_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRankStandard_D1_PSSK5
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRankStandard_D1_PSSK4
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRankStandard_D1_PSSK3
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRankStandard_D1_PSSK2
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRankStandard_D1_PSSK1
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
#include "test_pm_arrayRank@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Standard_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRank_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines