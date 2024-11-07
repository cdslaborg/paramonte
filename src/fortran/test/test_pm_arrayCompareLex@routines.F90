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

!>  \brief This file contains the implementations of the tests of module [pm_arrayCompareLex](@ref pm_arrayCompareLex).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

submodule (test_pm_arrayCompareLex) routines

    !use pm_distUnif, only: setUnifRand
    use pm_val2str, only: getStr
    use pm_kind, only: LK, SK
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_llt_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_llt_D0_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_llt_D0_D0_SK5_1
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_llt_D0_D0_SK4_1
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_llt_D0_D0_SK3_1
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_llt_D0_D0_SK2_1
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_llt_D0_D0_SK1_1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_llt_D0_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_llt_D1_D1_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_llt_D1_D1_SK5_1
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_llt_D1_D1_SK4_1
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_llt_D1_D1_SK3_1
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_llt_D1_D1_SK2_1
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_llt_D1_D1_SK1_1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_llt_D1_D1_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_llt_D1_D1_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_llt_D1_D1_IK5_1
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_llt_D1_D1_IK4_1
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_llt_D1_D1_IK3_1
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_llt_D1_D1_IK2_1
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_llt_D1_D1_IK1_1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_llt_D1_D1_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_llt_D1_D1_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_llt_D1_D1_LK5_1
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_llt_D1_D1_LK4_1
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_llt_D1_D1_LK3_1
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_llt_D1_D1_LK2_1
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_llt_D1_D1_LK1_1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_llt_D1_D1_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_llt_D1_D1_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_llt_D1_D1_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_llt_D1_D1_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_llt_D1_D1_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_llt_D1_D1_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_llt_D1_D1_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_llt_D1_D1_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_llt_D1_D1_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_llt_D1_D1_RK5_1
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_llt_D1_D1_RK4_1
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_llt_D1_D1_RK3_1
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_llt_D1_D1_RK2_1
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_llt_D1_D1_RK1_1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_llt_D1_D1_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef test_llt_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lle_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lle_D0_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_lle_D0_D0_SK5_1
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_lle_D0_D0_SK4_1
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_lle_D0_D0_SK3_1
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_lle_D0_D0_SK2_1
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_lle_D0_D0_SK1_1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lle_D0_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lle_D1_D1_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_lle_D1_D1_SK5_1
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_lle_D1_D1_SK4_1
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_lle_D1_D1_SK3_1
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_lle_D1_D1_SK2_1
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_lle_D1_D1_SK1_1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lle_D1_D1_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lle_D1_D1_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_lle_D1_D1_IK5_1
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_lle_D1_D1_IK4_1
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_lle_D1_D1_IK3_1
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_lle_D1_D1_IK2_1
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_lle_D1_D1_IK1_1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lle_D1_D1_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lle_D1_D1_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_lle_D1_D1_LK5_1
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_lle_D1_D1_LK4_1
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_lle_D1_D1_LK3_1
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_lle_D1_D1_LK2_1
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_lle_D1_D1_LK1_1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lle_D1_D1_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lle_D1_D1_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_lle_D1_D1_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_lle_D1_D1_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_lle_D1_D1_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_lle_D1_D1_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_lle_D1_D1_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lle_D1_D1_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lle_D1_D1_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_lle_D1_D1_RK5_1
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_lle_D1_D1_RK4_1
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_lle_D1_D1_RK3_1
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_lle_D1_D1_RK2_1
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_lle_D1_D1_RK1_1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lle_D1_D1_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef test_lle_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lge_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lge_D0_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_lge_D0_D0_SK5_1
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_lge_D0_D0_SK4_1
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_lge_D0_D0_SK3_1
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_lge_D0_D0_SK2_1
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_lge_D0_D0_SK1_1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lge_D0_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lge_D1_D1_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_lge_D1_D1_SK5_1
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_lge_D1_D1_SK4_1
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_lge_D1_D1_SK3_1
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_lge_D1_D1_SK2_1
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_lge_D1_D1_SK1_1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lge_D1_D1_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lge_D1_D1_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_lge_D1_D1_IK5_1
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_lge_D1_D1_IK4_1
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_lge_D1_D1_IK3_1
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_lge_D1_D1_IK2_1
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_lge_D1_D1_IK1_1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lge_D1_D1_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lge_D1_D1_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_lge_D1_D1_LK5_1
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_lge_D1_D1_LK4_1
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_lge_D1_D1_LK3_1
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_lge_D1_D1_LK2_1
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_lge_D1_D1_LK1_1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lge_D1_D1_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lge_D1_D1_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_lge_D1_D1_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_lge_D1_D1_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_lge_D1_D1_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_lge_D1_D1_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_lge_D1_D1_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lge_D1_D1_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lge_D1_D1_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_lge_D1_D1_RK5_1
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_lge_D1_D1_RK4_1
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_lge_D1_D1_RK3_1
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_lge_D1_D1_RK2_1
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_lge_D1_D1_RK1_1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lge_D1_D1_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef test_lge_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lgt_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lgt_D0_D0_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_lgt_D0_D0_SK5_1
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_lgt_D0_D0_SK4_1
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_lgt_D0_D0_SK3_1
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_lgt_D0_D0_SK2_1
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_lgt_D0_D0_SK1_1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lgt_D0_D0_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lgt_D1_D1_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_lgt_D1_D1_SK5_1
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_lgt_D1_D1_SK4_1
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_lgt_D1_D1_SK3_1
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_lgt_D1_D1_SK2_1
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_lgt_D1_D1_SK1_1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lgt_D1_D1_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lgt_D1_D1_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_lgt_D1_D1_IK5_1
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_lgt_D1_D1_IK4_1
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_lgt_D1_D1_IK3_1
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_lgt_D1_D1_IK2_1
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_lgt_D1_D1_IK1_1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lgt_D1_D1_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lgt_D1_D1_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_lgt_D1_D1_LK5_1
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_lgt_D1_D1_LK4_1
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_lgt_D1_D1_LK3_1
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_lgt_D1_D1_LK2_1
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_lgt_D1_D1_LK1_1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lgt_D1_D1_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lgt_D1_D1_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_lgt_D1_D1_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_lgt_D1_D1_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_lgt_D1_D1_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_lgt_D1_D1_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_lgt_D1_D1_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lgt_D1_D1_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_lgt_D1_D1_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_lgt_D1_D1_RK5_1
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_lgt_D1_D1_RK4_1
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_lgt_D1_D1_RK3_1
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_lgt_D1_D1_RK2_1
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_lgt_D1_D1_RK1_1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayCompareLex@routines.inc.F90"
    end procedure
#endif

#undef test_lgt_D1_D1_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef test_lgt_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
