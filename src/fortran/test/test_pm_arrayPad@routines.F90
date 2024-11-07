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
!>  This file contains procedure implementations of [test_pm_arrayPad](@ref test_pm_arrayPad).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_arrayPad) routines ! LCOV_EXCL_LINE

    use pm_kind, only: LK, SK
    use pm_val2str, only: getStr
    use pm_option, only: getOption
    use pm_distUnif, only: setUnifRand
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPadded_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getPadded_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getPadded_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getPadded_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getPadded_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getPadded_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
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
    module procedure test_getPadded_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getPadded_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getPadded_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getPadded_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getPadded_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getPadded_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getPadded_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getPadded_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getPadded_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getPadded_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getPadded_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getPadded_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getPadded_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getPadded_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getPadded_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getPadded_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getPadded_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getPadded_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getPadded_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getPadded_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPadded_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPadded_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPadded_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPadded_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPadded_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPadded_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPadded_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setPadded_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setPadded_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setPadded_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setPadded_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setPadded_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
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
    module procedure test_setPadded_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setPadded_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setPadded_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setPadded_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setPadded_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setPadded_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setPadded_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setPadded_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setPadded_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setPadded_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setPadded_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setPadded_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setPadded_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setPadded_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setPadded_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setPadded_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setPadded_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setPadded_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setPadded_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setPadded_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPadded_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPadded_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPadded_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPadded_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPadded_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPadded_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPaddedl_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getPaddedl_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getPaddedl_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getPaddedl_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getPaddedl_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getPaddedl_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
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
    module procedure test_getPaddedl_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getPaddedl_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getPaddedl_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getPaddedl_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getPaddedl_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getPaddedl_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getPaddedl_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getPaddedl_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getPaddedl_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getPaddedl_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getPaddedl_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getPaddedl_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getPaddedl_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getPaddedl_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getPaddedl_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getPaddedl_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getPaddedl_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getPaddedl_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getPaddedl_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getPaddedl_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPaddedl_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPaddedl_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPaddedl_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPaddedl_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPaddedl_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPaddedl_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPaddedl_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setPaddedl_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setPaddedl_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setPaddedl_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setPaddedl_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setPaddedl_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
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
    module procedure test_setPaddedl_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setPaddedl_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setPaddedl_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setPaddedl_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setPaddedl_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setPaddedl_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setPaddedl_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setPaddedl_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setPaddedl_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setPaddedl_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setPaddedl_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setPaddedl_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setPaddedl_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setPaddedl_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setPaddedl_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setPaddedl_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setPaddedl_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setPaddedl_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setPaddedl_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setPaddedl_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPaddedl_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPaddedl_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPaddedl_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPaddedl_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPaddedl_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPaddedl_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPaddedr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getPaddedr_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getPaddedr_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getPaddedr_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getPaddedr_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getPaddedr_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
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
    module procedure test_getPaddedr_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getPaddedr_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getPaddedr_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getPaddedr_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getPaddedr_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getPaddedr_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getPaddedr_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getPaddedr_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getPaddedr_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getPaddedr_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getPaddedr_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getPaddedr_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getPaddedr_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getPaddedr_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getPaddedr_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getPaddedr_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getPaddedr_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getPaddedr_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getPaddedr_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getPaddedr_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPaddedr_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPaddedr_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPaddedr_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPaddedr_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPaddedr_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPaddedr_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPaddedr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setPaddedr_D0_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setPaddedr_D0_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setPaddedr_D0_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setPaddedr_D0_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setPaddedr_D0_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
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
    module procedure test_setPaddedr_D1_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setPaddedr_D1_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setPaddedr_D1_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setPaddedr_D1_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setPaddedr_D1_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setPaddedr_D1_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setPaddedr_D1_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setPaddedr_D1_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setPaddedr_D1_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setPaddedr_D1_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setPaddedr_D1_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setPaddedr_D1_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setPaddedr_D1_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setPaddedr_D1_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setPaddedr_D1_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setPaddedr_D1_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setPaddedr_D1_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setPaddedr_D1_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setPaddedr_D1_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setPaddedr_D1_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPaddedr_D1_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPaddedr_D1_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPaddedr_D1_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPaddedr_D1_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPaddedr_D1_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_arrayPad@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPaddedr_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines