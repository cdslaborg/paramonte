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
!>  This file contains procedure implementations of [pm_mathMinMax](@ref pm_mathMinMax).
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

submodule (pm_mathMinMax) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMinMax_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Indi_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMinMaxIndi_SK5
        use pm_kind, only: SKG => SK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMinMaxIndi_SK4
        use pm_kind, only: SKG => SK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMinMaxIndi_SK3
        use pm_kind, only: SKG => SK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMinMaxIndi_SK2
        use pm_kind, only: SKG => SK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMinMaxIndi_SK1
        use pm_kind, only: SKG => SK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMinMaxIndi_IK5
        use pm_kind, only: IKG => IK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMinMaxIndi_IK4
        use pm_kind, only: IKG => IK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMinMaxIndi_IK3
        use pm_kind, only: IKG => IK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMinMaxIndi_IK2
        use pm_kind, only: IKG => IK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMinMaxIndi_IK1
        use pm_kind, only: IKG => IK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMinMaxIndi_LK5
        use pm_kind, only: LKG => LK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMinMaxIndi_LK4
        use pm_kind, only: LKG => LK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMinMaxIndi_LK3
        use pm_kind, only: LKG => LK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMinMaxIndi_LK2
        use pm_kind, only: LKG => LK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMinMaxIndi_LK1
        use pm_kind, only: LKG => LK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMinMaxIndi_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMinMaxIndi_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMinMaxIndi_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMinMaxIndi_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMinMaxIndi_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMinMaxIndi_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMinMaxIndi_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMinMaxIndi_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMinMaxIndi_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMinMaxIndi_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Indi_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Pair_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getMinMaxPair_SK5
        use pm_kind, only: SKG => SK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getMinMaxPair_SK4
        use pm_kind, only: SKG => SK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getMinMaxPair_SK3
        use pm_kind, only: SKG => SK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getMinMaxPair_SK2
        use pm_kind, only: SKG => SK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getMinMaxPair_SK1
        use pm_kind, only: SKG => SK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getMinMaxPair_IK5
        use pm_kind, only: IKG => IK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getMinMaxPair_IK4
        use pm_kind, only: IKG => IK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getMinMaxPair_IK3
        use pm_kind, only: IKG => IK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getMinMaxPair_IK2
        use pm_kind, only: IKG => IK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getMinMaxPair_IK1
        use pm_kind, only: IKG => IK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getMinMaxPair_LK5
        use pm_kind, only: LKG => LK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getMinMaxPair_LK4
        use pm_kind, only: LKG => LK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getMinMaxPair_LK3
        use pm_kind, only: LKG => LK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getMinMaxPair_LK2
        use pm_kind, only: LKG => LK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getMinMaxPair_LK1
        use pm_kind, only: LKG => LK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMinMaxPair_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMinMaxPair_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMinMaxPair_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMinMaxPair_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMinMaxPair_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMinMaxPair_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMinMaxPair_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMinMaxPair_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMinMaxPair_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMinMaxPair_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Pair_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMinMax_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMinMax_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Indi_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMinMaxIndi_SK5
        use pm_kind, only: SKG => SK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMinMaxIndi_SK4
        use pm_kind, only: SKG => SK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMinMaxIndi_SK3
        use pm_kind, only: SKG => SK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMinMaxIndi_SK2
        use pm_kind, only: SKG => SK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMinMaxIndi_SK1
        use pm_kind, only: SKG => SK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMinMaxIndi_IK5
        use pm_kind, only: IKG => IK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMinMaxIndi_IK4
        use pm_kind, only: IKG => IK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMinMaxIndi_IK3
        use pm_kind, only: IKG => IK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMinMaxIndi_IK2
        use pm_kind, only: IKG => IK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMinMaxIndi_IK1
        use pm_kind, only: IKG => IK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMinMaxIndi_LK5
        use pm_kind, only: LKG => LK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMinMaxIndi_LK4
        use pm_kind, only: LKG => LK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMinMaxIndi_LK3
        use pm_kind, only: LKG => LK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMinMaxIndi_LK2
        use pm_kind, only: LKG => LK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMinMaxIndi_LK1
        use pm_kind, only: LKG => LK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMinMaxIndi_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMinMaxIndi_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMinMaxIndi_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMinMaxIndi_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMinMaxIndi_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMinMaxIndi_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMinMaxIndi_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMinMaxIndi_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMinMaxIndi_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMinMaxIndi_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Indi_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Pair_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setMinMaxPair_SK5
        use pm_kind, only: SKG => SK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setMinMaxPair_SK4
        use pm_kind, only: SKG => SK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setMinMaxPair_SK3
        use pm_kind, only: SKG => SK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setMinMaxPair_SK2
        use pm_kind, only: SKG => SK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setMinMaxPair_SK1
        use pm_kind, only: SKG => SK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setMinMaxPair_IK5
        use pm_kind, only: IKG => IK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setMinMaxPair_IK4
        use pm_kind, only: IKG => IK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setMinMaxPair_IK3
        use pm_kind, only: IKG => IK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setMinMaxPair_IK2
        use pm_kind, only: IKG => IK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setMinMaxPair_IK1
        use pm_kind, only: IKG => IK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure setMinMaxPair_LK5
        use pm_kind, only: LKG => LK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure setMinMaxPair_LK4
        use pm_kind, only: LKG => LK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure setMinMaxPair_LK3
        use pm_kind, only: LKG => LK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure setMinMaxPair_LK2
        use pm_kind, only: LKG => LK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure setMinMaxPair_LK1
        use pm_kind, only: LKG => LK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMinMaxPair_CK5
        use pm_kind, only: CKG => CK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMinMaxPair_CK4
        use pm_kind, only: CKG => CK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMinMaxPair_CK3
        use pm_kind, only: CKG => CK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMinMaxPair_CK2
        use pm_kind, only: CKG => CK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMinMaxPair_CK1
        use pm_kind, only: CKG => CK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMinMaxPair_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMinMaxPair_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMinMaxPair_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMinMaxPair_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMinMaxPair_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathMinMax@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Pair_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMinMax_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines