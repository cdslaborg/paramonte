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
!>  This file contains the implementations of the tests of module [pm_distUnif](@ref pm_distUnif).
!>
!>  \final
!>
!>  \author 
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distUnif) routines

    use pm_sampleMean, only: getMean
    use pm_sampleVar, only: getVar
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getUnifCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getUnifCDF_IK5
        use pm_kind, only: RKG => RK, IKG => IK5
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getUnifCDF_IK4
        use pm_kind, only: RKG => RK, IKG => IK4
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getUnifCDF_IK3
        use pm_kind, only: RKG => RK, IKG => IK3
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getUnifCDF_IK2
        use pm_kind, only: RKG => RK, IKG => IK2
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getUnifCDF_IK1
        use pm_kind, only: RKG => RK, IKG => IK1
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getUnifCDF_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getUnifCDF_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getUnifCDF_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getUnifCDF_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getUnifCDF_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getUnifCDF_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getUnifCDF_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getUnifCDF_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getUnifCDF_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getUnifCDF_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getUnifCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setUnifCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && IK5_ENABLED
    module procedure test_setUnifCDF_RK5_IK5
        use pm_kind, only: RKG => RK5, IKG => IK5
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && IK4_ENABLED
    module procedure test_setUnifCDF_RK5_IK4
        use pm_kind, only: RKG => RK5, IKG => IK4
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && IK3_ENABLED
    module procedure test_setUnifCDF_RK5_IK3
        use pm_kind, only: RKG => RK5, IKG => IK3
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && IK2_ENABLED
    module procedure test_setUnifCDF_RK5_IK2
        use pm_kind, only: RKG => RK5, IKG => IK2
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && IK1_ENABLED
    module procedure test_setUnifCDF_RK5_IK1
        use pm_kind, only: RKG => RK5, IKG => IK1
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK4_ENABLED && IK5_ENABLED
    module procedure test_setUnifCDF_RK4_IK5
        use pm_kind, only: RKG => RK4, IKG => IK5
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK4_ENABLED
    module procedure test_setUnifCDF_RK4_IK4
        use pm_kind, only: RKG => RK4, IKG => IK4
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK3_ENABLED
    module procedure test_setUnifCDF_RK4_IK3
        use pm_kind, only: RKG => RK4, IKG => IK3
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK2_ENABLED
    module procedure test_setUnifCDF_RK4_IK2
        use pm_kind, only: RKG => RK4, IKG => IK2
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && IK1_ENABLED
    module procedure test_setUnifCDF_RK4_IK1
        use pm_kind, only: RKG => RK4, IKG => IK1
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED && IK5_ENABLED
    module procedure test_setUnifCDF_RK3_IK5
        use pm_kind, only: RKG => RK3, IKG => IK5
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK4_ENABLED
    module procedure test_setUnifCDF_RK3_IK4
        use pm_kind, only: RKG => RK3, IKG => IK4
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK3_ENABLED
    module procedure test_setUnifCDF_RK3_IK3
        use pm_kind, only: RKG => RK3, IKG => IK3
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK2_ENABLED
    module procedure test_setUnifCDF_RK3_IK2
        use pm_kind, only: RKG => RK3, IKG => IK2
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && IK1_ENABLED
    module procedure test_setUnifCDF_RK3_IK1
        use pm_kind, only: RKG => RK3, IKG => IK1
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK2_ENABLED && IK5_ENABLED
    module procedure test_setUnifCDF_RK2_IK5
        use pm_kind, only: RKG => RK2, IKG => IK5
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK4_ENABLED
    module procedure test_setUnifCDF_RK2_IK4
        use pm_kind, only: RKG => RK2, IKG => IK4
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK3_ENABLED
    module procedure test_setUnifCDF_RK2_IK3
        use pm_kind, only: RKG => RK2, IKG => IK3
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK2_ENABLED
    module procedure test_setUnifCDF_RK2_IK2
        use pm_kind, only: RKG => RK2, IKG => IK2
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && IK1_ENABLED
    module procedure test_setUnifCDF_RK2_IK1
        use pm_kind, only: RKG => RK2, IKG => IK1
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK1_ENABLED && IK5_ENABLED
    module procedure test_setUnifCDF_RK1_IK5
        use pm_kind, only: RKG => RK1, IKG => IK5
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK4_ENABLED
    module procedure test_setUnifCDF_RK1_IK4
        use pm_kind, only: RKG => RK1, IKG => IK4
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK3_ENABLED
    module procedure test_setUnifCDF_RK1_IK3
        use pm_kind, only: RKG => RK1, IKG => IK3
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK2_ENABLED
    module procedure test_setUnifCDF_RK1_IK2
        use pm_kind, only: RKG => RK1, IKG => IK2
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && IK1_ENABLED
    module procedure test_setUnifCDF_RK1_IK1
        use pm_kind, only: RKG => RK1, IKG => IK1
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module procedure test_setUnifCDF_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setUnifCDF_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setUnifCDF_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setUnifCDF_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setUnifCDF_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module procedure test_setUnifCDF_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setUnifCDF_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setUnifCDF_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setUnifCDF_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setUnifCDF_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distUnif@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setUnifCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getUnifRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getUnifRand_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getUnifRand_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getUnifRand_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getUnifRand_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getUnifRand_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getUnifRand_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getUnifRand_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getUnifRand_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getUnifRand_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getUnifRand_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getUnifRand_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getUnifRand_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getUnifRand_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getUnifRand_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getUnifRand_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getUnifRand_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getUnifRand_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getUnifRand_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getUnifRand_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getUnifRand_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getUnifRand_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getUnifRand_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getUnifRand_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getUnifRand_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getUnifRand_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distUnif@routines@getUnifRand.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getUnifRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setUnifRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setUnifRand_SK5
        use pm_kind, only: SKG => SK5
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setUnifRand_SK4
        use pm_kind, only: SKG => SK4
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setUnifRand_SK3
        use pm_kind, only: SKG => SK3
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setUnifRand_SK2
        use pm_kind, only: SKG => SK2
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setUnifRand_SK1
        use pm_kind, only: SKG => SK1
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setUnifRand_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setUnifRand_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setUnifRand_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setUnifRand_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setUnifRand_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setUnifRand_LK5
        use pm_kind, only: LKG => LK5
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setUnifRand_LK4
        use pm_kind, only: LKG => LK4
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setUnifRand_LK3
        use pm_kind, only: LKG => LK3
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setUnifRand_LK2
        use pm_kind, only: LKG => LK2
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setUnifRand_LK1
        use pm_kind, only: LKG => LK1
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setUnifRand_CK5
        use pm_kind, only: CKG => CK5
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setUnifRand_CK4
        use pm_kind, only: CKG => CK4
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setUnifRand_CK3
        use pm_kind, only: CKG => CK3
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setUnifRand_CK2
        use pm_kind, only: CKG => CK2
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setUnifRand_CK1
        use pm_kind, only: CKG => CK1
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setUnifRand_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setUnifRand_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setUnifRand_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setUnifRand_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setUnifRand_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distUnif@routines@setUnifRand.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setUnifRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE