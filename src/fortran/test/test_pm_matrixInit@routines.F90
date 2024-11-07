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
!>  This file contains test implementations of the procedure of [pm_matrixInit](@ref pm_matrixInit).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_matrixInit) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatInitULD_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatInitULD_D2_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getMatInitULD_D2_SK5_1
        use pm_kind, only: SKG => SK5
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getMatInitULD_D2_SK4_1
        use pm_kind, only: SKG => SK4
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getMatInitULD_D2_SK3_1
        use pm_kind, only: SKG => SK3
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getMatInitULD_D2_SK2_1
        use pm_kind, only: SKG => SK2
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getMatInitULD_D2_SK1_1
        use pm_kind, only: SKG => SK1
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef getMatInitULD_D2_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatInitULD_D2_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getMatInitULD_D2_IK5_1
        use pm_kind, only: IKG => IK5
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getMatInitULD_D2_IK4_1
        use pm_kind, only: IKG => IK4
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getMatInitULD_D2_IK3_1
        use pm_kind, only: IKG => IK3
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getMatInitULD_D2_IK2_1
        use pm_kind, only: IKG => IK2
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getMatInitULD_D2_IK1_1
        use pm_kind, only: IKG => IK1
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef getMatInitULD_D2_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatInitULD_D2_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_getMatInitULD_D2_LK5_1
        use pm_kind, only: LKG => LK5
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_getMatInitULD_D2_LK4_1
        use pm_kind, only: LKG => LK4
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_getMatInitULD_D2_LK3_1
        use pm_kind, only: LKG => LK3
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_getMatInitULD_D2_LK2_1
        use pm_kind, only: LKG => LK2
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_getMatInitULD_D2_LK1_1
        use pm_kind, only: LKG => LK1
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef getMatInitULD_D2_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatInitULD_D2_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getMatInitULD_D2_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getMatInitULD_D2_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getMatInitULD_D2_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getMatInitULD_D2_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getMatInitULD_D2_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef getMatInitULD_D2_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatInitULD_D2_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getMatInitULD_D2_RK5_1
        use pm_kind, only: RKG => RK5
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getMatInitULD_D2_RK4_1
        use pm_kind, only: RKG => RK4
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getMatInitULD_D2_RK3_1
        use pm_kind, only: RKG => RK3
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getMatInitULD_D2_RK2_1
        use pm_kind, only: RKG => RK2
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getMatInitULD_D2_RK1_1
        use pm_kind, only: RKG => RK1
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef getMatInitULD_D2_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatInitULD_D2_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#define setMatInitULD_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatInitULD_D2_SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setMatInitULD_D2_SK5_1
        use pm_kind, only: SKG => SK5
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setMatInitULD_D2_SK4_1
        use pm_kind, only: SKG => SK4
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setMatInitULD_D2_SK3_1
        use pm_kind, only: SKG => SK3
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setMatInitULD_D2_SK2_1
        use pm_kind, only: SKG => SK2
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setMatInitULD_D2_SK1_1
        use pm_kind, only: SKG => SK1
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef setMatInitULD_D2_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatInitULD_D2_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setMatInitULD_D2_IK5_1
        use pm_kind, only: IKG => IK5
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setMatInitULD_D2_IK4_1
        use pm_kind, only: IKG => IK4
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setMatInitULD_D2_IK3_1
        use pm_kind, only: IKG => IK3
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setMatInitULD_D2_IK2_1
        use pm_kind, only: IKG => IK2
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setMatInitULD_D2_IK1_1
        use pm_kind, only: IKG => IK1
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef setMatInitULD_D2_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatInitULD_D2_LK_ENABLED 1

#if LK5_ENABLED
    module procedure test_setMatInitULD_D2_LK5_1
        use pm_kind, only: LKG => LK5
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure test_setMatInitULD_D2_LK4_1
        use pm_kind, only: LKG => LK4
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure test_setMatInitULD_D2_LK3_1
        use pm_kind, only: LKG => LK3
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure test_setMatInitULD_D2_LK2_1
        use pm_kind, only: LKG => LK2
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure test_setMatInitULD_D2_LK1_1
        use pm_kind, only: LKG => LK1
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef setMatInitULD_D2_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatInitULD_D2_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setMatInitULD_D2_CK5_1
        use pm_kind, only: CKG => CK5
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setMatInitULD_D2_CK4_1
        use pm_kind, only: CKG => CK4
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setMatInitULD_D2_CK3_1
        use pm_kind, only: CKG => CK3
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setMatInitULD_D2_CK2_1
        use pm_kind, only: CKG => CK2
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setMatInitULD_D2_CK1_1
        use pm_kind, only: CKG => CK1
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef setMatInitULD_D2_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatInitULD_D2_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setMatInitULD_D2_RK5_1
        use pm_kind, only: RKG => RK5
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setMatInitULD_D2_RK4_1
        use pm_kind, only: RKG => RK4
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setMatInitULD_D2_RK3_1
        use pm_kind, only: RKG => RK3
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setMatInitULD_D2_RK2_1
        use pm_kind, only: RKG => RK2
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setMatInitULD_D2_RK1_1
        use pm_kind, only: RKG => RK1
#include "test_pm_matrixInit@routines.inc.F90"
    end procedure
#endif

#undef setMatInitULD_D2_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatInitULD_D2_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
