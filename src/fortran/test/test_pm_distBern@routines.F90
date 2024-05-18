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
!>  This include file contains procedure implementations of the tests of [pm_distBern](@ref pm_distBern).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distBern) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isHead_ENABLED 1

#if RK5_ENABLED
    module procedure test_isHead_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_isHead_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_isHead_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_isHead_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_isHead_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#undef isHead_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBernRand_ENABLED 1

#if RK5_ENABLED
    module procedure test_getBernRand_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getBernRand_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getBernRand_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getBernRand_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getBernRand_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#undef getBernRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBernRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED && RK5_ENABLED
    module procedure test_setBernRand_IK5_RK5_1
        use pm_kind, only: IKC => IK5, RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && RK4_ENABLED
    module procedure test_setBernRand_IK5_RK4_1
        use pm_kind, only: IKC => IK5, RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && RK3_ENABLED
    module procedure test_setBernRand_IK5_RK3_1
        use pm_kind, only: IKC => IK5, RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && RK2_ENABLED
    module procedure test_setBernRand_IK5_RK2_1
        use pm_kind, only: IKC => IK5, RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK5_ENABLED && RK1_ENABLED
    module procedure test_setBernRand_IK5_RK1_1
        use pm_kind, only: IKC => IK5, RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED && RK5_ENABLED
    module procedure test_setBernRand_IK4_RK5_1
        use pm_kind, only: IKC => IK4, RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && RK4_ENABLED
    module procedure test_setBernRand_IK4_RK4_1
        use pm_kind, only: IKC => IK4, RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && RK3_ENABLED
    module procedure test_setBernRand_IK4_RK3_1
        use pm_kind, only: IKC => IK4, RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && RK2_ENABLED
    module procedure test_setBernRand_IK4_RK2_1
        use pm_kind, only: IKC => IK4, RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && RK1_ENABLED
    module procedure test_setBernRand_IK4_RK1_1
        use pm_kind, only: IKC => IK4, RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK3_ENABLED && RK5_ENABLED
    module procedure test_setBernRand_IK3_RK5_1
        use pm_kind, only: IKC => IK3, RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && RK4_ENABLED
    module procedure test_setBernRand_IK3_RK4_1
        use pm_kind, only: IKC => IK3, RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && RK3_ENABLED
    module procedure test_setBernRand_IK3_RK3_1
        use pm_kind, only: IKC => IK3, RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && RK2_ENABLED
    module procedure test_setBernRand_IK3_RK2_1
        use pm_kind, only: IKC => IK3, RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && RK1_ENABLED
    module procedure test_setBernRand_IK3_RK1_1
        use pm_kind, only: IKC => IK3, RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK2_ENABLED && RK5_ENABLED
    module procedure test_setBernRand_IK2_RK5_1
        use pm_kind, only: IKC => IK2, RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && RK4_ENABLED
    module procedure test_setBernRand_IK2_RK4_1
        use pm_kind, only: IKC => IK2, RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && RK3_ENABLED
    module procedure test_setBernRand_IK2_RK3_1
        use pm_kind, only: IKC => IK2, RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && RK2_ENABLED
    module procedure test_setBernRand_IK2_RK2_1
        use pm_kind, only: IKC => IK2, RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && RK1_ENABLED
    module procedure test_setBernRand_IK2_RK1_1
        use pm_kind, only: IKC => IK2, RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK1_ENABLED && RK5_ENABLED
    module procedure test_setBernRand_IK1_RK5_1
        use pm_kind, only: IKC => IK1, RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && RK4_ENABLED
    module procedure test_setBernRand_IK1_RK4_1
        use pm_kind, only: IKC => IK1, RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && RK3_ENABLED
    module procedure test_setBernRand_IK1_RK3_1
        use pm_kind, only: IKC => IK1, RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && RK2_ENABLED
    module procedure test_setBernRand_IK1_RK2_1
        use pm_kind, only: IKC => IK1, RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && RK1_ENABLED
    module procedure test_setBernRand_IK1_RK1_1
        use pm_kind, only: IKC => IK1, RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED && RK5_ENABLED
    module procedure test_setBernRand_LK5_RK5_1
        use pm_kind, only: LKC => LK5, RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK5_ENABLED && RK4_ENABLED
    module procedure test_setBernRand_LK5_RK4_1
        use pm_kind, only: LKC => LK5, RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK5_ENABLED && RK3_ENABLED
    module procedure test_setBernRand_LK5_RK3_1
        use pm_kind, only: LKC => LK5, RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK5_ENABLED && RK2_ENABLED
    module procedure test_setBernRand_LK5_RK2_1
        use pm_kind, only: LKC => LK5, RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK5_ENABLED && RK1_ENABLED
    module procedure test_setBernRand_LK5_RK1_1
        use pm_kind, only: LKC => LK5, RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK4_ENABLED && RK5_ENABLED
    module procedure test_setBernRand_LK4_RK5_1
        use pm_kind, only: LKC => LK4, RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED && RK4_ENABLED
    module procedure test_setBernRand_LK4_RK4_1
        use pm_kind, only: LKC => LK4, RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED && RK3_ENABLED
    module procedure test_setBernRand_LK4_RK3_1
        use pm_kind, only: LKC => LK4, RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED && RK2_ENABLED
    module procedure test_setBernRand_LK4_RK2_1
        use pm_kind, only: LKC => LK4, RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED && RK1_ENABLED
    module procedure test_setBernRand_LK4_RK1_1
        use pm_kind, only: LKC => LK4, RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK3_ENABLED && RK5_ENABLED
    module procedure test_setBernRand_LK3_RK5_1
        use pm_kind, only: LKC => LK3, RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED && RK4_ENABLED
    module procedure test_setBernRand_LK3_RK4_1
        use pm_kind, only: LKC => LK3, RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED && RK3_ENABLED
    module procedure test_setBernRand_LK3_RK3_1
        use pm_kind, only: LKC => LK3, RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED && RK2_ENABLED
    module procedure test_setBernRand_LK3_RK2_1
        use pm_kind, only: LKC => LK3, RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED && RK1_ENABLED
    module procedure test_setBernRand_LK3_RK1_1
        use pm_kind, only: LKC => LK3, RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK2_ENABLED && RK5_ENABLED
    module procedure test_setBernRand_LK2_RK5_1
        use pm_kind, only: LKC => LK2, RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED && RK4_ENABLED
    module procedure test_setBernRand_LK2_RK4_1
        use pm_kind, only: LKC => LK2, RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED && RK3_ENABLED
    module procedure test_setBernRand_LK2_RK3_1
        use pm_kind, only: LKC => LK2, RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED && RK2_ENABLED
    module procedure test_setBernRand_LK2_RK2_1
        use pm_kind, only: LKC => LK2, RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED && RK1_ENABLED
    module procedure test_setBernRand_LK2_RK1_1
        use pm_kind, only: LKC => LK2, RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK1_ENABLED && RK5_ENABLED
    module procedure test_setBernRand_LK1_RK5_1
        use pm_kind, only: LKC => LK1, RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED && RK4_ENABLED
    module procedure test_setBernRand_LK1_RK4_1
        use pm_kind, only: LKC => LK1, RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED && RK3_ENABLED
    module procedure test_setBernRand_LK1_RK3_1
        use pm_kind, only: LKC => LK1, RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED && RK2_ENABLED
    module procedure test_setBernRand_LK1_RK2_1
        use pm_kind, only: LKC => LK1, RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED && RK1_ENABLED
    module procedure test_setBernRand_LK1_RK1_1
        use pm_kind, only: LKC => LK1, RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module procedure test_setBernRand_RK5_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setBernRand_RK4_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setBernRand_RK3_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setBernRand_RK2_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setBernRand_RK1_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distBern@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBernRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
