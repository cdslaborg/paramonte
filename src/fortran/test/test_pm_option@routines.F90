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
!>  This file contains procedure implementations of [test_pm_option](@ref test_pm_option).
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_option) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define option_CK_ENABLED 1

#if CK3_ENABLED
    module procedure test_getOption_CK3_1
        use pm_kind, only: IK, LK, CK => CK3
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getOption_CK2_1
        use pm_kind, only: IK, LK, CK => CK2
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getOption_CK1_1
        use pm_kind, only: IK, LK, CK => CK1
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#if CK5_ENABLED
    module procedure test_getOption_CK5_1
        use pm_kind, only: IK, LK, CK => CK5
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#undef option_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define option_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getOption_RK5_1
        use pm_kind, only: IK, LK, RK => RK5
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getOption_RK4_1
        use pm_kind, only: IK, LK, RK => RK4
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getOption_RK3_1
        use pm_kind, only: IK, LK, RK => RK3
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getOption_RK2_1
        use pm_kind, only: IK, LK, RK => RK2
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getOption_RK1_1
        use pm_kind, only: IK, LK, RK => RK1
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#undef option_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define option_IK_ENABLED 1

#if IK4_ENABLED
    module procedure test_getOption_IK4_1
        use pm_kind, only: IK, LK, IKC => IK4
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getOption_IK3_1
        use pm_kind, only: IK, LK, IKC => IK3
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getOption_IK2_1
        use pm_kind, only: IK, LK, IKC => IK2
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getOption_IK1_1
        use pm_kind, only: IK, LK, IKC => IK1
#include "test_pm_option@routines.inc.F90"
    end procedure
#endif

#undef option_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define option_LK_ENABLED 1

    module procedure test_getOption_LK_1
        use pm_kind, only: IK, LK
#include "test_pm_option@routines.inc.F90"
    end procedure

#undef option_LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define option_SK_ENABLED 1

    module procedure test_getOption_SK_1
        use pm_kind, only: SK, IK, LK
#include "test_pm_option@routines.inc.F90"
    end procedure

#undef option_SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
