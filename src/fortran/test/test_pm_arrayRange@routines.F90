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
!>  This file contains the implementations of the tests of module [pm_arrayRange](@ref pm_arrayRange).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_arrayRange) routines

    use pm_kind, only: LK
    use pm_val2str, only: getStr
    use pm_mathCompare, only: isClose
    use pm_io, only: display_type
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRange_SK5
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRange_SK4
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRange_SK3
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRange_SK2
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRange_SK1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getRange_IK5
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getRange_IK4
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getRange_IK3
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getRange_IK2
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getRange_IK1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRange_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRange_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRange_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRange_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRange_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRange_SK5
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRange_SK4
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRange_SK3
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRange_SK2
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRange_SK1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setRange_IK5
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setRange_IK4
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setRange_IK3
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setRange_IK2
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setRange_IK1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRange_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRange_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRange_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRange_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRange_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
