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
!>  This file contains the implementations of the tests of module [pm_arrayChange](@ref pm_arrayChange).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_arrayChange) routines

    use pm_kind, only: LK, SK
    use pm_arrayRange, only: getRange
    use pm_arrayUnique, only: getUnique
    use pm_arrayResize, only: setResized
    use pm_arrayUnique, only: isUniqueAll
    use pm_arrayVerbose, only: getVerbose
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_arrayMembership, only: operator(.in.)
    use pm_arrayMembership, only: operator(.allin.)
    use pm_io, only: display_type
    use pm_distUnif, only: rngf
    use pm_option, only: getOption
    use pm_val2str, only: getStr
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getChange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getChange_SK5
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getChange_SK4
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getChange_SK3
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getChange_SK2
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getChange_SK1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayChange@routines.inc.F90"
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

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getChange_IK5
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getChange_IK4
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getChange_IK3
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getChange_IK2
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getChange_IK1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getChange_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getChange_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getChange_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getChange_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getChange_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getChange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setChange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setChange_SK5
        use pm_kind, only: SKC => SK5
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setChange_SK4
        use pm_kind, only: SKC => SK4
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setChange_SK3
        use pm_kind, only: SKC => SK3
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setChange_SK2
        use pm_kind, only: SKC => SK2
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setChange_SK1
        use pm_kind, only: SKC => SK1
#include "test_pm_arrayChange@routines.inc.F90"
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

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setChange_IK5
        use pm_kind, only: IKC => IK5
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setChange_IK4
        use pm_kind, only: IKC => IK4
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setChange_IK3
        use pm_kind, only: IKC => IK3
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setChange_IK2
        use pm_kind, only: IKC => IK2
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setChange_IK1
        use pm_kind, only: IKC => IK1
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setChange_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setChange_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setChange_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setChange_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setChange_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_arrayChange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setChange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE