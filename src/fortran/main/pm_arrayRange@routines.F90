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
!>  This file contains procedure implementations of [pm_arrayRange](@ref pm_arrayRange).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_arrayRange) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arrayResize, only: setResized
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRangeUnit_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRangeUnit_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRangeUnit_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRangeUnit_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRangeUnit_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Step_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRangeStep_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRangeStep_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRangeStep_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRangeStep_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRangeStep_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Step_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRangeUnit_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRangeUnit_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRangeUnit_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRangeUnit_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRangeUnit_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRangeUnit_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRangeUnit_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRangeUnit_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRangeUnit_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRangeUnit_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Step_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRangeStep_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRangeStep_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRangeStep_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRangeStep_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRangeStep_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRangeStep_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRangeStep_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRangeStep_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRangeStep_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRangeStep_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Step_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#define D2_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define Unit_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure getRangeUnit_D2_IK5
!        use pm_kind, only: IKC => IK5
!#include "pm_arrayRange@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure getRangeUnit_D2_IK4
!        use pm_kind, only: IKC => IK4
!#include "pm_arrayRange@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure getRangeUnit_D2_IK3
!        use pm_kind, only: IKC => IK3
!#include "pm_arrayRange@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure getRangeUnit_D2_IK2
!        use pm_kind, only: IKC => IK2
!#include "pm_arrayRange@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure getRangeUnit_D2_IK1
!        use pm_kind, only: IKC => IK1
!#include "pm_arrayRange@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef Unit_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define Step_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure getRangeStep_D2_IK5
!        use pm_kind, only: IKC => IK5
!#include "pm_arrayRange@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure getRangeStep_D2_IK4
!        use pm_kind, only: IKC => IK4
!#include "pm_arrayRange@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure getRangeStep_D2_IK3
!        use pm_kind, only: IKC => IK3
!#include "pm_arrayRange@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure getRangeStep_D2_IK2
!        use pm_kind, only: IKC => IK2
!#include "pm_arrayRange@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure getRangeStep_D2_IK1
!        use pm_kind, only: IKC => IK1
!#include "pm_arrayRange@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef Step_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef D2_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRange_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRangeUnit_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRangeUnit_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRangeUnit_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRangeUnit_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRangeUnit_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Step_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRangeStep_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRangeStep_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRangeStep_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRangeStep_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRangeStep_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Step_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Unit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRangeUnit_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRangeUnit_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRangeUnit_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRangeUnit_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRangeUnit_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRangeUnit_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRangeUnit_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRangeUnit_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRangeUnit_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRangeUnit_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Unit_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Step_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRangeStep_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRangeStep_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRangeStep_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRangeStep_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRangeStep_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRangeStep_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRangeStep_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRangeStep_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRangeStep_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRangeStep_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_arrayRange@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Step_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRange_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines