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
!>  This file contains procedure implementations of [pm_optimization](@ref pm_optimization).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_optimization) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_matrixInit, only: setMatInit, uppLowDia
    use pm_mathConst, only: GOLDEN_RATIO
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isBracketMax_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isBracketMax_RK5
        use pm_kind, only: RKC => RK5
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isBracketMax_RK4
        use pm_kind, only: RKC => RK4
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isBracketMax_RK3
        use pm_kind, only: RKC => RK3
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isBracketMax_RK2
        use pm_kind, only: RKC => RK2
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isBracketMax_RK1
        use pm_kind, only: RKC => RK1
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isBracketMax_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isBracketMin_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isBracketMin_RK5
        use pm_kind, only: RKC => RK5
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isBracketMin_RK4
        use pm_kind, only: RKC => RK4
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isBracketMin_RK3
        use pm_kind, only: RKC => RK3
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isBracketMin_RK2
        use pm_kind, only: RKC => RK2
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isBracketMin_RK1
        use pm_kind, only: RKC => RK1
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isBracketMin_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBracketMax_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBracketMax_RK5
        use pm_kind, only: RKC => RK5
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBracketMax_RK4
        use pm_kind, only: RKC => RK4
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBracketMax_RK3
        use pm_kind, only: RKC => RK3
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBracketMax_RK2
        use pm_kind, only: RKC => RK2
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBracketMax_RK1
        use pm_kind, only: RKC => RK1
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBracketMax_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBracketMin_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBracketMin_RK5
        use pm_kind, only: RKC => RK5
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBracketMin_RK4
        use pm_kind, only: RKC => RK4
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBracketMin_RK3
        use pm_kind, only: RKC => RK3
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBracketMin_RK2
        use pm_kind, only: RKC => RK2
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBracketMin_RK1
        use pm_kind, only: RKC => RK1
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBracketMin_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMinBrent_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMinBrent_RK5
        use pm_kind, only: RKC => RK5
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMinBrent_RK4
        use pm_kind, only: RKC => RK4
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMinBrent_RK3
        use pm_kind, only: RKC => RK3
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMinBrent_RK2
        use pm_kind, only: RKC => RK2
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMinBrent_RK1
        use pm_kind, only: RKC => RK1
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMinBrent_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMinBrent_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMinBrent_RK5
        use pm_kind, only: RKC => RK5
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMinBrent_RK4
        use pm_kind, only: RKC => RK4
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMinBrent_RK3
        use pm_kind, only: RKC => RK3
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMinBrent_RK2
        use pm_kind, only: RKC => RK2
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMinBrent_RK1
        use pm_kind, only: RKC => RK1
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMinBrent_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isFailedMinPowell_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure isFailedMinPowell_RK5
        use pm_kind, only: RKC => RK5
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure isFailedMinPowell_RK4
        use pm_kind, only: RKC => RK4
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure isFailedMinPowell_RK3
        use pm_kind, only: RKC => RK3
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure isFailedMinPowell_RK2
        use pm_kind, only: RKC => RK2
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure isFailedMinPowell_RK1
        use pm_kind, only: RKC => RK1
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isFailedMinPowell_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMinPowell_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMinPowell_RK5
        use pm_kind, only: RKC => RK5
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMinPowell_RK4
        use pm_kind, only: RKC => RK4
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMinPowell_RK3
        use pm_kind, only: RKC => RK3
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMinPowell_RK2
        use pm_kind, only: RKC => RK2
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMinPowell_RK1
        use pm_kind, only: RKC => RK1
#include "pm_optimization@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMinPowell_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines