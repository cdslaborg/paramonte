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
!>  This file contains procedure implementations of [pm_val2Real](@ref pm_val2real).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_val2real) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getReal_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRealDef_SK5
        use pm_kind, only: RKG => RK, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRealDef_SK4
        use pm_kind, only: RKG => RK, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRealDef_SK3
        use pm_kind, only: RKG => RK, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRealDef_SK2
        use pm_kind, only: RKG => RK, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRealDef_SK1
        use pm_kind, only: RKG => RK, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

#if LK5_ENABLED
    module procedure getRealDef_LK5
        use pm_kind, only: RKG => RK, LKG => LK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if LK4_ENABLED
    module procedure getRealDef_LK4
        use pm_kind, only: RKG => RK, LKG => LK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if LK3_ENABLED
    module procedure getRealDef_LK3
        use pm_kind, only: RKG => RK, LKG => LK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if LK2_ENABLED
    module procedure getRealDef_LK2
        use pm_kind, only: RKG => RK, LKG => LK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if LK1_ENABLED
    module procedure getRealDef_LK1
        use pm_kind, only: RKG => RK, LKG => LK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getReal_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setReal_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Err_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && SK5_ENABLED
    module procedure setRealErr_RK5_SK5
        use pm_kind, only: RKG => RK5, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && SK4_ENABLED
    module procedure setRealErr_RK5_SK4
        use pm_kind, only: RKG => RK5, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && SK3_ENABLED
    module procedure setRealErr_RK5_SK3
        use pm_kind, only: RKG => RK5, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && SK2_ENABLED
    module procedure setRealErr_RK5_SK2
        use pm_kind, only: RKG => RK5, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && SK1_ENABLED
    module procedure setRealErr_RK5_SK1
        use pm_kind, only: RKG => RK5, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK4_ENABLED && SK5_ENABLED
    module procedure setRealErr_RK4_SK5
        use pm_kind, only: RKG => RK4, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && SK4_ENABLED
    module procedure setRealErr_RK4_SK4
        use pm_kind, only: RKG => RK4, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && SK3_ENABLED
    module procedure setRealErr_RK4_SK3
        use pm_kind, only: RKG => RK4, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && SK2_ENABLED
    module procedure setRealErr_RK4_SK2
        use pm_kind, only: RKG => RK4, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && SK1_ENABLED
    module procedure setRealErr_RK4_SK1
        use pm_kind, only: RKG => RK4, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED && SK5_ENABLED
    module procedure setRealErr_RK3_SK5
        use pm_kind, only: RKG => RK3, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && SK4_ENABLED
    module procedure setRealErr_RK3_SK4
        use pm_kind, only: RKG => RK3, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && SK3_ENABLED
    module procedure setRealErr_RK3_SK3
        use pm_kind, only: RKG => RK3, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && SK2_ENABLED
    module procedure setRealErr_RK3_SK2
        use pm_kind, only: RKG => RK3, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && SK1_ENABLED
    module procedure setRealErr_RK3_SK1
        use pm_kind, only: RKG => RK3, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK2_ENABLED && SK5_ENABLED
    module procedure setRealErr_RK2_SK5
        use pm_kind, only: RKG => RK2, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && SK4_ENABLED
    module procedure setRealErr_RK2_SK4
        use pm_kind, only: RKG => RK2, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && SK3_ENABLED
    module procedure setRealErr_RK2_SK3
        use pm_kind, only: RKG => RK2, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && SK2_ENABLED
    module procedure setRealErr_RK2_SK2
        use pm_kind, only: RKG => RK2, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && SK1_ENABLED
    module procedure setRealErr_RK2_SK1
        use pm_kind, only: RKG => RK2, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK1_ENABLED && SK5_ENABLED
    module procedure setRealErr_RK1_SK5
        use pm_kind, only: RKG => RK1, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && SK4_ENABLED
    module procedure setRealErr_RK1_SK4
        use pm_kind, only: RKG => RK1, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && SK3_ENABLED
    module procedure setRealErr_RK1_SK3
        use pm_kind, only: RKG => RK1, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && SK2_ENABLED
    module procedure setRealErr_RK1_SK2
        use pm_kind, only: RKG => RK1, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && SK1_ENABLED
    module procedure setRealErr_RK1_SK1
        use pm_kind, only: RKG => RK1, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Err_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && SK5_ENABLED
    module procedure setRealDef_RK5_SK5
        use pm_kind, only: RKG => RK5, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && SK4_ENABLED
    module procedure setRealDef_RK5_SK4
        use pm_kind, only: RKG => RK5, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && SK3_ENABLED
    module procedure setRealDef_RK5_SK3
        use pm_kind, only: RKG => RK5, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && SK2_ENABLED
    module procedure setRealDef_RK5_SK2
        use pm_kind, only: RKG => RK5, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && SK1_ENABLED
    module procedure setRealDef_RK5_SK1
        use pm_kind, only: RKG => RK5, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK4_ENABLED && SK5_ENABLED
    module procedure setRealDef_RK4_SK5
        use pm_kind, only: RKG => RK4, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && SK4_ENABLED
    module procedure setRealDef_RK4_SK4
        use pm_kind, only: RKG => RK4, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && SK3_ENABLED
    module procedure setRealDef_RK4_SK3
        use pm_kind, only: RKG => RK4, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && SK2_ENABLED
    module procedure setRealDef_RK4_SK2
        use pm_kind, only: RKG => RK4, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && SK1_ENABLED
    module procedure setRealDef_RK4_SK1
        use pm_kind, only: RKG => RK4, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED && SK5_ENABLED
    module procedure setRealDef_RK3_SK5
        use pm_kind, only: RKG => RK3, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && SK4_ENABLED
    module procedure setRealDef_RK3_SK4
        use pm_kind, only: RKG => RK3, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && SK3_ENABLED
    module procedure setRealDef_RK3_SK3
        use pm_kind, only: RKG => RK3, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && SK2_ENABLED
    module procedure setRealDef_RK3_SK2
        use pm_kind, only: RKG => RK3, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && SK1_ENABLED
    module procedure setRealDef_RK3_SK1
        use pm_kind, only: RKG => RK3, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK2_ENABLED && SK5_ENABLED
    module procedure setRealDef_RK2_SK5
        use pm_kind, only: RKG => RK2, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && SK4_ENABLED
    module procedure setRealDef_RK2_SK4
        use pm_kind, only: RKG => RK2, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && SK3_ENABLED
    module procedure setRealDef_RK2_SK3
        use pm_kind, only: RKG => RK2, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && SK2_ENABLED
    module procedure setRealDef_RK2_SK2
        use pm_kind, only: RKG => RK2, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && SK1_ENABLED
    module procedure setRealDef_RK2_SK1
        use pm_kind, only: RKG => RK2, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK1_ENABLED && SK5_ENABLED
    module procedure setRealDef_RK1_SK5
        use pm_kind, only: RKG => RK1, SKG => SK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && SK4_ENABLED
    module procedure setRealDef_RK1_SK4
        use pm_kind, only: RKG => RK1, SKG => SK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && SK3_ENABLED
    module procedure setRealDef_RK1_SK3
        use pm_kind, only: RKG => RK1, SKG => SK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && SK2_ENABLED
    module procedure setRealDef_RK1_SK2
        use pm_kind, only: RKG => RK1, SKG => SK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && SK1_ENABLED
    module procedure setRealDef_RK1_SK1
        use pm_kind, only: RKG => RK1, SKG => SK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED && LK5_ENABLED
    module procedure setRealDef_RK5_LK5
        use pm_kind, only: RKG => RK5, LKG => LK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && LK4_ENABLED
    module procedure setRealDef_RK5_LK4
        use pm_kind, only: RKG => RK5, LKG => LK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && LK3_ENABLED
    module procedure setRealDef_RK5_LK3
        use pm_kind, only: RKG => RK5, LKG => LK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && LK2_ENABLED
    module procedure setRealDef_RK5_LK2
        use pm_kind, only: RKG => RK5, LKG => LK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK5_ENABLED && LK1_ENABLED
    module procedure setRealDef_RK5_LK1
        use pm_kind, only: RKG => RK5, LKG => LK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK4_ENABLED && LK5_ENABLED
    module procedure setRealDef_RK4_LK5
        use pm_kind, only: RKG => RK4, LKG => LK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && LK4_ENABLED
    module procedure setRealDef_RK4_LK4
        use pm_kind, only: RKG => RK4, LKG => LK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && LK3_ENABLED
    module procedure setRealDef_RK4_LK3
        use pm_kind, only: RKG => RK4, LKG => LK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && LK2_ENABLED
    module procedure setRealDef_RK4_LK2
        use pm_kind, only: RKG => RK4, LKG => LK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED && LK1_ENABLED
    module procedure setRealDef_RK4_LK1
        use pm_kind, only: RKG => RK4, LKG => LK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK3_ENABLED && LK5_ENABLED
    module procedure setRealDef_RK3_LK5
        use pm_kind, only: RKG => RK3, LKG => LK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && LK4_ENABLED
    module procedure setRealDef_RK3_LK4
        use pm_kind, only: RKG => RK3, LKG => LK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && LK3_ENABLED
    module procedure setRealDef_RK3_LK3
        use pm_kind, only: RKG => RK3, LKG => LK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && LK2_ENABLED
    module procedure setRealDef_RK3_LK2
        use pm_kind, only: RKG => RK3, LKG => LK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED && LK1_ENABLED
    module procedure setRealDef_RK3_LK1
        use pm_kind, only: RKG => RK3, LKG => LK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK2_ENABLED && LK5_ENABLED
    module procedure setRealDef_RK2_LK5
        use pm_kind, only: RKG => RK2, LKG => LK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && LK4_ENABLED
    module procedure setRealDef_RK2_LK4
        use pm_kind, only: RKG => RK2, LKG => LK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && LK3_ENABLED
    module procedure setRealDef_RK2_LK3
        use pm_kind, only: RKG => RK2, LKG => LK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && LK2_ENABLED
    module procedure setRealDef_RK2_LK2
        use pm_kind, only: RKG => RK2, LKG => LK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED && LK1_ENABLED
    module procedure setRealDef_RK2_LK1
        use pm_kind, only: RKG => RK2, LKG => LK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK1_ENABLED && LK5_ENABLED
    module procedure setRealDef_RK1_LK5
        use pm_kind, only: RKG => RK1, LKG => LK5
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && LK4_ENABLED
    module procedure setRealDef_RK1_LK4
        use pm_kind, only: RKG => RK1, LKG => LK4
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && LK3_ENABLED
    module procedure setRealDef_RK1_LK3
        use pm_kind, only: RKG => RK1, LKG => LK3
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && LK2_ENABLED
    module procedure setRealDef_RK1_LK2
        use pm_kind, only: RKG => RK1, LKG => LK2
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED && LK1_ENABLED
    module procedure setRealDef_RK1_LK1
        use pm_kind, only: RKG => RK1, LKG => LK1
#include "pm_val2real@routines.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setReal_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines