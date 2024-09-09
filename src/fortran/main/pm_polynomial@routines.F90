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
!>  This file contains procedure implementations of [pm_polynomial](@ref pm_polynomial).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_polynomial) routines ! LCOV_EXCL_LINE

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
    use pm_str, only: getTrimmedTZ
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPolyVal_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyVal_D0_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyVal_D0_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyVal_D0_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyVal_D0_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyVal_D0_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyVal_D0_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyVal_D0_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyVal_D0_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyVal_D0_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyVal_D0_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyVal_D0_RK5_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyVal_D0_RK4_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyVal_D0_RK3_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyVal_D0_RK2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyVal_D0_RK1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyVal_D1_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyVal_D1_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyVal_D1_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyVal_D1_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyVal_D1_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyVal_D1_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyVal_D1_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyVal_D1_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyVal_D1_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyVal_D1_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyVal_D1_RK5_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyVal_D1_RK4_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyVal_D1_RK3_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyVal_D1_RK2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyVal_D1_RK1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPolyVal_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPolyAdd_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyAdd_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyAdd_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyAdd_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyAdd_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyAdd_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyAdd_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyAdd_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyAdd_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyAdd_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyAdd_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPolyAdd_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPolyAdd_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolyAdd_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolyAdd_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolyAdd_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolyAdd_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolyAdd_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyAdd_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyAdd_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyAdd_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyAdd_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyAdd_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPolyAdd_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPolySub_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolySub_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolySub_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolySub_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolySub_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolySub_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolySub_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolySub_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolySub_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolySub_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolySub_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPolySub_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPolySub_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolySub_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolySub_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolySub_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolySub_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolySub_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolySub_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolySub_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolySub_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolySub_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolySub_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPolySub_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPolyMul_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyMul_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyMul_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyMul_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyMul_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyMul_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyMul_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyMul_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyMul_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyMul_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyMul_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPolyMul_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPolyMul_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolyMul_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolyMul_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolyMul_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolyMul_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolyMul_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyMul_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyMul_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyMul_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyMul_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyMul_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPolyMul_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPolyDiv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolyDiv_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolyDiv_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolyDiv_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolyDiv_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolyDiv_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyDiv_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyDiv_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyDiv_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyDiv_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyDiv_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPolyDiv_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPolyDiff_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyDiffDef_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyDiffDef_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyDiffDef_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyDiffDef_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyDiffDef_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyDiffDef_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyDiffDef_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyDiffDef_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyDiffDef_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyDiffDef_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ord_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyDiffOrd_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyDiffOrd_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyDiffOrd_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyDiffOrd_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyDiffOrd_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyDiffOrd_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyDiffOrd_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyDiffOrd_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyDiffOrd_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyDiffOrd_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ord_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPolyDiff_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPolyDiff_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolyDiffDef_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolyDiffDef_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolyDiffDef_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolyDiffDef_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolyDiffDef_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyDiffDef_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyDiffDef_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyDiffDef_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyDiffDef_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyDiffDef_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ord_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolyDiffOrd_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolyDiffOrd_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolyDiffOrd_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolyDiffOrd_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolyDiffOrd_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyDiffOrd_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyDiffOrd_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyDiffOrd_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyDiffOrd_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyDiffOrd_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ord_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPolyDiff_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPolyStr_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyStrDef_CK5
        use pm_kind, only: SKG => SK, TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyStrDef_CK4
        use pm_kind, only: SKG => SK, TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyStrDef_CK3
        use pm_kind, only: SKG => SK, TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyStrDef_CK2
        use pm_kind, only: SKG => SK, TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyStrDef_CK1
        use pm_kind, only: SKG => SK, TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyStrDef_RK5
        use pm_kind, only: SKG => SK, TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyStrDef_RK4
        use pm_kind, only: SKG => SK, TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyStrDef_RK3
        use pm_kind, only: SKG => SK, TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyStrDef_RK2
        use pm_kind, only: SKG => SK, TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyStrDef_RK1
        use pm_kind, only: SKG => SK, TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPolyStr_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPolyRoot_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyRootDef_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyRootDef_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyRootDef_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyRootDef_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyRootDef_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyRootDef_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyRootDef_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyRootDef_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyRootDef_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyRootDef_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Eig_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyRootEig_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyRootEig_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyRootEig_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyRootEig_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyRootEig_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyRootEig_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyRootEig_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyRootEig_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyRootEig_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyRootEig_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Eig_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Jen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyRootJen_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyRootJen_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyRootJen_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyRootJen_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyRootJen_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyRootJen_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyRootJen_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyRootJen_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyRootJen_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyRootJen_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Jen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Lag_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyRootLag_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyRootLag_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyRootLag_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyRootLag_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyRootLag_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyRootLag_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyRootLag_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyRootLag_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyRootLag_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyRootLag_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Lag_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SGL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure getPolyRootSGL_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPolyRootSGL_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPolyRootSGL_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPolyRootSGL_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPolyRootSGL_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure getPolyRootSGL_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPolyRootSGL_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPolyRootSGL_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPolyRootSGL_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPolyRootSGL_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SGL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPolyRoot_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPolyRoot_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Eig_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolyRootEig_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolyRootEig_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolyRootEig_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolyRootEig_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolyRootEig_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyRootEig_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyRootEig_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyRootEig_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyRootEig_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyRootEig_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Eig_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Jen_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolyRootJen_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolyRootJen_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolyRootJen_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolyRootJen_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolyRootJen_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyRootJen_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyRootJen_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyRootJen_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyRootJen_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyRootJen_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Jen_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Lag_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolyRootLag_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolyRootLag_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolyRootLag_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolyRootLag_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolyRootLag_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyRootLag_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyRootLag_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyRootLag_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyRootLag_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyRootLag_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Lag_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SGL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolyRootSGL_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolyRootSGL_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolyRootSGL_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolyRootSGL_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolyRootSGL_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyRootSGL_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyRootSGL_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyRootSGL_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyRootSGL_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyRootSGL_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SGL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPolyRoot_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPolyRootPolished_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Lag_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_CK_ENABLED 1

#if CK5_ENABLED
    module procedure setPolyRootPolishedLag_CK5_CK5
        use pm_kind, only: TKG => CK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPolyRootPolishedLag_CK4_CK4
        use pm_kind, only: TKG => CK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPolyRootPolishedLag_CK3_CK3
        use pm_kind, only: TKG => CK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPolyRootPolishedLag_CK2_CK2
        use pm_kind, only: TKG => CK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPolyRootPolishedLag_CK1_CK1
        use pm_kind, only: TKG => CK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef CK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_CK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyRootPolishedLag_RK5_CK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyRootPolishedLag_RK4_CK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyRootPolishedLag_RK3_CK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyRootPolishedLag_RK2_CK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyRootPolishedLag_RK1_CK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPolyRootPolishedLag_RK5_RK5
        use pm_kind, only: TKG => RK5
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPolyRootPolishedLag_RK4_RK4
        use pm_kind, only: TKG => RK4
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPolyRootPolishedLag_RK3_RK3
        use pm_kind, only: TKG => RK3
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPolyRootPolishedLag_RK2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPolyRootPolishedLag_RK1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_polynomial@routines.inc.F90"
    end procedure
#endif

#undef RK_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Lag_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPolyRootPolished_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines