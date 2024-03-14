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
!>  This file contains procedure implementations of [pm_sampleAffinity](@ref pm_sampleAffinity).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sampleAffinity) routines ! LCOV_EXCL_LINE

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

#define getAffinity_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getAffinity_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getAffinity_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getAffinity_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getAffinity_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getAffinity_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getAffinity_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getAffinity_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getAffinity_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getAffinity_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getAffinity_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getAffinity_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setAffinity_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ATL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CGR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CGR_ATL_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CGR_ATL_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CGR_ATL_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CGR_ATL_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CGR_ATL_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CGR_ATL_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CGR_ATL_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CGR_ATL_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CGR_ATL_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CGR_ATL_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CUD_ATL_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CUD_ATL_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CUD_ATL_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CUD_ATL_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CUD_ATL_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CUD_ATL_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CUD_ATL_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CUD_ATL_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CUD_ATL_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CUD_ATL_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CLD_ATL_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CLD_ATL_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CLD_ATL_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CLD_ATL_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CLD_ATL_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CLD_ATL_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CLD_ATL_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CLD_ATL_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CLD_ATL_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CLD_ATL_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CUU_ATL_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CUU_ATL_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CUU_ATL_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CUU_ATL_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CUU_ATL_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CUU_ATL_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CUU_ATL_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CUU_ATL_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CUU_ATL_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CUU_ATL_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CLU_ATL_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CLU_ATL_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CLU_ATL_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CLU_ATL_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CLU_ATL_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CLU_ATL_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CLU_ATL_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CLU_ATL_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CLU_ATL_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CLU_ATL_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ATL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DTL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CGR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CGR_DTL_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CGR_DTL_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CGR_DTL_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CGR_DTL_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CGR_DTL_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CGR_DTL_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CGR_DTL_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CGR_DTL_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CGR_DTL_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CGR_DTL_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CGR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CUD_DTL_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CUD_DTL_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CUD_DTL_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CUD_DTL_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CUD_DTL_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CUD_DTL_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CUD_DTL_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CUD_DTL_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CUD_DTL_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CUD_DTL_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CLD_DTL_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CLD_DTL_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CLD_DTL_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CLD_DTL_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CLD_DTL_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CLD_DTL_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CLD_DTL_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CLD_DTL_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CLD_DTL_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CLD_DTL_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CUU_DTL_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CUU_DTL_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CUU_DTL_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CUU_DTL_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CUU_DTL_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CUU_DTL_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CUU_DTL_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CUU_DTL_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CUU_DTL_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CUU_DTL_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CUU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CLU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CLU_DTL_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CLU_DTL_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CLU_DTL_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CLU_DTL_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CLU_DTL_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setAffinity_CLU_DTL_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setAffinity_CLU_DTL_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setAffinity_CLU_DTL_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setAffinity_CLU_DTL_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setAffinity_CLU_DTL_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_sampleAffinity@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CLU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DTL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setAffinity_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines