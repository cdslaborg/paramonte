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
!>  This file contains procedure implementations of [pm_sampleMean](@ref pm_sampleMean).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sampleMean) routines ! LCOV_EXCL_LINE

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

#define getMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanALL_WNO_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanALL_WNO_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanALL_WNO_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanALL_WNO_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanALL_WNO_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanALL_WNO_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanALL_WNO_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanALL_WNO_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanALL_WNO_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanALL_WNO_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanALL_WNO_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanALL_WNO_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanALL_WNO_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanALL_WNO_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanALL_WNO_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanALL_WNO_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanALL_WNO_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanALL_WNO_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanALL_WNO_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanALL_WNO_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanALL_WTI_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanALL_WTI_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanALL_WTI_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanALL_WTI_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanALL_WTI_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanALL_WTI_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanALL_WTI_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanALL_WTI_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanALL_WTI_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanALL_WTI_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanALL_WTI_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanALL_WTI_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanALL_WTI_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanALL_WTI_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanALL_WTI_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanALL_WTI_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanALL_WTI_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanALL_WTI_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanALL_WTI_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanALL_WTI_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanALL_WTR_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanALL_WTR_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanALL_WTR_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanALL_WTR_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanALL_WTR_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanALL_WTR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanALL_WTR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanALL_WTR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanALL_WTR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanALL_WTR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanALL_WTR_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanALL_WTR_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanALL_WTR_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanALL_WTR_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanALL_WTR_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanALL_WTR_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanALL_WTR_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanALL_WTR_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanALL_WTR_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanALL_WTR_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanDIM_WNO_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanDIM_WNO_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanDIM_WNO_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanDIM_WNO_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanDIM_WNO_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanDIM_WNO_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanDIM_WNO_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanDIM_WNO_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanDIM_WNO_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanDIM_WNO_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanDIM_WNO_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanDIM_WNO_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanDIM_WNO_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanDIM_WNO_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanDIM_WNO_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanDIM_WNO_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanDIM_WNO_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanDIM_WNO_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanDIM_WNO_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanDIM_WNO_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanDIM_WTI_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanDIM_WTI_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanDIM_WTI_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanDIM_WTI_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanDIM_WTI_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanDIM_WTI_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanDIM_WTI_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanDIM_WTI_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanDIM_WTI_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanDIM_WTI_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanDIM_WTI_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanDIM_WTI_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanDIM_WTI_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanDIM_WTI_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanDIM_WTI_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanDIM_WTI_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanDIM_WTI_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanDIM_WTI_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanDIM_WTI_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanDIM_WTI_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanDIM_WTR_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanDIM_WTR_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanDIM_WTR_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanDIM_WTR_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanDIM_WTR_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanDIM_WTR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanDIM_WTR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanDIM_WTR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanDIM_WTR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanDIM_WTR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanDIM_WTR_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanDIM_WTR_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanDIM_WTR_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanDIM_WTR_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanDIM_WTR_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanDIM_WTR_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanDIM_WTR_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanDIM_WTR_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanDIM_WTR_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanDIM_WTR_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DIM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XY_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanALL_WNO_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanALL_WNO_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanALL_WNO_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanALL_WNO_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanALL_WNO_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanALL_WNO_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanALL_WNO_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanALL_WNO_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanALL_WNO_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanALL_WNO_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanALL_WTI_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanALL_WTI_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanALL_WTI_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanALL_WTI_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanALL_WTI_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanALL_WTI_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanALL_WTI_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanALL_WTI_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanALL_WTI_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanALL_WTI_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanALL_WTR_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanALL_WTR_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanALL_WTR_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanALL_WTR_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanALL_WTR_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanALL_WTR_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanALL_WTR_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanALL_WTR_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanALL_WTR_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanALL_WTR_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XY_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanALL_WNO_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanALL_WNO_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanALL_WNO_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanALL_WNO_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanALL_WNO_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanALL_WNO_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanALL_WNO_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanALL_WNO_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanALL_WNO_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanALL_WNO_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanALL_WNO_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanALL_WNO_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanALL_WNO_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanALL_WNO_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanALL_WNO_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanALL_WNO_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanALL_WNO_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanALL_WNO_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanALL_WNO_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanALL_WNO_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanALL_WTI_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanALL_WTI_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanALL_WTI_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanALL_WTI_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanALL_WTI_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanALL_WTI_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanALL_WTI_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanALL_WTI_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanALL_WTI_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanALL_WTI_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanALL_WTI_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanALL_WTI_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanALL_WTI_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanALL_WTI_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanALL_WTI_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanALL_WTI_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanALL_WTI_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanALL_WTI_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanALL_WTI_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanALL_WTI_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanALL_WTR_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanALL_WTR_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanALL_WTR_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanALL_WTR_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanALL_WTR_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanALL_WTR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanALL_WTR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanALL_WTR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanALL_WTR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanALL_WTR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanALL_WTR_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanALL_WTR_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanALL_WTR_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanALL_WTR_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanALL_WTR_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanALL_WTR_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanALL_WTR_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanALL_WTR_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanALL_WTR_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanALL_WTR_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanDIM_WNO_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanDIM_WNO_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanDIM_WNO_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanDIM_WNO_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanDIM_WNO_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanDIM_WNO_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanDIM_WNO_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanDIM_WNO_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanDIM_WNO_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanDIM_WNO_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanDIM_WNO_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanDIM_WNO_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanDIM_WNO_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanDIM_WNO_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanDIM_WNO_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanDIM_WNO_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanDIM_WNO_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanDIM_WNO_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanDIM_WNO_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanDIM_WNO_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanDIM_WTI_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanDIM_WTI_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanDIM_WTI_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanDIM_WTI_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanDIM_WTI_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanDIM_WTI_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanDIM_WTI_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanDIM_WTI_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanDIM_WTI_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanDIM_WTI_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanDIM_WTI_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanDIM_WTI_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanDIM_WTI_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanDIM_WTI_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanDIM_WTI_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanDIM_WTI_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanDIM_WTI_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanDIM_WTI_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanDIM_WTI_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanDIM_WTI_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanDIM_WTR_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanDIM_WTR_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanDIM_WTR_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanDIM_WTR_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanDIM_WTR_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanDIM_WTR_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanDIM_WTR_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanDIM_WTR_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanDIM_WTR_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanDIM_WTR_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
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

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanDIM_WTR_D2_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanDIM_WTR_D2_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanDIM_WTR_D2_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanDIM_WTR_D2_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanDIM_WTR_D2_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanDIM_WTR_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanDIM_WTR_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanDIM_WTR_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanDIM_WTR_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanDIM_WTR_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DIM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMeanMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanMergedNew_D0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanMergedNew_D0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanMergedNew_D0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanMergedNew_D0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanMergedNew_D0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanMergedNew_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanMergedNew_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanMergedNew_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanMergedNew_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanMergedNew_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMeanMergedNew_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMeanMergedNew_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMeanMergedNew_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMeanMergedNew_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMeanMergedNew_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMeanMergedNew_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMeanMergedNew_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMeanMergedNew_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMeanMergedNew_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMeanMergedNew_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMeanMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMeanMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanMergedNew_D0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanMergedNew_D0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanMergedNew_D0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanMergedNew_D0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanMergedNew_D0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanMergedNew_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanMergedNew_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanMergedNew_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanMergedNew_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanMergedNew_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanMergedNew_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanMergedNew_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanMergedNew_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanMergedNew_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanMergedNew_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanMergedNew_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanMergedNew_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanMergedNew_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanMergedNew_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanMergedNew_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanMergedOld_D0_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanMergedOld_D0_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanMergedOld_D0_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanMergedOld_D0_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanMergedOld_D0_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanMergedOld_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanMergedOld_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanMergedOld_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanMergedOld_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanMergedOld_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMeanMergedOld_D1_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMeanMergedOld_D1_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMeanMergedOld_D1_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMeanMergedOld_D1_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMeanMergedOld_D1_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMeanMergedOld_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMeanMergedOld_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMeanMergedOld_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMeanMergedOld_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMeanMergedOld_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleMean@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMeanMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines