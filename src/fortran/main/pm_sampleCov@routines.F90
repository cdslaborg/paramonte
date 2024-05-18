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
!>  This file contains procedure implementations of [pm_sampleCov](@ref pm_sampleCov).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 2:48 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sampleCov) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_sampleMean, only: getMean, setMean
    use pm_complexMinMax, only: minval, maxval
    use pm_complexCompareAll, only: operator(<=)
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCov_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CorStd_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCovCorStd_ULD_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovCorStd_ULD_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovCorStd_ULD_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovCorStd_ULD_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovCorStd_ULD_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovCorStd_ULD_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovCorStd_ULD_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovCorStd_ULD_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovCorStd_ULD_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovCorStd_ULD_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCovCorStd_ULD_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovCorStd_ULD_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovCorStd_ULD_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovCorStd_ULD_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovCorStd_ULD_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovCorStd_ULD_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovCorStd_ULD_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovCorStd_ULD_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovCorStd_ULD_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovCorStd_ULD_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CorStd_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCov_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCov_ENABLED 1

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
    module procedure getCovWNO_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovWNO_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovWNO_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovWNO_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovWNO_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovWNO_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovWNO_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovWNO_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovWNO_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovWNO_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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
    module procedure getCovWTI_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovWTI_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovWTI_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovWTI_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovWTI_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovWTI_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovWTI_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovWTI_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovWTI_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovWTI_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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
    module procedure getCovWTR_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovWTR_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovWTR_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovWTR_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovWTR_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovWTR_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovWTR_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovWTR_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovWTR_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovWTR_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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

#undef getCov_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCov_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCovWNO_ULD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovWNO_ULD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovWNO_ULD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovWNO_ULD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovWNO_ULD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovWNO_ULD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovWNO_ULD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovWNO_ULD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovWNO_ULD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovWNO_ULD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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
    module procedure getCovWTI_ULD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovWTI_ULD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovWTI_ULD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovWTI_ULD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovWTI_ULD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovWTI_ULD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovWTI_ULD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovWTI_ULD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovWTI_ULD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovWTI_ULD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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
    module procedure getCovWTR_ULD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovWTR_ULD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovWTR_ULD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovWTR_ULD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovWTR_ULD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovWTR_ULD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovWTR_ULD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovWTR_ULD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovWTR_ULD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovWTR_ULD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCov_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCov_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CorStd_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovCorStd_UXD_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovCorStd_UXD_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovCorStd_UXD_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovCorStd_UXD_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovCorStd_UXD_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovCorStd_UXD_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovCorStd_UXD_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovCorStd_UXD_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovCorStd_UXD_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovCorStd_UXD_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovCorStd_UXD_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovCorStd_UXD_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovCorStd_UXD_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovCorStd_UXD_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovCorStd_UXD_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovCorStd_UXD_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovCorStd_UXD_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovCorStd_UXD_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovCorStd_UXD_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovCorStd_UXD_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovCorStd_XLD_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovCorStd_XLD_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovCorStd_XLD_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovCorStd_XLD_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovCorStd_XLD_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovCorStd_XLD_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovCorStd_XLD_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovCorStd_XLD_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovCorStd_XLD_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovCorStd_XLD_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovCorStd_XLD_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovCorStd_XLD_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovCorStd_XLD_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovCorStd_XLD_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovCorStd_XLD_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovCorStd_XLD_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovCorStd_XLD_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovCorStd_XLD_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovCorStd_XLD_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovCorStd_XLD_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CorStd_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCov_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCov_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XY_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovAvgWNO_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovAvgWNO_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovAvgWNO_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovAvgWNO_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovAvgWNO_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovAvgWNO_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovAvgWNO_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovAvgWNO_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovAvgWNO_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovAvgWNO_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovOrgWNO_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovOrgWNO_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovOrgWNO_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovOrgWNO_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovOrgWNO_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovOrgWNO_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovOrgWNO_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovOrgWNO_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovOrgWNO_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovOrgWNO_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

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

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovAvgWTI_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovAvgWTI_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovAvgWTI_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovAvgWTI_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovAvgWTI_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovAvgWTI_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovAvgWTI_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovAvgWTI_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovAvgWTI_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovAvgWTI_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovOrgWTI_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovOrgWTI_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovOrgWTI_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovOrgWTI_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovOrgWTI_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovOrgWTI_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovOrgWTI_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovOrgWTI_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovOrgWTI_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovOrgWTI_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

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

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovAvgWTR_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovAvgWTR_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovAvgWTR_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovAvgWTR_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovAvgWTR_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovAvgWTR_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovAvgWTR_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovAvgWTR_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovAvgWTR_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovAvgWTR_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovOrgWTR_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovOrgWTR_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovOrgWTR_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovOrgWTR_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovOrgWTR_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovOrgWTR_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovOrgWTR_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovOrgWTR_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovOrgWTR_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovOrgWTR_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XY_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCov_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCov_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovAvgWNO_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovAvgWNO_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovAvgWNO_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovAvgWNO_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovAvgWNO_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovAvgWNO_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovAvgWNO_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovAvgWNO_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovAvgWNO_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovAvgWNO_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovOrgWNO_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovOrgWNO_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovOrgWNO_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovOrgWNO_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovOrgWNO_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovOrgWNO_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovOrgWNO_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovOrgWNO_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovOrgWNO_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovOrgWNO_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

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

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovAvgWTI_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovAvgWTI_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovAvgWTI_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovAvgWTI_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovAvgWTI_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovAvgWTI_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovAvgWTI_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovAvgWTI_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovAvgWTI_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovAvgWTI_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovOrgWTI_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovOrgWTI_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovOrgWTI_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovOrgWTI_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovOrgWTI_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovOrgWTI_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovOrgWTI_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovOrgWTI_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovOrgWTI_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovOrgWTI_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

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

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovAvgWTR_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovAvgWTR_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovAvgWTR_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovAvgWTR_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovAvgWTR_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovAvgWTR_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovAvgWTR_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovAvgWTR_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovAvgWTR_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovAvgWTR_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovOrgWTR_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovOrgWTR_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovOrgWTR_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovOrgWTR_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovOrgWTR_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovOrgWTR_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovOrgWTR_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovOrgWTR_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovOrgWTR_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovOrgWTR_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCov_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCov_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovAvgWNO_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovAvgWNO_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovAvgWNO_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovAvgWNO_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovAvgWNO_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovAvgWNO_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovAvgWNO_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovAvgWNO_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovAvgWNO_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovAvgWNO_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovOrgWNO_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovOrgWNO_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovOrgWNO_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovOrgWNO_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovOrgWNO_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovOrgWNO_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovOrgWNO_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovOrgWNO_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovOrgWNO_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovOrgWNO_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

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

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovAvgWTI_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovAvgWTI_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovAvgWTI_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovAvgWTI_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovAvgWTI_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovAvgWTI_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovAvgWTI_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovAvgWTI_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovAvgWTI_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovAvgWTI_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovOrgWTI_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovOrgWTI_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovOrgWTI_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovOrgWTI_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovOrgWTI_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovOrgWTI_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovOrgWTI_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovOrgWTI_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovOrgWTI_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovOrgWTI_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

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

#define Avg_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovAvgWTR_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovAvgWTR_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovAvgWTR_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovAvgWTR_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovAvgWTR_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovAvgWTR_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovAvgWTR_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovAvgWTR_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovAvgWTR_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovAvgWTR_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Avg_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Org_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovOrgWTR_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovOrgWTR_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovOrgWTR_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovOrgWTR_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovOrgWTR_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovOrgWTR_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovOrgWTR_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovOrgWTR_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovOrgWTR_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovOrgWTR_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Org_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCov_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovMean_ENABLED 1

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
    module procedure setCovMeanWNO_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanWNO_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanWNO_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanWNO_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanWNO_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanWNO_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanWNO_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanWNO_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanWNO_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanWNO_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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
    module procedure setCovMeanWTI_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanWTI_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanWTI_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanWTI_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanWTI_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanWTI_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanWTI_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanWTI_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanWTI_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanWTI_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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
    module procedure setCovMeanWTR_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanWTR_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanWTR_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanWTR_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanWTR_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanWTR_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanWTR_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanWTR_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanWTR_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanWTR_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XY_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMeanWNO_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanWNO_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanWNO_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanWNO_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanWNO_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanWNO_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanWNO_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanWNO_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanWNO_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanWNO_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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
    module procedure setCovMeanWTI_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanWTI_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanWTI_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanWTI_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanWTI_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanWTI_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanWTI_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanWTI_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanWTI_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanWTI_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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
    module procedure setCovMeanWTR_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanWTR_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanWTR_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanWTR_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanWTR_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanWTR_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanWTR_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanWTR_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanWTR_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanWTR_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMeanWNO_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanWNO_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanWNO_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanWNO_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanWNO_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanWNO_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanWNO_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanWNO_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanWNO_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanWNO_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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
    module procedure setCovMeanWTI_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanWTI_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanWTI_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanWTI_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanWTI_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanWTI_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanWTI_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanWTI_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanWTI_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanWTI_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
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
    module procedure setCovMeanWTR_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanWTR_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanWTR_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanWTR_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanWTR_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanWTR_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanWTR_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanWTR_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanWTR_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanWTR_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCovMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCovMergedNew_RDP_ULD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCovMergedNew_RDP_ULD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCovMergedNew_RDP_ULD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCovMergedNew_RDP_ULD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCovMergedNew_RDP_ULD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCovMergedNew_RDP_ULD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCovMergedNew_RDP_ULD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCovMergedNew_RDP_ULD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCovMergedNew_RDP_ULD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCovMergedNew_RDP_ULD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCovMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMergedNew_RDP_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMergedNew_RDP_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMergedNew_RDP_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMergedNew_RDP_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMergedNew_RDP_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMergedNew_RDP_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMergedNew_RDP_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMergedNew_RDP_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMergedNew_RDP_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMergedNew_RDP_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMergedNew_RDP_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMergedNew_RDP_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMergedNew_RDP_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMergedNew_RDP_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMergedNew_RDP_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMergedNew_RDP_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMergedNew_RDP_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMergedNew_RDP_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMergedNew_RDP_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMergedNew_RDP_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

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

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMergedOld_RDP_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMergedOld_RDP_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMergedOld_RDP_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMergedOld_RDP_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMergedOld_RDP_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMergedOld_RDP_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMergedOld_RDP_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMergedOld_RDP_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMergedOld_RDP_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMergedOld_RDP_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMergedOld_RDP_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMergedOld_RDP_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMergedOld_RDP_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMergedOld_RDP_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMergedOld_RDP_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMergedOld_RDP_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMergedOld_RDP_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMergedOld_RDP_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMergedOld_RDP_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMergedOld_RDP_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovMeanMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMeanMergedNew_RDP_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanMergedNew_RDP_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanMergedNew_RDP_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanMergedNew_RDP_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanMergedNew_RDP_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanMergedNew_RDP_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanMergedNew_RDP_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanMergedNew_RDP_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanMergedNew_RDP_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanMergedNew_RDP_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMeanMergedNew_RDP_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanMergedNew_RDP_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanMergedNew_RDP_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanMergedNew_RDP_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanMergedNew_RDP_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanMergedNew_RDP_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanMergedNew_RDP_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanMergedNew_RDP_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanMergedNew_RDP_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanMergedNew_RDP_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

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

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMeanMergedOld_RDP_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanMergedOld_RDP_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanMergedOld_RDP_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanMergedOld_RDP_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanMergedOld_RDP_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanMergedOld_RDP_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanMergedOld_RDP_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanMergedOld_RDP_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanMergedOld_RDP_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanMergedOld_RDP_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMeanMergedOld_RDP_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanMergedOld_RDP_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanMergedOld_RDP_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanMergedOld_RDP_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanMergedOld_RDP_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanMergedOld_RDP_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanMergedOld_RDP_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanMergedOld_RDP_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanMergedOld_RDP_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanMergedOld_RDP_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovMeanMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovUpdated_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovUpdatedOld_RDP_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovUpdatedOld_RDP_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovUpdatedOld_RDP_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovUpdatedOld_RDP_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovUpdatedOld_RDP_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovUpdatedOld_RDP_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovUpdatedOld_RDP_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovUpdatedOld_RDP_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovUpdatedOld_RDP_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovUpdatedOld_RDP_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovUpdatedOld_RDP_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovUpdatedOld_RDP_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovUpdatedOld_RDP_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovUpdatedOld_RDP_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovUpdatedOld_RDP_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovUpdatedOld_RDP_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovUpdatedOld_RDP_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovUpdatedOld_RDP_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovUpdatedOld_RDP_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovUpdatedOld_RDP_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovUpdated_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovMeanUpdated_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Old_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCovMeanUpdatedOld_RDP_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Old_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovMeanUpdated_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines