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
!>  This file contains procedure implementations of [pm_matrixDet](@ref pm_matrixDet).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_matrixDet) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_matrixChol, only: setMatChol, iteration
    use pm_matrixCopy, only: setMatCopy, rdpack
    use pm_matrixLUP, only: setMatLUP
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatDet_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatDet_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatDet_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatDet_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatDet_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatDet_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatDet_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatDet_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatDet_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatDet_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatDet_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatDet_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatDet_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatDet_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatDet_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatDet_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatDet_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatDet_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatDet_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatDet_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatDet_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatDet_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatDet_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatDet_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatDetSqrt_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatDetSqrt_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatDetSqrt_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatDetSqrt_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatDetSqrt_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatDetSqrt_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatDetSqrt_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatDetSqrt_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatDetSqrt_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatDetSqrt_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatDetSqrt_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatDetSqrt_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatDetSqrt_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatDetSqrt_UXD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatDetSqrt_UXD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatDetSqrt_UXD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatDetSqrt_UXD_ONO_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatDetSqrt_UXD_ONO_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatDetSqrt_UXD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatDetSqrt_UXD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatDetSqrt_UXD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatDetSqrt_UXD_ONO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatDetSqrt_UXD_ONO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
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
    module procedure setMatDetSqrt_XLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatDetSqrt_XLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatDetSqrt_XLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatDetSqrt_XLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatDetSqrt_XLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatDetSqrt_XLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatDetSqrt_XLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatDetSqrt_XLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatDetSqrt_XLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatDetSqrt_XLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTH_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatDetSqrt_UXD_OTH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatDetSqrt_UXD_OTH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatDetSqrt_UXD_OTH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatDetSqrt_UXD_OTH_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatDetSqrt_UXD_OTH_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatDetSqrt_UXD_OTH_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatDetSqrt_UXD_OTH_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatDetSqrt_UXD_OTH_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatDetSqrt_UXD_OTH_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatDetSqrt_UXD_OTH_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
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
    module procedure setMatDetSqrt_XLD_OTH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatDetSqrt_XLD_OTH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatDetSqrt_XLD_OTH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatDetSqrt_XLD_OTH_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatDetSqrt_XLD_OTH_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatDetSqrt_XLD_OTH_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatDetSqrt_XLD_OTH_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatDetSqrt_XLD_OTH_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatDetSqrt_XLD_OTH_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatDetSqrt_XLD_OTH_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatDetSqrt_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatDetSqrtLog_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getMatDetSqrtLog_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getMatDetSqrtLog_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getMatDetSqrtLog_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getMatDetSqrtLog_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getMatDetSqrtLog_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMatDetSqrtLog_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMatDetSqrtLog_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMatDetSqrtLog_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMatDetSqrtLog_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMatDetSqrtLog_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatDetSqrtLog_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatDetSqrtLog_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatDetSqrtLog_UXD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatDetSqrtLog_UXD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatDetSqrtLog_UXD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatDetSqrtLog_UXD_ONO_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatDetSqrtLog_UXD_ONO_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatDetSqrtLog_UXD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatDetSqrtLog_UXD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatDetSqrtLog_UXD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatDetSqrtLog_UXD_ONO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatDetSqrtLog_UXD_ONO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
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
    module procedure setMatDetSqrtLog_XLD_ONO_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatDetSqrtLog_XLD_ONO_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatDetSqrtLog_XLD_ONO_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatDetSqrtLog_XLD_ONO_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatDetSqrtLog_XLD_ONO_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatDetSqrtLog_XLD_ONO_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatDetSqrtLog_XLD_ONO_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatDetSqrtLog_XLD_ONO_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatDetSqrtLog_XLD_ONO_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatDetSqrtLog_XLD_ONO_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTH_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setMatDetSqrtLog_UXD_OTH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatDetSqrtLog_UXD_OTH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatDetSqrtLog_UXD_OTH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatDetSqrtLog_UXD_OTH_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatDetSqrtLog_UXD_OTH_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatDetSqrtLog_UXD_OTH_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatDetSqrtLog_UXD_OTH_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatDetSqrtLog_UXD_OTH_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatDetSqrtLog_UXD_OTH_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatDetSqrtLog_UXD_OTH_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
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
    module procedure setMatDetSqrtLog_XLD_OTH_CK5
        use pm_kind, only: TKG => CK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setMatDetSqrtLog_XLD_OTH_CK4
        use pm_kind, only: TKG => CK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setMatDetSqrtLog_XLD_OTH_CK3
        use pm_kind, only: TKG => CK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setMatDetSqrtLog_XLD_OTH_CK2
        use pm_kind, only: TKG => CK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setMatDetSqrtLog_XLD_OTH_CK1
        use pm_kind, only: TKG => CK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMatDetSqrtLog_XLD_OTH_RK5
        use pm_kind, only: TKG => RK5
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMatDetSqrtLog_XLD_OTH_RK4
        use pm_kind, only: TKG => RK4
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMatDetSqrtLog_XLD_OTH_RK3
        use pm_kind, only: TKG => RK3
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMatDetSqrtLog_XLD_OTH_RK2
        use pm_kind, only: TKG => RK2
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMatDetSqrtLog_XLD_OTH_RK1
        use pm_kind, only: TKG => RK1
#include "pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatDetSqrtLog_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines