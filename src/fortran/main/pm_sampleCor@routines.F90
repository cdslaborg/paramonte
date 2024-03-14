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
!>  This file contains procedure implementations of [pm_sampleCor](@ref pm_sampleCor).
!>
!>  final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sampleCor) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_arrayRank, only: setRankFractional
    use pm_sampleMean, only: setMean
    use pm_sampleCov, only: setCov
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCor_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CFC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VUXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCFC_RULD_VUXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCFC_RULD_VUXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCFC_RULD_VUXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCFC_RULD_VUXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCFC_RULD_VUXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCFC_RULD_VUXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCFC_RULD_VUXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCFC_RULD_VUXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCFC_RULD_VUXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCFC_RULD_VUXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VUXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VXLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCFC_RULD_VXLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCFC_RULD_VXLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCFC_RULD_VXLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCFC_RULD_VXLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCFC_RULD_VXLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCFC_RULD_VXLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCFC_RULD_VXLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCFC_RULD_VXLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCFC_RULD_VXLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCFC_RULD_VXLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VXLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VUXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCFC_RULD_VUXX_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCFC_RULD_VUXX_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCFC_RULD_VUXX_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCFC_RULD_VUXX_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCFC_RULD_VUXX_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCFC_RULD_VUXX_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCFC_RULD_VUXX_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCFC_RULD_VUXX_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCFC_RULD_VUXX_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCFC_RULD_VUXX_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VUXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VXLX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getCFC_RULD_VXLX_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getCFC_RULD_VXLX_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getCFC_RULD_VXLX_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getCFC_RULD_VXLX_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getCFC_RULD_VXLX_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCFC_RULD_VXLX_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCFC_RULD_VXLX_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCFC_RULD_VXLX_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCFC_RULD_VXLX_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCFC_RULD_VXLX_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VXLX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CFC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCor_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCor_ENABLED 1

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
    module procedure getPrsWNO_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPrsWNO_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPrsWNO_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPrsWNO_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPrsWNO_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPrsWNO_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPrsWNO_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPrsWNO_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPrsWNO_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPrsWNO_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure getPrsWTI_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPrsWTI_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPrsWTI_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPrsWTI_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPrsWTI_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPrsWTI_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPrsWTI_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPrsWTI_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPrsWTI_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPrsWTI_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure getPrsWTR_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPrsWTR_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPrsWTR_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPrsWTR_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPrsWTR_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPrsWTR_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPrsWTR_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPrsWTR_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPrsWTR_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPrsWTR_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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

#undef getCor_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCor_ENABLED 1

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
    module procedure getPrsWNO_ULD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPrsWNO_ULD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPrsWNO_ULD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPrsWNO_ULD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPrsWNO_ULD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPrsWNO_ULD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPrsWNO_ULD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPrsWNO_ULD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPrsWNO_ULD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPrsWNO_ULD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure getPrsWTI_ULD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPrsWTI_ULD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPrsWTI_ULD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPrsWTI_ULD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPrsWTI_ULD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPrsWTI_ULD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPrsWTI_ULD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPrsWTI_ULD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPrsWTI_ULD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPrsWTI_ULD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure getPrsWTR_ULD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getPrsWTR_ULD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getPrsWTR_ULD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getPrsWTR_ULD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getPrsWTR_ULD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPrsWTR_ULD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPrsWTR_ULD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPrsWTR_ULD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPrsWTR_ULD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPrsWTR_ULD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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

#undef getCor_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCor_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CFC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RUXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VUXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RUXX_VUXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RUXX_VUXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RUXX_VUXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RUXX_VUXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RUXX_VUXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RUXX_VUXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RUXX_VUXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RUXX_VUXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RUXX_VUXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RUXX_VUXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VUXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VUXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RUXX_VUXX_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RUXX_VUXX_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RUXX_VUXX_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RUXX_VUXX_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RUXX_VUXX_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RUXX_VUXX_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RUXX_VUXX_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RUXX_VUXX_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RUXX_VUXX_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RUXX_VUXX_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VUXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VXLX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RUXX_VXLX_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RUXX_VXLX_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RUXX_VXLX_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RUXX_VXLX_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RUXX_VXLX_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RUXX_VXLX_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RUXX_VXLX_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RUXX_VXLX_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RUXX_VXLX_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RUXX_VXLX_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VXLX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VXLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RUXX_VXLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RUXX_VXLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RUXX_VXLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RUXX_VXLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RUXX_VXLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RUXX_VXLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RUXX_VXLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RUXX_VXLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RUXX_VXLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RUXX_VXLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VXLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RUXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RUXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VUXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RUXD_VUXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RUXD_VUXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RUXD_VUXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RUXD_VUXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RUXD_VUXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RUXD_VUXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RUXD_VUXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RUXD_VUXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RUXD_VUXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RUXD_VUXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VUXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VUXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RUXD_VUXX_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RUXD_VUXX_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RUXD_VUXX_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RUXD_VUXX_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RUXD_VUXX_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RUXD_VUXX_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RUXD_VUXX_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RUXD_VUXX_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RUXD_VUXX_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RUXD_VUXX_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VUXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VXLX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RUXD_VXLX_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RUXD_VXLX_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RUXD_VXLX_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RUXD_VXLX_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RUXD_VXLX_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RUXD_VXLX_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RUXD_VXLX_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RUXD_VXLX_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RUXD_VXLX_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RUXD_VXLX_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VXLX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VXLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RUXD_VXLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RUXD_VXLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RUXD_VXLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RUXD_VXLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RUXD_VXLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RUXD_VXLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RUXD_VXLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RUXD_VXLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RUXD_VXLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RUXD_VXLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VXLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RUXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RXLX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VUXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RXLX_VUXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RXLX_VUXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RXLX_VUXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RXLX_VUXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RXLX_VUXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RXLX_VUXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RXLX_VUXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RXLX_VUXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RXLX_VUXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RXLX_VUXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VUXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VUXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RXLX_VUXX_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RXLX_VUXX_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RXLX_VUXX_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RXLX_VUXX_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RXLX_VUXX_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RXLX_VUXX_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RXLX_VUXX_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RXLX_VUXX_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RXLX_VUXX_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RXLX_VUXX_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VUXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VXLX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RXLX_VXLX_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RXLX_VXLX_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RXLX_VXLX_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RXLX_VXLX_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RXLX_VXLX_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RXLX_VXLX_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RXLX_VXLX_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RXLX_VXLX_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RXLX_VXLX_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RXLX_VXLX_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VXLX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VXLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RXLX_VXLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RXLX_VXLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RXLX_VXLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RXLX_VXLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RXLX_VXLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RXLX_VXLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RXLX_VXLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RXLX_VXLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RXLX_VXLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RXLX_VXLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VXLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RXLX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RXLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VUXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RXLD_VUXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RXLD_VUXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RXLD_VUXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RXLD_VUXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RXLD_VUXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RXLD_VUXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RXLD_VUXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RXLD_VUXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RXLD_VUXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RXLD_VUXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VUXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VUXX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RXLD_VUXX_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RXLD_VUXX_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RXLD_VUXX_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RXLD_VUXX_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RXLD_VUXX_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RXLD_VUXX_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RXLD_VUXX_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RXLD_VUXX_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RXLD_VUXX_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RXLD_VUXX_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VUXX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VXLX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RXLD_VXLX_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RXLD_VXLX_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RXLD_VXLX_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RXLD_VXLX_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RXLD_VXLX_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RXLD_VXLX_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RXLD_VXLX_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RXLD_VXLX_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RXLD_VXLX_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RXLD_VXLX_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VXLX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define VXLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setCFC_RXLD_VXLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setCFC_RXLD_VXLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setCFC_RXLD_VXLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setCFC_RXLD_VXLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setCFC_RXLD_VXLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCFC_RXLD_VXLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCFC_RXLD_VXLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCFC_RXLD_VXLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCFC_RXLD_VXLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCFC_RXLD_VXLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef VXLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RXLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CFC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCor_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCor_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Prs_ENABLED 1

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
    module procedure setPrsAvgWNO_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsAvgWNO_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsAvgWNO_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsAvgWNO_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsAvgWNO_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsAvgWNO_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsAvgWNO_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsAvgWNO_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsAvgWNO_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsAvgWNO_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsOrgWNO_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsOrgWNO_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsOrgWNO_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsOrgWNO_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsOrgWNO_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsOrgWNO_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsOrgWNO_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsOrgWNO_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsOrgWNO_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsOrgWNO_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsAvgWTI_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsAvgWTI_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsAvgWTI_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsAvgWTI_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsAvgWTI_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsAvgWTI_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsAvgWTI_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsAvgWTI_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsAvgWTI_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsAvgWTI_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsOrgWTI_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsOrgWTI_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsOrgWTI_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsOrgWTI_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsOrgWTI_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsOrgWTI_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsOrgWTI_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsOrgWTI_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsOrgWTI_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsOrgWTI_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsAvgWTR_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsAvgWTR_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsAvgWTR_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsAvgWTR_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsAvgWTR_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsAvgWTR_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsAvgWTR_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsAvgWTR_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsAvgWTR_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsAvgWTR_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsOrgWTR_XY_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsOrgWTR_XY_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsOrgWTR_XY_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsOrgWTR_XY_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsOrgWTR_XY_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsOrgWTR_XY_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsOrgWTR_XY_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsOrgWTR_XY_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsOrgWTR_XY_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsOrgWTR_XY_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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

#undef XY_ENABLED

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
    module procedure setPrsAvgWNO_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsAvgWNO_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsAvgWNO_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsAvgWNO_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsAvgWNO_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsAvgWNO_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsAvgWNO_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsAvgWNO_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsAvgWNO_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsAvgWNO_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsOrgWNO_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsOrgWNO_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsOrgWNO_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsOrgWNO_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsOrgWNO_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsOrgWNO_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsOrgWNO_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsOrgWNO_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsOrgWNO_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsOrgWNO_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsAvgWTI_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsAvgWTI_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsAvgWTI_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsAvgWTI_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsAvgWTI_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsAvgWTI_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsAvgWTI_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsAvgWTI_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsAvgWTI_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsAvgWTI_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsOrgWTI_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsOrgWTI_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsOrgWTI_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsOrgWTI_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsOrgWTI_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsOrgWTI_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsOrgWTI_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsOrgWTI_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsOrgWTI_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsOrgWTI_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsAvgWTR_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsAvgWTR_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsAvgWTR_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsAvgWTR_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsAvgWTR_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsAvgWTR_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsAvgWTR_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsAvgWTR_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsAvgWTR_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsAvgWTR_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsOrgWTR_UXD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsOrgWTR_UXD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsOrgWTR_UXD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsOrgWTR_UXD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsOrgWTR_UXD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsOrgWTR_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsOrgWTR_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsOrgWTR_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsOrgWTR_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsOrgWTR_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsAvgWNO_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsAvgWNO_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsAvgWNO_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsAvgWNO_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsAvgWNO_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsAvgWNO_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsAvgWNO_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsAvgWNO_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsAvgWNO_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsAvgWNO_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsOrgWNO_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsOrgWNO_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsOrgWNO_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsOrgWNO_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsOrgWNO_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsOrgWNO_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsOrgWNO_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsOrgWNO_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsOrgWNO_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsOrgWNO_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsAvgWTI_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsAvgWTI_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsAvgWTI_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsAvgWTI_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsAvgWTI_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsAvgWTI_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsAvgWTI_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsAvgWTI_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsAvgWTI_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsAvgWTI_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsOrgWTI_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsOrgWTI_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsOrgWTI_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsOrgWTI_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsOrgWTI_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsOrgWTI_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsOrgWTI_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsOrgWTI_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsOrgWTI_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsOrgWTI_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsAvgWTR_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsAvgWTR_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsAvgWTR_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsAvgWTR_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsAvgWTR_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsAvgWTR_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsAvgWTR_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsAvgWTR_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsAvgWTR_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsAvgWTR_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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
    module procedure setPrsOrgWTR_XLD_CK5
        use pm_kind, only: TKC => CK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setPrsOrgWTR_XLD_CK4
        use pm_kind, only: TKC => CK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setPrsOrgWTR_XLD_CK3
        use pm_kind, only: TKC => CK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setPrsOrgWTR_XLD_CK2
        use pm_kind, only: TKC => CK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setPrsOrgWTR_XLD_CK1
        use pm_kind, only: TKC => CK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPrsOrgWTR_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPrsOrgWTR_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPrsOrgWTR_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPrsOrgWTR_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPrsOrgWTR_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
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

#undef Prs_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCor_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRho_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XY_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWNO_XY_D0_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWNO_XY_D0_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWNO_XY_D0_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWNO_XY_D0_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWNO_XY_D0_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWTI_XY_D0_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWTI_XY_D0_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWTI_XY_D0_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWTI_XY_D0_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWTI_XY_D0_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWTR_XY_D0_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWTR_XY_D0_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWTR_XY_D0_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWTR_XY_D0_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWTR_XY_D0_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XY_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRho_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRho_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XY_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWNO_XY_D1_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWNO_XY_D1_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWNO_XY_D1_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWNO_XY_D1_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWNO_XY_D1_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRhoWNO_XY_D1_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRhoWNO_XY_D1_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRhoWNO_XY_D1_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRhoWNO_XY_D1_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRhoWNO_XY_D1_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRhoWNO_XY_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRhoWNO_XY_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRhoWNO_XY_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRhoWNO_XY_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRhoWNO_XY_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWNO_XY_D1_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWNO_XY_D1_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWNO_XY_D1_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWNO_XY_D1_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWNO_XY_D1_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRhoWNO_XY_D1_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWTI_XY_D1_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWTI_XY_D1_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWTI_XY_D1_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWTI_XY_D1_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWTI_XY_D1_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRhoWTI_XY_D1_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRhoWTI_XY_D1_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRhoWTI_XY_D1_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRhoWTI_XY_D1_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRhoWTI_XY_D1_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRhoWTI_XY_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRhoWTI_XY_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRhoWTI_XY_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRhoWTI_XY_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRhoWTI_XY_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWTI_XY_D1_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWTI_XY_D1_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWTI_XY_D1_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWTI_XY_D1_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWTI_XY_D1_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRhoWTI_XY_D1_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWTR_XY_D1_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWTR_XY_D1_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWTR_XY_D1_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWTR_XY_D1_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWTR_XY_D1_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRhoWTR_XY_D1_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRhoWTR_XY_D1_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRhoWTR_XY_D1_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRhoWTR_XY_D1_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRhoWTR_XY_D1_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRhoWTR_XY_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRhoWTR_XY_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRhoWTR_XY_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRhoWTR_XY_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRhoWTR_XY_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWTR_XY_D1_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWTR_XY_D1_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWTR_XY_D1_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWTR_XY_D1_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWTR_XY_D1_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRhoWTR_XY_D1_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XY_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRho_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRho_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWNO_ULD_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWNO_ULD_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWNO_ULD_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWNO_ULD_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWNO_ULD_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRhoWNO_ULD_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRhoWNO_ULD_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRhoWNO_ULD_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRhoWNO_ULD_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRhoWNO_ULD_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRhoWNO_ULD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRhoWNO_ULD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRhoWNO_ULD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRhoWNO_ULD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRhoWNO_ULD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWNO_ULD_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWNO_ULD_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWNO_ULD_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWNO_ULD_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWNO_ULD_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRhoWNO_ULD_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWTI_ULD_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWTI_ULD_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWTI_ULD_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWTI_ULD_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWTI_ULD_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRhoWTI_ULD_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRhoWTI_ULD_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRhoWTI_ULD_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRhoWTI_ULD_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRhoWTI_ULD_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRhoWTI_ULD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRhoWTI_ULD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRhoWTI_ULD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRhoWTI_ULD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRhoWTI_ULD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWTI_ULD_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWTI_ULD_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWTI_ULD_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWTI_ULD_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWTI_ULD_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRhoWTI_ULD_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWTR_ULD_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWTR_ULD_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWTR_ULD_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWTR_ULD_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWTR_ULD_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getRhoWTR_ULD_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getRhoWTR_ULD_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getRhoWTR_ULD_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getRhoWTR_ULD_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getRhoWTR_ULD_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRhoWTR_ULD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRhoWTR_ULD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRhoWTR_ULD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRhoWTR_ULD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRhoWTR_ULD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure getRhoWTR_ULD_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getRhoWTR_ULD_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getRhoWTR_ULD_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getRhoWTR_ULD_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getRhoWTR_ULD_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure getRhoWTR_ULD_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRho_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRho_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XY_D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWNO_XY_D0_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWNO_XY_D0_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWNO_XY_D0_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWNO_XY_D0_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWNO_XY_D0_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTI_XY_D0_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTI_XY_D0_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTI_XY_D0_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTI_XY_D0_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTI_XY_D0_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTR_XY_D0_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTR_XY_D0_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTR_XY_D0_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTR_XY_D0_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTR_XY_D0_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XY_D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRho_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRho_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XY_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWNO_XY_D1_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWNO_XY_D1_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWNO_XY_D1_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWNO_XY_D1_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWNO_XY_D1_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRhoWNO_XY_D1_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRhoWNO_XY_D1_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRhoWNO_XY_D1_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRhoWNO_XY_D1_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRhoWNO_XY_D1_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRhoWNO_XY_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRhoWNO_XY_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRhoWNO_XY_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRhoWNO_XY_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRhoWNO_XY_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWNO_XY_D1_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWNO_XY_D1_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWNO_XY_D1_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWNO_XY_D1_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWNO_XY_D1_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRhoWNO_XY_D1_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTI_XY_D1_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTI_XY_D1_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTI_XY_D1_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTI_XY_D1_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTI_XY_D1_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRhoWTI_XY_D1_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRhoWTI_XY_D1_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRhoWTI_XY_D1_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRhoWTI_XY_D1_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRhoWTI_XY_D1_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRhoWTI_XY_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRhoWTI_XY_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRhoWTI_XY_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRhoWTI_XY_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRhoWTI_XY_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTI_XY_D1_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTI_XY_D1_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTI_XY_D1_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTI_XY_D1_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTI_XY_D1_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRhoWTI_XY_D1_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTR_XY_D1_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTR_XY_D1_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTR_XY_D1_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTR_XY_D1_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTR_XY_D1_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRhoWTR_XY_D1_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRhoWTR_XY_D1_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRhoWTR_XY_D1_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRhoWTR_XY_D1_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRhoWTR_XY_D1_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRhoWTR_XY_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRhoWTR_XY_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRhoWTR_XY_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRhoWTR_XY_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRhoWTR_XY_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTR_XY_D1_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTR_XY_D1_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTR_XY_D1_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTR_XY_D1_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTR_XY_D1_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRhoWTR_XY_D1_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XY_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRho_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRho_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWNO_UXD_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWNO_UXD_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWNO_UXD_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWNO_UXD_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWNO_UXD_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRhoWNO_UXD_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRhoWNO_UXD_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRhoWNO_UXD_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRhoWNO_UXD_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRhoWNO_UXD_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRhoWNO_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRhoWNO_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRhoWNO_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRhoWNO_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRhoWNO_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWNO_UXD_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWNO_UXD_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWNO_UXD_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWNO_UXD_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWNO_UXD_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRhoWNO_UXD_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTI_UXD_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTI_UXD_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTI_UXD_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTI_UXD_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTI_UXD_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRhoWTI_UXD_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRhoWTI_UXD_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRhoWTI_UXD_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRhoWTI_UXD_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRhoWTI_UXD_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRhoWTI_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRhoWTI_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRhoWTI_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRhoWTI_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRhoWTI_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTI_UXD_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTI_UXD_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTI_UXD_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTI_UXD_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTI_UXD_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRhoWTI_UXD_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTR_UXD_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTR_UXD_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTR_UXD_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTR_UXD_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTR_UXD_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRhoWTR_UXD_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRhoWTR_UXD_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRhoWTR_UXD_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRhoWTR_UXD_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRhoWTR_UXD_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRhoWTR_UXD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRhoWTR_UXD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRhoWTR_UXD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRhoWTR_UXD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRhoWTR_UXD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTR_UXD_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTR_UXD_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTR_UXD_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTR_UXD_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTR_UXD_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRhoWTR_UXD_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRho_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRho_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WNO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWNO_XLD_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWNO_XLD_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWNO_XLD_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWNO_XLD_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWNO_XLD_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRhoWNO_XLD_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRhoWNO_XLD_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRhoWNO_XLD_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRhoWNO_XLD_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRhoWNO_XLD_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRhoWNO_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRhoWNO_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRhoWNO_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRhoWNO_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRhoWNO_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWNO_XLD_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWNO_XLD_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWNO_XLD_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWNO_XLD_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWNO_XLD_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRhoWNO_XLD_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WNO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTI_XLD_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTI_XLD_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTI_XLD_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTI_XLD_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTI_XLD_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRhoWTI_XLD_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRhoWTI_XLD_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRhoWTI_XLD_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRhoWTI_XLD_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRhoWTI_XLD_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRhoWTI_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRhoWTI_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRhoWTI_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRhoWTI_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRhoWTI_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTI_XLD_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTI_XLD_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTI_XLD_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTI_XLD_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTI_XLD_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRhoWTI_XLD_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define WTR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTR_XLD_SK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTR_XLD_SK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTR_XLD_SK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTR_XLD_SK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTR_XLD_SK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setRhoWTR_XLD_IK5
        use pm_kind, only: TKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setRhoWTR_XLD_IK4
        use pm_kind, only: TKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setRhoWTR_XLD_IK3
        use pm_kind, only: TKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setRhoWTR_XLD_IK2
        use pm_kind, only: TKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setRhoWTR_XLD_IK1
        use pm_kind, only: TKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRhoWTR_XLD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRhoWTR_XLD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRhoWTR_XLD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRhoWTR_XLD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRhoWTR_XLD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setRhoWTR_XLD_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setRhoWTR_XLD_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setRhoWTR_XLD_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setRhoWTR_XLD_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setRhoWTR_XLD_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setRhoWTR_XLD_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef WTR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRho_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getTau_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define XY_D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WNO_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWNO_XY_D0_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWNO_XY_D0_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWNO_XY_D0_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWNO_XY_D0_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWNO_XY_D0_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WNO_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTI_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWTI_XY_D0_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWTI_XY_D0_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWTI_XY_D0_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWTI_XY_D0_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWTI_XY_D0_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTI_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTR_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWTR_XY_D0_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWTR_XY_D0_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWTR_XY_D0_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWTR_XY_D0_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWTR_XY_D0_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTR_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef XY_D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef getTau_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getTau_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define XY_D1_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WNO_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWNO_XY_D1_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWNO_XY_D1_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWNO_XY_D1_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWNO_XY_D1_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWNO_XY_D1_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure getTauWNO_XY_D1_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure getTauWNO_XY_D1_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure getTauWNO_XY_D1_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure getTauWNO_XY_D1_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure getTauWNO_XY_D1_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getTauWNO_XY_D1_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getTauWNO_XY_D1_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getTauWNO_XY_D1_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getTauWNO_XY_D1_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getTauWNO_XY_D1_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWNO_XY_D1_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWNO_XY_D1_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWNO_XY_D1_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWNO_XY_D1_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWNO_XY_D1_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure getTauWNO_XY_D1_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WNO_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTI_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWTI_XY_D1_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWTI_XY_D1_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWTI_XY_D1_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWTI_XY_D1_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWTI_XY_D1_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure getTauWTI_XY_D1_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure getTauWTI_XY_D1_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure getTauWTI_XY_D1_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure getTauWTI_XY_D1_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure getTauWTI_XY_D1_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getTauWTI_XY_D1_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getTauWTI_XY_D1_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getTauWTI_XY_D1_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getTauWTI_XY_D1_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getTauWTI_XY_D1_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWTI_XY_D1_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWTI_XY_D1_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWTI_XY_D1_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWTI_XY_D1_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWTI_XY_D1_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure getTauWTI_XY_D1_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTI_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTR_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWTR_XY_D1_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWTR_XY_D1_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWTR_XY_D1_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWTR_XY_D1_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWTR_XY_D1_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure getTauWTR_XY_D1_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure getTauWTR_XY_D1_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure getTauWTR_XY_D1_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure getTauWTR_XY_D1_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure getTauWTR_XY_D1_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getTauWTR_XY_D1_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getTauWTR_XY_D1_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getTauWTR_XY_D1_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getTauWTR_XY_D1_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getTauWTR_XY_D1_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWTR_XY_D1_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWTR_XY_D1_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWTR_XY_D1_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWTR_XY_D1_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWTR_XY_D1_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure getTauWTR_XY_D1_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTR_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef XY_D1_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef getTau_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getTau_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define ULD_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WNO_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWNO_ULD_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWNO_ULD_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWNO_ULD_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWNO_ULD_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWNO_ULD_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure getTauWNO_ULD_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure getTauWNO_ULD_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure getTauWNO_ULD_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure getTauWNO_ULD_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure getTauWNO_ULD_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getTauWNO_ULD_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getTauWNO_ULD_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getTauWNO_ULD_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getTauWNO_ULD_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getTauWNO_ULD_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWNO_ULD_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWNO_ULD_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWNO_ULD_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWNO_ULD_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWNO_ULD_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure getTauWNO_ULD_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WNO_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTI_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWTI_ULD_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWTI_ULD_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWTI_ULD_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWTI_ULD_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWTI_ULD_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure getTauWTI_ULD_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure getTauWTI_ULD_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure getTauWTI_ULD_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure getTauWTI_ULD_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure getTauWTI_ULD_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getTauWTI_ULD_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getTauWTI_ULD_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getTauWTI_ULD_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getTauWTI_ULD_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getTauWTI_ULD_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWTI_ULD_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWTI_ULD_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWTI_ULD_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWTI_ULD_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWTI_ULD_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure getTauWTI_ULD_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTI_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTR_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWTR_ULD_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWTR_ULD_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWTR_ULD_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWTR_ULD_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWTR_ULD_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure getTauWTR_ULD_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure getTauWTR_ULD_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure getTauWTR_ULD_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure getTauWTR_ULD_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure getTauWTR_ULD_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getTauWTR_ULD_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getTauWTR_ULD_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getTauWTR_ULD_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getTauWTR_ULD_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getTauWTR_ULD_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure getTauWTR_ULD_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure getTauWTR_ULD_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure getTauWTR_ULD_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure getTauWTR_ULD_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure getTauWTR_ULD_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure getTauWTR_ULD_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTR_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef ULD_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef getTau_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define setTau_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define XY_D0_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WNO_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWNO_XY_D0_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWNO_XY_D0_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWNO_XY_D0_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWNO_XY_D0_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWNO_XY_D0_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WNO_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTI_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTI_XY_D0_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTI_XY_D0_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTI_XY_D0_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTI_XY_D0_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTI_XY_D0_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTI_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTR_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTR_XY_D0_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTR_XY_D0_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTR_XY_D0_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTR_XY_D0_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTR_XY_D0_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTR_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef XY_D0_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef setTau_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define setTau_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define XY_D1_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WNO_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWNO_XY_D1_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWNO_XY_D1_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWNO_XY_D1_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWNO_XY_D1_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWNO_XY_D1_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setTauWNO_XY_D1_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setTauWNO_XY_D1_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setTauWNO_XY_D1_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setTauWNO_XY_D1_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setTauWNO_XY_D1_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setTauWNO_XY_D1_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setTauWNO_XY_D1_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setTauWNO_XY_D1_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setTauWNO_XY_D1_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setTauWNO_XY_D1_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWNO_XY_D1_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWNO_XY_D1_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWNO_XY_D1_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWNO_XY_D1_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWNO_XY_D1_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure setTauWNO_XY_D1_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WNO_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTI_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTI_XY_D1_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTI_XY_D1_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTI_XY_D1_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTI_XY_D1_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTI_XY_D1_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setTauWTI_XY_D1_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setTauWTI_XY_D1_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setTauWTI_XY_D1_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setTauWTI_XY_D1_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setTauWTI_XY_D1_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setTauWTI_XY_D1_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setTauWTI_XY_D1_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setTauWTI_XY_D1_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setTauWTI_XY_D1_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setTauWTI_XY_D1_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTI_XY_D1_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTI_XY_D1_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTI_XY_D1_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTI_XY_D1_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTI_XY_D1_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure setTauWTI_XY_D1_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTI_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTR_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTR_XY_D1_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTR_XY_D1_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTR_XY_D1_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTR_XY_D1_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTR_XY_D1_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setTauWTR_XY_D1_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setTauWTR_XY_D1_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setTauWTR_XY_D1_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setTauWTR_XY_D1_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setTauWTR_XY_D1_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setTauWTR_XY_D1_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setTauWTR_XY_D1_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setTauWTR_XY_D1_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setTauWTR_XY_D1_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setTauWTR_XY_D1_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTR_XY_D1_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTR_XY_D1_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTR_XY_D1_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTR_XY_D1_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTR_XY_D1_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure setTauWTR_XY_D1_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTR_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef XY_D1_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef setTau_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define setTau_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define UXD_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WNO_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWNO_UXD_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWNO_UXD_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWNO_UXD_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWNO_UXD_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWNO_UXD_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setTauWNO_UXD_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setTauWNO_UXD_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setTauWNO_UXD_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setTauWNO_UXD_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setTauWNO_UXD_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setTauWNO_UXD_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setTauWNO_UXD_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setTauWNO_UXD_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setTauWNO_UXD_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setTauWNO_UXD_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWNO_UXD_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWNO_UXD_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWNO_UXD_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWNO_UXD_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWNO_UXD_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure setTauWNO_UXD_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WNO_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTI_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTI_UXD_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTI_UXD_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTI_UXD_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTI_UXD_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTI_UXD_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setTauWTI_UXD_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setTauWTI_UXD_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setTauWTI_UXD_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setTauWTI_UXD_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setTauWTI_UXD_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setTauWTI_UXD_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setTauWTI_UXD_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setTauWTI_UXD_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setTauWTI_UXD_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setTauWTI_UXD_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTI_UXD_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTI_UXD_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTI_UXD_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTI_UXD_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTI_UXD_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure setTauWTI_UXD_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTI_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTR_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTR_UXD_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTR_UXD_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTR_UXD_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTR_UXD_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTR_UXD_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setTauWTR_UXD_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setTauWTR_UXD_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setTauWTR_UXD_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setTauWTR_UXD_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setTauWTR_UXD_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setTauWTR_UXD_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setTauWTR_UXD_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setTauWTR_UXD_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setTauWTR_UXD_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setTauWTR_UXD_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTR_UXD_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTR_UXD_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTR_UXD_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTR_UXD_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTR_UXD_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure setTauWTR_UXD_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTR_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef UXD_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef setTau_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define setTau_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define XLD_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WNO_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWNO_XLD_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWNO_XLD_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWNO_XLD_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWNO_XLD_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWNO_XLD_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setTauWNO_XLD_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setTauWNO_XLD_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setTauWNO_XLD_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setTauWNO_XLD_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setTauWNO_XLD_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setTauWNO_XLD_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setTauWNO_XLD_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setTauWNO_XLD_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setTauWNO_XLD_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setTauWNO_XLD_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWNO_XLD_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWNO_XLD_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWNO_XLD_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWNO_XLD_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWNO_XLD_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure setTauWNO_XLD_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WNO_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTI_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTI_XLD_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTI_XLD_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTI_XLD_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTI_XLD_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTI_XLD_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setTauWTI_XLD_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setTauWTI_XLD_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setTauWTI_XLD_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setTauWTI_XLD_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setTauWTI_XLD_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setTauWTI_XLD_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setTauWTI_XLD_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setTauWTI_XLD_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setTauWTI_XLD_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setTauWTI_XLD_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTI_XLD_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTI_XLD_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTI_XLD_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTI_XLD_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTI_XLD_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure setTauWTI_XLD_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTI_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define WTR_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define SK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTR_XLD_SK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTR_XLD_SK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTR_XLD_SK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTR_XLD_SK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTR_XLD_SK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef SK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define IK_ENABLED 1
!
!#if IK5_ENABLED
!    module procedure setTauWTR_XLD_IK5
!        use pm_kind, only: TKC => IK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK4_ENABLED
!    module procedure setTauWTR_XLD_IK4
!        use pm_kind, only: TKC => IK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK3_ENABLED
!    module procedure setTauWTR_XLD_IK3
!        use pm_kind, only: TKC => IK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK2_ENABLED
!    module procedure setTauWTR_XLD_IK2
!        use pm_kind, only: TKC => IK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if IK1_ENABLED
!    module procedure setTauWTR_XLD_IK1
!        use pm_kind, only: TKC => IK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef IK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure setTauWTR_XLD_RK5
!        use pm_kind, only: TKC => RK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure setTauWTR_XLD_RK4
!        use pm_kind, only: TKC => RK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure setTauWTR_XLD_RK3
!        use pm_kind, only: TKC => RK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure setTauWTR_XLD_RK2
!        use pm_kind, only: TKC => RK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure setTauWTR_XLD_RK1
!        use pm_kind, only: TKC => RK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if !__GFORTRAN__
!#define PSSK_ENABLED 1
!
!#if SK5_ENABLED
!    module procedure setTauWTR_XLD_PSSK5
!        use pm_kind, only: TKC => SK5
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK4_ENABLED
!    module procedure setTauWTR_XLD_PSSK4
!        use pm_kind, only: TKC => SK4
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK3_ENABLED
!    module procedure setTauWTR_XLD_PSSK3
!        use pm_kind, only: TKC => SK3
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK2_ENABLED
!    module procedure setTauWTR_XLD_PSSK2
!        use pm_kind, only: TKC => SK2
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#if SK1_ENABLED
!    module procedure setTauWTR_XLD_PSSK1
!        use pm_kind, only: TKC => SK1
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!#endif
!
!#undef PSSK_ENABLED
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define BSSK_ENABLED 1
!
!    module procedure setTauWTR_XLD_BSSK
!        use pm_kind, only: TKC => SK
!#include "pm_sampleCor@routines.inc.F90"
!    end procedure
!
!#undef BSSK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef WTR_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef XLD_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef setTau_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCordance_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Sum_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setCordanceSum_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setCordanceSum_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setCordanceSum_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setCordanceSum_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setCordanceSum_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sampleCor@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setCordanceSum_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setCordanceSum_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setCordanceSum_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setCordanceSum_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setCordanceSum_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCordanceSum_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCordanceSum_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCordanceSum_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCordanceSum_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCordanceSum_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCordanceSum_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCordanceSum_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCordanceSum_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCordanceSum_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCordanceSum_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setCordanceSum_D1_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setCordanceSum_D1_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setCordanceSum_D1_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setCordanceSum_D1_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setCordanceSum_D1_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setCordanceSum_D1_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Sum_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define All_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setCordanceAll_D0_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setCordanceAll_D0_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setCordanceAll_D0_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setCordanceAll_D0_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setCordanceAll_D0_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sampleCor@routines.inc.F90"
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

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setCordanceAll_D1_SK5
        use pm_kind, only: SKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setCordanceAll_D1_SK4
        use pm_kind, only: SKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setCordanceAll_D1_SK3
        use pm_kind, only: SKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setCordanceAll_D1_SK2
        use pm_kind, only: SKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setCordanceAll_D1_SK1
        use pm_kind, only: SKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setCordanceAll_D1_IK5
        use pm_kind, only: IKC => IK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setCordanceAll_D1_IK4
        use pm_kind, only: IKC => IK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setCordanceAll_D1_IK3
        use pm_kind, only: IKC => IK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setCordanceAll_D1_IK2
        use pm_kind, only: IKC => IK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setCordanceAll_D1_IK1
        use pm_kind, only: IKC => IK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCordanceAll_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCordanceAll_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCordanceAll_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCordanceAll_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCordanceAll_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !__GFORTRAN__
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure setCordanceAll_D1_PSSK5
        use pm_kind, only: TKC => SK5
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setCordanceAll_D1_PSSK4
        use pm_kind, only: TKC => SK4
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setCordanceAll_D1_PSSK3
        use pm_kind, only: TKC => SK3
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setCordanceAll_D1_PSSK2
        use pm_kind, only: TKC => SK2
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setCordanceAll_D1_PSSK1
        use pm_kind, only: TKC => SK1
#include "pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure setCordanceAll_D1_BSSK
        use pm_kind, only: TKC => SK
#include "pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef All_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCordance_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines