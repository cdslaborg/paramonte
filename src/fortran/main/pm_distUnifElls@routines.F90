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
!>  This file contains procedure implementations of [pm_distUnifElls](@ref pm_distUnifElls).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distUnifElls) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif
    use pm_matrixTrace, only: getMatMulTraceLog
    use pm_mathCumPropExp, only: setCumPropExp, sequence
    use pm_ellipsoid, only: getLogVolUnitBall
    use pm_distanceMahal, only: setDisMahalSq
    use pm_mathLogSumExp, only: getLogSumExp
    use pm_distUnifEll, only: setUnifEllRand
    use pm_arraySpace, only: setLinSpace
    use pm_ellipsoid, only: isMemberEll
    use pm_distUnif, only: setUnifRand

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getUnifEllsLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMMUP_RNGF_AM_DC_XXX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMMUP_RNGF_AM_DC_XXX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMMUP_RNGF_AM_DC_XXX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMMUP_RNGF_AM_DC_XXX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMMUP_RNGF_AM_DC_XXX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMMUP_RNGF_AM_AC_UXD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMMUP_RNGF_AM_AC_UXD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMMUP_RNGF_AM_AC_UXD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMMUP_RNGF_AM_AC_UXD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMMUP_RNGF_AM_AC_UXD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
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

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMMUP_RNGF_AM_AC_XLD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMMUP_RNGF_AM_AC_XLD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMMUP_RNGF_AM_AC_XLD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMMUP_RNGF_AM_AC_XLD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMMUP_RNGF_AM_AC_XLD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMMUP_RNGX_AM_DC_XXX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMMUP_RNGX_AM_DC_XXX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMMUP_RNGX_AM_DC_XXX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMMUP_RNGX_AM_DC_XXX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMMUP_RNGX_AM_DC_XXX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMMUP_RNGX_AM_AC_UXD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMMUP_RNGX_AM_AC_UXD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMMUP_RNGX_AM_AC_UXD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMMUP_RNGX_AM_AC_UXD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMMUP_RNGX_AM_AC_UXD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
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

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMMUP_RNGX_AM_AC_XLD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMMUP_RNGX_AM_AC_XLD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMMUP_RNGX_AM_AC_XLD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMMUP_RNGX_AM_AC_XLD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMMUP_RNGX_AM_AC_XLD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getUnifEllsLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMMUR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGF_AM_DC_XXX_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGF_AM_DC_XXX_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGF_AM_DC_XXX_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGF_AM_DC_XXX_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGF_AM_DC_XXX_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGF_AM_AC_UXD_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGF_AM_AC_UXD_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGF_AM_AC_UXD_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGF_AM_AC_UXD_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGF_AM_AC_UXD_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
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

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGF_AM_AC_XLD_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGF_AM_AC_XLD_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGF_AM_AC_XLD_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGF_AM_AC_XLD_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGF_AM_AC_XLD_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGX_AM_DC_XXX_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGX_AM_DC_XXX_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGX_AM_DC_XXX_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGX_AM_DC_XXX_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGX_AM_DC_XXX_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGX_AM_AC_UXD_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGX_AM_AC_UXD_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGX_AM_AC_UXD_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGX_AM_AC_UXD_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGX_AM_AC_UXD_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
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

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGX_AM_AC_XLD_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGX_AM_AC_XLD_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGX_AM_AC_XLD_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGX_AM_AC_XLD_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGX_AM_AC_XLD_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGF_AM_DC_XXX_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGF_AM_DC_XXX_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGF_AM_DC_XXX_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGF_AM_DC_XXX_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGF_AM_DC_XXX_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGF_AM_AC_UXD_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGF_AM_AC_UXD_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGF_AM_AC_UXD_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGF_AM_AC_UXD_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGF_AM_AC_UXD_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
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

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGF_AM_AC_XLD_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGF_AM_AC_XLD_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGF_AM_AC_XLD_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGF_AM_AC_XLD_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGF_AM_AC_XLD_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGX_AM_DC_XXX_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGX_AM_DC_XXX_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGX_AM_DC_XXX_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGX_AM_DC_XXX_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGX_AM_DC_XXX_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGX_AM_AC_UXD_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGX_AM_AC_UXD_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGX_AM_AC_UXD_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGX_AM_AC_UXD_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGX_AM_AC_UXD_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
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

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMMUR_RNGX_AM_AC_XLD_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMMUR_RNGX_AM_AC_XLD_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMMUR_RNGX_AM_AC_XLD_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMMUR_RNGX_AM_AC_XLD_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMMUR_RNGX_AM_AC_XLD_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distUnifElls@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMMUR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines