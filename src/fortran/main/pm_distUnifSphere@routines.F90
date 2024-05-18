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
!>  This file contains procedure implementations of [pm_distUnifSphere](@ref pm_distUnifSphere).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distUnifSphere) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_matrixDet, only: transHerm
    use pm_matrixDet, only: setMatDetSqrtLog
    use pm_matrixCopy, only: setMatCopy, rdpack
    use pm_matrixTrace, only: getMatMulTraceLog
    use pm_ellipsoid, only: getLogVolUnitBall
    use pm_distUnif, only: setUnifRand
    use pm_distNorm, only: setNormRand

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getUnifSphereLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getUnifSphereLogPDF_D0_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUnifSphereLogPDF_D0_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUnifSphereLogPDF_D0_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUnifSphereLogPDF_D0_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUnifSphereLogPDF_D0_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getUnifSphereLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getUSR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

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
    module procedure getUSR_RNGD_DM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGD_DM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGD_DM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGD_DM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGD_DM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGD_DM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGD_DM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGD_DM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGD_DM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGD_DM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure getUSR_RNGD_AM_DC_XXX_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGD_AM_DC_XXX_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGD_AM_DC_XXX_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGD_AM_DC_XXX_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGD_AM_DC_XXX_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGD_AM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGD_AM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGD_AM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGD_AM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGD_AM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGD_AM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGD_AM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGD_AM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGD_AM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGD_AM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

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
    module procedure getUSR_RNGF_DM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGF_DM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGF_DM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGF_DM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGF_DM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGF_DM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGF_DM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGF_DM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGF_DM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGF_DM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure getUSR_RNGF_AM_DC_XXX_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGF_AM_DC_XXX_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGF_AM_DC_XXX_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGF_AM_DC_XXX_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGF_AM_DC_XXX_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGF_AM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGF_AM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGF_AM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGF_AM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGF_AM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGF_AM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGF_AM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGF_AM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGF_AM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGF_AM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#define DM_ENABLED 1

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
    module procedure getUSR_RNGX_DM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGX_DM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGX_DM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGX_DM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGX_DM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGX_DM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGX_DM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGX_DM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGX_DM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGX_DM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure getUSR_RNGX_AM_DC_XXX_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGX_AM_DC_XXX_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGX_AM_DC_XXX_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGX_AM_DC_XXX_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGX_AM_DC_XXX_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGX_AM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGX_AM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGX_AM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGX_AM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGX_AM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGX_AM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGX_AM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGX_AM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGX_AM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGX_AM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

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
    module procedure getUSR_RNGD_DM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGD_DM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGD_DM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGD_DM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGD_DM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGD_DM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGD_DM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGD_DM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGD_DM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGD_DM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure getUSR_RNGD_AM_DC_XXX_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGD_AM_DC_XXX_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGD_AM_DC_XXX_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGD_AM_DC_XXX_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGD_AM_DC_XXX_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGD_AM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGD_AM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGD_AM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGD_AM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGD_AM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGD_AM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGD_AM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGD_AM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGD_AM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGD_AM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

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
    module procedure getUSR_RNGF_DM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGF_DM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGF_DM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGF_DM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGF_DM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGF_DM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGF_DM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGF_DM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGF_DM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGF_DM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure getUSR_RNGF_AM_DC_XXX_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGF_AM_DC_XXX_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGF_AM_DC_XXX_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGF_AM_DC_XXX_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGF_AM_DC_XXX_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGF_AM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGF_AM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGF_AM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGF_AM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGF_AM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGF_AM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGF_AM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGF_AM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGF_AM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGF_AM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#define DM_ENABLED 1

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
    module procedure getUSR_RNGX_DM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGX_DM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGX_DM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGX_DM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGX_DM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGX_DM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGX_DM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGX_DM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGX_DM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGX_DM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure getUSR_RNGX_AM_DC_XXX_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGX_AM_DC_XXX_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGX_AM_DC_XXX_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGX_AM_DC_XXX_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGX_AM_DC_XXX_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGX_AM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGX_AM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGX_AM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGX_AM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGX_AM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure getUSR_RNGX_AM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getUSR_RNGX_AM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getUSR_RNGX_AM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getUSR_RNGX_AM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getUSR_RNGX_AM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef getUSR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setUSR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUSR_RNGD_DM_DC_XXX_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_DM_DC_XXX_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_DM_DC_XXX_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_DM_DC_XXX_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_DM_DC_XXX_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGD_DM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_DM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_DM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_DM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_DM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGD_DM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_DM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_DM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_DM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_DM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure setUSR_RNGD_AM_DC_XXX_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_AM_DC_XXX_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_AM_DC_XXX_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_AM_DC_XXX_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_AM_DC_XXX_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGD_AM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_AM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_AM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_AM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_AM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGD_AM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_AM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_AM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_AM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_AM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUSR_RNGF_DM_DC_XXX_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_DM_DC_XXX_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_DM_DC_XXX_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_DM_DC_XXX_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_DM_DC_XXX_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGF_DM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_DM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_DM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_DM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_DM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGF_DM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_DM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_DM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_DM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_DM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure setUSR_RNGF_AM_DC_XXX_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_AM_DC_XXX_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_AM_DC_XXX_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_AM_DC_XXX_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_AM_DC_XXX_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGF_AM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_AM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_AM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_AM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_AM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGF_AM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_AM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_AM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_AM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_AM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUSR_RNGX_DM_DC_XXX_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_DM_DC_XXX_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_DM_DC_XXX_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_DM_DC_XXX_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_DM_DC_XXX_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGX_DM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_DM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_DM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_DM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_DM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGX_DM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_DM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_DM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_DM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_DM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure setUSR_RNGX_AM_DC_XXX_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_AM_DC_XXX_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_AM_DC_XXX_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_AM_DC_XXX_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_AM_DC_XXX_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGX_AM_AC_UXD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_AM_AC_UXD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_AM_AC_UXD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_AM_AC_UXD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_AM_AC_UXD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGX_AM_AC_XLD_D1_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_AM_AC_XLD_D1_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_AM_AC_XLD_D1_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_AM_AC_XLD_D1_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_AM_AC_XLD_D1_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUSR_RNGD_DM_DC_XXX_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_DM_DC_XXX_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_DM_DC_XXX_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_DM_DC_XXX_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_DM_DC_XXX_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGD_DM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_DM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_DM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_DM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_DM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGD_DM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_DM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_DM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_DM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_DM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure setUSR_RNGD_AM_DC_XXX_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_AM_DC_XXX_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_AM_DC_XXX_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_AM_DC_XXX_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_AM_DC_XXX_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGD_AM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_AM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_AM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_AM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_AM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGD_AM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGD_AM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGD_AM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGD_AM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGD_AM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUSR_RNGF_DM_DC_XXX_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_DM_DC_XXX_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_DM_DC_XXX_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_DM_DC_XXX_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_DM_DC_XXX_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGF_DM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_DM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_DM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_DM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_DM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGF_DM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_DM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_DM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_DM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_DM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure setUSR_RNGF_AM_DC_XXX_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_AM_DC_XXX_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_AM_DC_XXX_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_AM_DC_XXX_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_AM_DC_XXX_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGF_AM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_AM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_AM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_AM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_AM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGF_AM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGF_AM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGF_AM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGF_AM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGF_AM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setUSR_RNGX_DM_DC_XXX_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_DM_DC_XXX_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_DM_DC_XXX_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_DM_DC_XXX_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_DM_DC_XXX_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGX_DM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_DM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_DM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_DM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_DM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGX_DM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_DM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_DM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_DM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_DM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef DM_ENABLED

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
    module procedure setUSR_RNGX_AM_DC_XXX_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_AM_DC_XXX_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_AM_DC_XXX_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_AM_DC_XXX_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_AM_DC_XXX_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGX_AM_AC_UXD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_AM_AC_UXD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_AM_AC_UXD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_AM_AC_UXD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_AM_AC_UXD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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
    module procedure setUSR_RNGX_AM_AC_XLD_D2_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setUSR_RNGX_AM_AC_XLD_D2_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setUSR_RNGX_AM_AC_XLD_D2_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setUSR_RNGX_AM_AC_XLD_D2_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distUnifSphere@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setUSR_RNGX_AM_AC_XLD_D2_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distUnifSphere@routines.inc.F90"
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

#undef setUSR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines