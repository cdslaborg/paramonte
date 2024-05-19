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
!>  This file contains procedure implementations of [pm_distGeom](@ref pm_distGeom).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distGeom) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_distUnif, only: setUnifRand
    use pm_math1mexp, only: get1mexp
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGeomLogPMF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGeomLogPMF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGeomLogPMF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGeomLogPMF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGeomLogPMF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGeomLogPMF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGeomLogPMF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGeomLogPMF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomLogPMFDef_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomLogPMFDef_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomLogPMFDef_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomLogPMFDef_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomLogPMFDef_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Log_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomLogPMFLog_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomLogPMFLog_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomLogPMFLog_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomLogPMFLog_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomLogPMFLog_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Log_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGeomLogPMF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGeomCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGeomCDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGeomCDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGeomCDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGeomCDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGeomCDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGeomCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGeomCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomCDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomCDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomCDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomCDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomCDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGeomCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGeomRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGeomRand_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGeomRand_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGeomRand_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGeomRand_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGeomRand_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGeomRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGeomRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomRandRNGD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomRandRNGD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomRandRNGD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomRandRNGD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomRandRNGD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomRandRNGF_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomRandRNGF_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomRandRNGF_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomRandRNGF_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomRandRNGF_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomRandRNGX_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomRandRNGX_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomRandRNGX_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomRandRNGX_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomRandRNGX_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomRandRNGD_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomRandRNGD_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomRandRNGD_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomRandRNGD_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomRandRNGD_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomRandRNGF_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomRandRNGF_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomRandRNGF_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomRandRNGF_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomRandRNGF_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGeomRandRNGX_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGeomRandRNGX_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGeomRandRNGX_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGeomRandRNGX_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGeomRandRNGX_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGeom@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGeomRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
