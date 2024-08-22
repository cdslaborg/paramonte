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
!>  This file contains procedure implementations of [pm_distPois](@ref pm_distPois).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distPois) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_mathGamma, only: setGammaInc!Upp
    use pm_distUnif, only: setUnifRand
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPoisLogPMF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPoisLogPMF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPoisLogPMF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPoisLogPMF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPoisLogPMF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPoisLogPMF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPoisLogPMF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPoisLogPMF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPoisLogPMFDef_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisLogPMFDef_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisLogPMFDef_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisLogPMFDef_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisLogPMFDef_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisLogPMFLog_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisLogPMFLog_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisLogPMFLog_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisLogPMFLog_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisLogPMFLog_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Log_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPoisLogPMF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPoisCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPoisCDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPoisCDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPoisCDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPoisCDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPoisCDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPoisCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPoisCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Log_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPoisCDFLog_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisCDFLog_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisCDFLog_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisCDFLog_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisCDFLog_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Log_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPoisCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPoisRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPoisRand_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPoisRand_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPoisRand_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPoisRand_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPoisRand_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPoisRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPoisRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Exp_ENABLED 1

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
    module procedure setPoisRandExpRNGD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandExpRNGD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandExpRNGD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandExpRNGD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandExpRNGD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisRandExpRNGF_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandExpRNGF_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandExpRNGF_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandExpRNGF_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandExpRNGF_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisRandExpRNGX_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandExpRNGX_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandExpRNGX_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandExpRNGX_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandExpRNGX_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisRandExpRNGD_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandExpRNGD_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandExpRNGD_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandExpRNGD_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandExpRNGD_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisRandExpRNGF_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandExpRNGF_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandExpRNGF_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandExpRNGF_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandExpRNGF_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisRandExpRNGX_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandExpRNGX_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandExpRNGX_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandExpRNGX_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandExpRNGX_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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

#undef Exp_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rej_ENABLED 1

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
    module procedure setPoisRandRejRNGD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandRejRNGD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandRejRNGD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandRejRNGD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandRejRNGD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisRandRejRNGF_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandRejRNGF_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandRejRNGF_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandRejRNGF_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandRejRNGF_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisRandRejRNGX_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandRejRNGX_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandRejRNGX_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandRejRNGX_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandRejRNGX_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisRandRejRNGD_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandRejRNGD_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandRejRNGD_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandRejRNGD_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandRejRNGD_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisRandRejRNGF_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandRejRNGF_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandRejRNGF_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandRejRNGF_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandRejRNGF_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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
    module procedure setPoisRandRejRNGX_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPoisRandRejRNGX_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPoisRandRejRNGX_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPoisRandRejRNGX_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPois@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPoisRandRejRNGX_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPois@routines.inc.F90"
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

#undef Rej_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPoisRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines