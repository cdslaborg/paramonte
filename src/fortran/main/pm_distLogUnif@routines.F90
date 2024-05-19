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
!>  This file contains procedure implementations of [pm_distLogUnif](@ref pm_distLogUnif).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distLogUnif) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_mathLogSubExp, only: getLogSubExp
    use pm_distUnif, only: getUnifRand
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogUnifPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogUnifPDFNF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogUnifPDFNF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogUnifPDFNF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogUnifPDFNF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogUnifPDFNF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogUnifPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogUnifPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogUnifPDFMM_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogUnifPDFMM_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogUnifPDFMM_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogUnifPDFMM_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogUnifPDFMM_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogUnifPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setLogUnifPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogUnifPDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogUnifPDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogUnifPDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogUnifPDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogUnifPDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLogUnifPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogUnifCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogUnifCDFLL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogUnifCDFLL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogUnifCDFLL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogUnifCDFLL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogUnifCDFLL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogUnifCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setLogUnifCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogUnifCDFLL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogUnifCDFLL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogUnifCDFLL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogUnifCDFLL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogUnifCDFLL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLogUnifCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogUnifLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogUnifLogQuanLL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogUnifLogQuanLL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogUnifLogQuanLL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogUnifLogQuanLL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogUnifLogQuanLL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogUnifLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setLogUnifLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LLLP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogUnifLogQuanLLLP_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogUnifLogQuanLLLP_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogUnifLogQuanLLLP_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogUnifLogQuanLLLP_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogUnifLogQuanLLLP_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LLLP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLogUnifLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogUnifRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getLogUnifRandMM_IK5
        use pm_kind, only: IKG => IK5
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getLogUnifRandMM_IK4
        use pm_kind, only: IKG => IK4
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getLogUnifRandMM_IK3
        use pm_kind, only: IKG => IK3
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getLogUnifRandMM_IK2
        use pm_kind, only: IKG => IK2
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getLogUnifRandMM_IK1
        use pm_kind, only: IKG => IK1
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogUnifRandMM_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogUnifRandMM_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogUnifRandMM_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogUnifRandMM_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogUnifRandMM_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogUnifRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setLogUnifLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LLLP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogUnifLogRandLLLP_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogUnifLogRandLLLP_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogUnifLogRandLLLP_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogUnifLogRandLLLP_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogUnifLogRandLLLP_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogUnif@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LLLP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLogUnifLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines