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
!>  This file contains procedure implementations of [pm_distPareto](@ref pm_distPareto).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distPareto) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_distNegExp, only: getNegExpRand
    use pm_distUnif, only: getUnifRand
    use pm_mathLogSubExp, only: getLogSubExp
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getParetoLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogPDFNFALD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogPDFNFALD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogPDFNFALD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogPDFNFALD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogPDFNFALD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogPDFNFALL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogPDFNFALL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogPDFNFALL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogPDFNFALL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogPDFNFALL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getParetoLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getParetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogPDFALD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogPDFALD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogPDFALD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogPDFALD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogPDFALD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogPDFALL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogPDFALL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogPDFALL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogPDFALL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogPDFALL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getParetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setParetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setParetoLogPDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setParetoLogPDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setParetoLogPDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setParetoLogPDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setParetoLogPDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setParetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getParetoLogCDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogCDFNFALD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogCDFNFALD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogCDFNFALD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogCDFNFALD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogCDFNFALD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogCDFNFALL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogCDFNFALL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogCDFNFALL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogCDFNFALL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogCDFNFALL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getParetoLogCDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getParetoLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogCDFALD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogCDFALD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogCDFALD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogCDFALD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogCDFALD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogCDFALL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogCDFALL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogCDFALL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogCDFALL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogCDFALL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getParetoLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setParetoLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setParetoLogCDFALD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setParetoLogCDFALD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setParetoLogCDFALD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setParetoLogCDFALD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setParetoLogCDFALD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setParetoLogCDFALL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setParetoLogCDFALL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setParetoLogCDFALL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setParetoLogCDFALL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setParetoLogCDFALL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setParetoLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getParetoLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogQuanALD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogQuanALD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogQuanALD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogQuanALD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogQuanALD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogQuanALL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogQuanALL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogQuanALL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogQuanALL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogQuanALL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getParetoLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setParetoLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LLALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setParetoLogQuanLLALD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setParetoLogQuanLLALD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setParetoLogQuanLLALD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setParetoLogQuanLLALD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setParetoLogQuanLLALD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LLALD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LLALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setParetoLogQuanLLALL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setParetoLogQuanLLALL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setParetoLogQuanLLALL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setParetoLogQuanLLALL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setParetoLogQuanLLALL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LLALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setParetoLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getParetoLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogRandALD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogRandALD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogRandALD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogRandALD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogRandALD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getParetoLogRandALL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getParetoLogRandALL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getParetoLogRandALL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getParetoLogRandALL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getParetoLogRandALL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getParetoLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setParetoLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LNALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setParetoLogRandLNALD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setParetoLogRandLNALD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setParetoLogRandLNALD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setParetoLogRandLNALD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setParetoLogRandLNALD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LNALD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LNALL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setParetoLogRandLNALL_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setParetoLogRandLNALL_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setParetoLogRandLNALL_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setParetoLogRandLNALL_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setParetoLogRandLNALL_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LNALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setParetoLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
