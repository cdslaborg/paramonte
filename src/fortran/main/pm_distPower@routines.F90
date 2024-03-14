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
!>  This file contains procedure implementations of [pm_distPower](@ref pm_distPower).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distPower) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_mathMinMax, only: getMinMax
    use pm_mathLogSubExp, only: getLogSubExp
    use pm_mathLogAddExp, only: getLogAddExp
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowerLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowerLogPDFNFALD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogPDFNFALD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogPDFNFALD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogPDFNFALD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogPDFNFALD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
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
    module procedure getPowerLogPDFNFALL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogPDFNFALL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogPDFNFALL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogPDFNFALL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogPDFNFALL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowerLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowerLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowerLogPDFALD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogPDFALD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogPDFALD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogPDFALD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogPDFALD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
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
    module procedure getPowerLogPDFALL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogPDFALL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogPDFALL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogPDFALL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogPDFALL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowerLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowerLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPowerLogPDF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowerLogPDF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowerLogPDF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowerLogPDF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowerLogPDF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowerLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowerLogCDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowerLogCDFNFALD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogCDFNFALD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogCDFNFALD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogCDFNFALD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogCDFNFALD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
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
    module procedure getPowerLogCDFNFALL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogCDFNFALL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogCDFNFALL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogCDFNFALL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogCDFNFALL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowerLogCDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowerLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowerLogCDFALD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogCDFALD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogCDFALD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogCDFALD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogCDFALD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
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
    module procedure getPowerLogCDFALL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogCDFALL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogCDFALL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogCDFALL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogCDFALL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowerLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowerLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPowerLogCDFALD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowerLogCDFALD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowerLogCDFALD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowerLogCDFALD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowerLogCDFALD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
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
    module procedure setPowerLogCDFALL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowerLogCDFALL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowerLogCDFALL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowerLogCDFALL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowerLogCDFALL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowerLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowerLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowerLogQuanALD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogQuanALD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogQuanALD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogQuanALD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogQuanALD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
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
    module procedure getPowerLogQuanALL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogQuanALL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogQuanALL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogQuanALL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogQuanALL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowerLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowerLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LLALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPowerLogQuanLLALD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowerLogQuanLLALD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowerLogQuanLLALD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowerLogQuanLLALD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowerLogQuanLLALD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
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
    module procedure setPowerLogQuanLLALL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowerLogQuanLLALL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowerLogQuanLLALL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowerLogQuanLLALL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowerLogQuanLLALL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LLALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowerLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowerLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getPowerLogRandALD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogRandALD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogRandALD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogRandALD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogRandALD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
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
    module procedure getPowerLogRandALL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getPowerLogRandALL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getPowerLogRandALL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getPowerLogRandALL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getPowerLogRandALL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowerLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowerLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LNALD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setPowerLogRandLNALD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowerLogRandLNALD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowerLogRandLNALD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowerLogRandLNALD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowerLogRandLNALD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
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
    module procedure setPowerLogRandLNALL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setPowerLogRandLNALL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setPowerLogRandLNALL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setPowerLogRandLNALL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setPowerLogRandLNALL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LNALL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowerLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines