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
!>  This file contains procedure implementations of [pm_distNegExp](@ref pm_distNegExp).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distNegExp) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNegExpLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XMI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNegExpLogPDFXMI_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNegExpLogPDFXMI_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNegExpLogPDFXMI_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNegExpLogPDFXMI_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNegExpLogPDFXMI_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNegExpLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNegExpLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNegExpLogPDFDDD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNegExpLogPDFDDD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNegExpLogPDFDDD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNegExpLogPDFDDD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNegExpLogPDFDDD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNegExpLogPDFXMDD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNegExpLogPDFXMDD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNegExpLogPDFXMDD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNegExpLogPDFXMDD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNegExpLogPDFXMDD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNegExpLogPDFDIL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNegExpLogPDFDIL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNegExpLogPDFDIL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNegExpLogPDFDIL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNegExpLogPDFDIL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DIL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MIL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNegExpLogPDFMIL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNegExpLogPDFMIL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNegExpLogPDFMIL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNegExpLogPDFMIL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNegExpLogPDFMIL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MIL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNegExpLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNegExpCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XMI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNegExpCDFXMI_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNegExpCDFXMI_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNegExpCDFXMI_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNegExpCDFXMI_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNegExpCDFXMI_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNegExpCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNegExpCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNegExpCDFXDD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNegExpCDFXDD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNegExpCDFXDD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNegExpCDFXDD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNegExpCDFXDD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XDI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNegExpCDFXDI_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNegExpCDFXDI_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNegExpCDFXDI_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNegExpCDFXDI_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNegExpCDFXDI_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XDI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XMI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNegExpCDFXMI_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNegExpCDFXMI_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNegExpCDFXMI_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNegExpCDFXMI_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNegExpCDFXMI_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNegExpCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNegExpRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNegExpRandSM_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNegExpRandSM_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNegExpRandSM_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNegExpRandSM_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNegExpRandSM_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNegExpRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNegExpRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNegExpRandUDD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNegExpRandUDD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNegExpRandUDD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNegExpRandUDD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNegExpRandUDD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define USD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNegExpRandUSD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNegExpRandUSD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNegExpRandUSD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNegExpRandUSD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNegExpRandUSD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef USD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define USM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNegExpRandUSM_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNegExpRandUSM_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNegExpRandUSM_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNegExpRandUSM_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNegExpRandUSM_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distNegExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef USM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNegExpRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines