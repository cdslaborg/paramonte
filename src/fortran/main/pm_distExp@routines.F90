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
!>  This file contains procedure implementations of [pm_distExp](@ref pm_distExp).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distExp) routines ! LCOV_EXCL_LINE

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
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XMI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getExpLogPDFXMI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExpLogPDFXMI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExpLogPDFXMI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExpLogPDFXMI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExpLogPDFXMI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setExpLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExpLogPDFXDD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpLogPDFXDD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpLogPDFXDD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpLogPDFXDD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpLogPDFXDD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XMD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExpLogPDFXMD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpLogPDFXMD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpLogPDFXMD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpLogPDFXMD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpLogPDFXMD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XDI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExpLogPDFXDI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpLogPDFXDI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpLogPDFXDI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpLogPDFXDI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpLogPDFXDI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
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
    module procedure setExpLogPDFXMI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpLogPDFXMI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpLogPDFXMI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpLogPDFXMI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpLogPDFXMI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setExpLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XMI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getExpCDFXMI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExpCDFXMI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExpCDFXMI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExpCDFXMI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExpCDFXMI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setExpCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExpCDFXDD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpCDFXDD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpCDFXDD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpCDFXDD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpCDFXDD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
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
    module procedure setExpCDFXDI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpCDFXDI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpCDFXDI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpCDFXDI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpCDFXDI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
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
    module procedure setExpCDFXMI_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpCDFXMI_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpCDFXMI_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpCDFXMI_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpCDFXMI_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setExpCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getExpRandSM_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExpRandSM_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExpRandSM_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExpRandSM_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExpRandSM_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setExpRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExpRandDD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpRandDD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpRandDD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpRandDD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpRandDD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExpRandSD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpRandSD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpRandSD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpRandSD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpRandSD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExpRandSM_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpRandSM_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpRandSM_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpRandSM_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpRandSM_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setExpRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines