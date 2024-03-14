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
!>  This file contains procedure implementations of [pm_distCosRaised](@ref pm_distCosRaised).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distCosRaised) routines ! LCOV_EXCL_LINE

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

#define getCosRaisedPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCosRaisedPDF_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCosRaisedPDF_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCosRaisedPDF_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCosRaisedPDF_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCosRaisedPDF_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCosRaisedPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCosRaisedPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCosRaisedPDFXDD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCosRaisedPDFXDD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCosRaisedPDFXDD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCosRaisedPDFXDD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCosRaisedPDFXDD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCosRaised@routines.inc.F90"
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
    module procedure setCosRaisedPDFXMD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCosRaisedPDFXMD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCosRaisedPDFXMD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCosRaisedPDFXMD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCosRaisedPDFXMD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XMI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCosRaisedPDFXMI_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCosRaisedPDFXMI_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCosRaisedPDFXMI_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCosRaisedPDFXMI_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCosRaisedPDFXMI_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCosRaisedPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCosRaisedCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getCosRaisedCDF_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getCosRaisedCDF_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getCosRaisedCDF_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getCosRaisedCDF_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getCosRaisedCDF_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCosRaisedCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCosRaisedCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCosRaisedCDFXDD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCosRaisedCDFXDD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCosRaisedCDFXDD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCosRaisedCDFXDD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCosRaisedCDFXDD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCosRaised@routines.inc.F90"
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
    module procedure setCosRaisedCDFXMD_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCosRaisedCDFXMD_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCosRaisedCDFXMD_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCosRaisedCDFXMD_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCosRaisedCDFXMD_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XMI_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setCosRaisedCDFXMI_RK5
        use pm_kind, only: TKC => RK5
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setCosRaisedCDFXMI_RK4
        use pm_kind, only: TKC => RK4
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setCosRaisedCDFXMI_RK3
        use pm_kind, only: TKC => RK3
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setCosRaisedCDFXMI_RK2
        use pm_kind, only: TKC => RK2
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setCosRaisedCDFXMI_RK1
        use pm_kind, only: TKC => RK1
#include "pm_distCosRaised@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XMI_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCosRaisedCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines