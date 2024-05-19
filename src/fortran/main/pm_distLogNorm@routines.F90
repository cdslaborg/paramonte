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
!>  This file contains procedure implementations of [pm_distLogNorm](@ref pm_distLogNorm).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distLogNorm) routines ! LCOV_EXCL_LINE

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

#define getLogNormLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogNormLogPDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogNormLogPDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogNormLogPDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogNormLogPDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogNormLogPDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogNormLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setLogNormLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogNormLogPDFDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogNormLogPDFDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogNormLogPDFDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogNormLogPDFDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogNormLogPDFDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogNormLogPDFMD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogNormLogPDFMD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogNormLogPDFMD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogNormLogPDFMD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogNormLogPDFMD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogNormLogPDFDS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogNormLogPDFDS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogNormLogPDFDS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogNormLogPDFDS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogNormLogPDFDS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogNormLogPDFMS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogNormLogPDFMS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogNormLogPDFMS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogNormLogPDFMS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogNormLogPDFMS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLogNormLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogNormCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogNormCDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogNormCDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogNormCDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogNormCDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogNormCDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogNormCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setLogNormCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogNormCDFDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogNormCDFDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogNormCDFDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogNormCDFDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogNormCDFDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogNormCDFMD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogNormCDFMD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogNormCDFMD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogNormCDFMD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogNormCDFMD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setLogNormCDFMS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setLogNormCDFMS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setLogNormCDFMS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setLogNormCDFMS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setLogNormCDFMS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distLogNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setLogNormCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines