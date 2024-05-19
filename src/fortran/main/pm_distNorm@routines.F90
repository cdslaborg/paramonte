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
!>  This file contains procedure implementations of [pm_distNorm](@ref pm_distNorm).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distNorm) routines ! LCOV_EXCL_LINE

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
    use pm_ziggurat, only: getZig
    use pm_mathErf, only: setErfInv
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNormLogPDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNormLogPDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNormLogPDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNormLogPDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNormLogPDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNormLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNormLogPDFDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormLogPDFDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormLogPDFDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormLogPDFDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormLogPDFDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormLogPDFMD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormLogPDFMD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormLogPDFMD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormLogPDFMD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormLogPDFMD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormLogPDFDS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormLogPDFDS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormLogPDFDS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormLogPDFDS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormLogPDFDS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormLogPDFMS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormLogPDFMS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormLogPDFMS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormLogPDFMS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormLogPDFMS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNormLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNormCDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNormCDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNormCDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNormCDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNormCDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNormCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNormCDFDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormCDFDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormCDFDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormCDFDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormCDFDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormCDFMD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormCDFMD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormCDFMD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormCDFMD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormCDFMD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormCDFMS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormCDFMS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormCDFMS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormCDFMS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormCDFMS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNormCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNormQuan_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNormQuan_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNormQuan_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNormQuan_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNormQuan_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNormQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNormQuanDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormQuanDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormQuanDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormQuanDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormQuanDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormQuanMD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormQuanMD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormQuanMD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormQuanMD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormQuanMD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormQuanMS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormQuanMS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormQuanMS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormQuanMS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormQuanMS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNormQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getZigNorm_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getZigNorm_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getZigNorm_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getZigNorm_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getZigNorm_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getZigNorm_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getZigNorm_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getFuncNorm_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getFuncNorm_RK5
!        use pm_kind, only: RKG => RK5
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getFuncNorm_RK4
!        use pm_kind, only: RKG => RK4
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getFuncNorm_RK3
!        use pm_kind, only: RKG => RK3
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getFuncNorm_RK2
!        use pm_kind, only: RKG => RK2
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getFuncNorm_RK1
!        use pm_kind, only: RKG => RK1
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef getFuncNorm_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getGradNorm_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getGradNorm_RK5
!        use pm_kind, only: RKG => RK5
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getGradNorm_RK4
!        use pm_kind, only: RKG => RK4
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getGradNorm_RK3
!        use pm_kind, only: RKG => RK3
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getGradNorm_RK2
!        use pm_kind, only: RKG => RK2
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getGradNorm_RK1
!        use pm_kind, only: RKG => RK1
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef getGradNorm_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getFuncInvNorm_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getFuncInvNorm_RK5
!        use pm_kind, only: RKG => RK5
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getFuncInvNorm_RK4
!        use pm_kind, only: RKG => RK4
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getFuncInvNorm_RK3
!        use pm_kind, only: RKG => RK3
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getFuncInvNorm_RK2
!        use pm_kind, only: RKG => RK2
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getFuncInvNorm_RK1
!        use pm_kind, only: RKG => RK1
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef getFuncInvNorm_ENABLED
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define getZigAreaNorm_ENABLED 1
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#define RK_ENABLED 1
!
!#if RK5_ENABLED
!    module procedure getZigAreaNorm_RK5
!        use pm_kind, only: RKG => RK5
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK4_ENABLED
!    module procedure getZigAreaNorm_RK4
!        use pm_kind, only: RKG => RK4
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK3_ENABLED
!    module procedure getZigAreaNorm_RK3
!        use pm_kind, only: RKG => RK3
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK2_ENABLED
!    module procedure getZigAreaNorm_RK2
!        use pm_kind, only: RKG => RK2
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#if RK1_ENABLED
!    module procedure getZigAreaNorm_RK1
!        use pm_kind, only: RKG => RK1
!#include "pm_distNorm@routines.inc.F90"
!    end procedure
!#endif
!
!#undef RK_ENABLED
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#undef getZigAreaNorm_ENABLED
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNormRandRDMASA_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNormRandRDMASA_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNormRandRDMASA_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNormRandRDMASA_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNormRandRDMASA_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNormRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNormRandUDZD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUDZD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUDZD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUDZD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUDZD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormRandUFZD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUFZD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUFZD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUFZD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUFZD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormRandUXZD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUXZD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUXZD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUXZD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUXZD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNormRandUDZA_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUDZA_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUDZA_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUDZA_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUDZA_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormRandUFZA_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUFZA_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUFZA_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUFZA_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUFZA_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormRandUXZA_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUXZA_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUXZA_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUXZA_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUXZA_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZA_ENABLED

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

#define ZD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNormRandUDZD_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUDZD_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUDZD_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUDZD_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUDZD_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormRandUFZD_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUFZD_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUFZD_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUFZD_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUFZD_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormRandUXZD_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUXZD_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUXZD_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUXZD_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUXZD_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ZA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNormRandUDZA_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUDZA_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUDZA_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUDZA_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUDZA_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormRandUFZA_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUFZA_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUFZA_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUFZA_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUFZA_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
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
    module procedure setNormRandUXZA_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandUXZA_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandUXZA_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandUXZA_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandUXZA_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ZA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNormRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNormRandBox_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Basic_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNormRandBoxBasicDD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandBoxBasicDD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandBoxBasicDD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandBoxBasicDD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandBoxBasicDD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Basic_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Polar_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setNormRandBoxPolarDD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setNormRandBoxPolarDD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setNormRandBoxPolarDD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setNormRandBoxPolarDD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setNormRandBoxPolarDD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Polar_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNormRandBox_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormEntropy_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNormEntropy_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNormEntropy_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNormEntropy_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNormEntropy_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNormEntropy_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormEntropy_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormFisher_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNormFisher_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNormFisher_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNormFisher_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNormFisher_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNormFisher_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormFisher_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormKLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNormKLDMD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNormKLDMD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNormKLDMD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNormKLDMD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNormKLDMD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNormKLDDV_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNormKLDDV_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNormKLDDV_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNormKLDDV_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNormKLDDV_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getNormKLDMV_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getNormKLDMV_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getNormKLDMV_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getNormKLDMV_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getNormKLDMV_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormKLD_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines