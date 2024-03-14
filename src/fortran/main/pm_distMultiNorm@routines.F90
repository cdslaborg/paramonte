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
!>  This file contains procedure implementations of [pm_distMultiNorm](@ref pm_distMultiNorm).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distMultiNorm) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif
    
    use pm_distNorm, only: setNormRand
    use pm_distanceMahal, only: getMahalSq
    use pm_matrixCopy, only: setMatCopy, rdpack
    use pm_matrixDet, only: getMatDetSqrtLog, setMatDetSqrtLog, nothing
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMultiNormLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define I_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFNFI_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFNFI_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFNFI_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFNFI_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFNFI_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef I_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IF_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFNFIF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFNFIF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFNFIF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFNFIF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFNFIF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef IF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ND_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFNFND_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFNFND_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFNFND_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFNFND_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFNFND_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef ND_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMultiNormLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMultiNormLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DDD_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFDDD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFDDD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFDDD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFDDD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFDDD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef DDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MDD_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFMDD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFMDD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFMDD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFMDD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFMDD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef MDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DID_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFDID_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFDID_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFDID_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFDID_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFDID_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef DID_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MID_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFMID_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFMID_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFMID_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFMID_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFMID_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef MID_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DDN_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFDDN_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFDDN_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFDDN_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFDDN_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFDDN_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef DDN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MDN_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFMDN_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFMDN_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFMDN_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFMDN_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFMDN_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef MDN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DIN_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFDIN_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFDIN_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFDIN_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFDIN_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFDIN_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef DIN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MIN_ENABLED 1

#if RK5_ENABLED
    module procedure getMultiNormLogPDFMIN_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMultiNormLogPDFMIN_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMultiNormLogPDFMIN_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMultiNormLogPDFMIN_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMultiNormLogPDFMIN_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef MIN_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMultiNormLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMNR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGD_DM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGD_DM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGD_DM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGD_DM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGD_DM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGD_DM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGD_DM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGD_DM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGD_DM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGD_DM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGD_AM_DC_XXX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGD_AM_DC_XXX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGD_AM_DC_XXX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGD_AM_DC_XXX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGD_AM_DC_XXX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGD_AM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGD_AM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGD_AM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGD_AM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGD_AM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGD_AM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGD_AM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGD_AM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGD_AM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGD_AM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGF_DM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGF_DM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGF_DM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGF_DM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGF_DM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGF_DM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGF_DM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGF_DM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGF_DM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGF_DM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGF_AM_DC_XXX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGF_AM_DC_XXX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGF_AM_DC_XXX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGF_AM_DC_XXX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGF_AM_DC_XXX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGF_AM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGF_AM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGF_AM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGF_AM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGF_AM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGF_AM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGF_AM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGF_AM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGF_AM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGF_AM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGX_DM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGX_DM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGX_DM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGX_DM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGX_DM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGX_DM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGX_DM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGX_DM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGX_DM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGX_DM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGX_AM_DC_XXX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGX_AM_DC_XXX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGX_AM_DC_XXX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGX_AM_DC_XXX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGX_AM_DC_XXX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGX_AM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGX_AM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGX_AM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGX_AM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGX_AM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGX_AM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGX_AM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGX_AM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGX_AM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGX_AM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGD_DM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGD_DM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGD_DM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGD_DM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGD_DM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGD_DM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGD_DM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGD_DM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGD_DM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGD_DM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGD_AM_DC_XXX_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGD_AM_DC_XXX_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGD_AM_DC_XXX_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGD_AM_DC_XXX_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGD_AM_DC_XXX_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGD_AM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGD_AM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGD_AM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGD_AM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGD_AM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGD_AM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGD_AM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGD_AM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGD_AM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGD_AM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGF_DM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGF_DM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGF_DM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGF_DM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGF_DM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGF_DM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGF_DM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGF_DM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGF_DM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGF_DM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGF_AM_DC_XXX_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGF_AM_DC_XXX_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGF_AM_DC_XXX_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGF_AM_DC_XXX_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGF_AM_DC_XXX_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGF_AM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGF_AM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGF_AM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGF_AM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGF_AM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGF_AM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGF_AM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGF_AM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGF_AM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGF_AM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGX_DM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGX_DM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGX_DM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGX_DM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGX_DM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGX_DM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGX_DM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGX_DM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGX_DM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGX_DM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGX_AM_DC_XXX_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGX_AM_DC_XXX_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGX_AM_DC_XXX_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGX_AM_DC_XXX_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGX_AM_DC_XXX_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGX_AM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGX_AM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGX_AM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGX_AM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGX_AM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getMNR_RNGX_AM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getMNR_RNGX_AM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getMNR_RNGX_AM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getMNR_RNGX_AM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getMNR_RNGX_AM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMNR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMNR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_DM_DC_XXX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_DM_DC_XXX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_DM_DC_XXX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_DM_DC_XXX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_DM_DC_XXX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_DM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_DM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_DM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_DM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_DM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_DM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_DM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_DM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_DM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_DM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_AM_DC_XXX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_AM_DC_XXX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_AM_DC_XXX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_AM_DC_XXX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_AM_DC_XXX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_AM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_AM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_AM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_AM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_AM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_AM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_AM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_AM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_AM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_AM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_DM_DC_XXX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_DM_DC_XXX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_DM_DC_XXX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_DM_DC_XXX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_DM_DC_XXX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_DM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_DM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_DM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_DM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_DM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_DM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_DM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_DM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_DM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_DM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_AM_DC_XXX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_AM_DC_XXX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_AM_DC_XXX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_AM_DC_XXX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_AM_DC_XXX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_AM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_AM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_AM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_AM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_AM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_AM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_AM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_AM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_AM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_AM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_DM_DC_XXX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_DM_DC_XXX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_DM_DC_XXX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_DM_DC_XXX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_DM_DC_XXX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_DM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_DM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_DM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_DM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_DM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_DM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_DM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_DM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_DM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_DM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_AM_DC_XXX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_AM_DC_XXX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_AM_DC_XXX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_AM_DC_XXX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_AM_DC_XXX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_AM_AC_UXD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_AM_AC_UXD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_AM_AC_UXD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_AM_AC_UXD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_AM_AC_UXD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_AM_AC_XLD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_AM_AC_XLD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_AM_AC_XLD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_AM_AC_XLD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_AM_AC_XLD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_DM_DC_XXX_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_DM_DC_XXX_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_DM_DC_XXX_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_DM_DC_XXX_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_DM_DC_XXX_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_DM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_DM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_DM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_DM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_DM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_DM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_DM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_DM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_DM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_DM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_AM_DC_XXX_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_AM_DC_XXX_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_AM_DC_XXX_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_AM_DC_XXX_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_AM_DC_XXX_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_AM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_AM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_AM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_AM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_AM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGD_AM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGD_AM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGD_AM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGD_AM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGD_AM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_DM_DC_XXX_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_DM_DC_XXX_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_DM_DC_XXX_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_DM_DC_XXX_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_DM_DC_XXX_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_DM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_DM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_DM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_DM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_DM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_DM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_DM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_DM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_DM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_DM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_AM_DC_XXX_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_AM_DC_XXX_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_AM_DC_XXX_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_AM_DC_XXX_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_AM_DC_XXX_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_AM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_AM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_AM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_AM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_AM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGF_AM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGF_AM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGF_AM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGF_AM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGF_AM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RNGX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_DM_DC_XXX_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_DM_DC_XXX_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_DM_DC_XXX_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_DM_DC_XXX_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_DM_DC_XXX_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_DM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_DM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_DM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_DM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_DM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_DM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_DM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_DM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_DM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_DM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_AM_DC_XXX_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_AM_DC_XXX_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_AM_DC_XXX_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_AM_DC_XXX_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_AM_DC_XXX_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define AC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define UXD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_AM_AC_UXD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_AM_AC_UXD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_AM_AC_UXD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_AM_AC_UXD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_AM_AC_UXD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef UXD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setMNR_RNGX_AM_AC_XLD_D2_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setMNR_RNGX_AM_AC_XLD_D2_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setMNR_RNGX_AM_AC_XLD_D2_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setMNR_RNGX_AM_AC_XLD_D2_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setMNR_RNGX_AM_AC_XLD_D2_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XLD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef AM_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RNGX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMNR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines