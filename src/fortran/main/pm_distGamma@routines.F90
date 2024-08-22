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
!>  This file contains procedure implementations of [pm_distGamma](@ref pm_distGamma).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distGamma) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_mathGamma, only: setGammaInc
    use pm_distUnif, only: setUnifRand
    use pm_distNorm, only: setNormRand
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGammaLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaLogPDFNFKD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaLogPDFNFKD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaLogPDFNFKD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaLogPDFNFKD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaLogPDFNFKD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaLogPDFNFKS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaLogPDFNFKS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaLogPDFNFKS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaLogPDFNFKS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaLogPDFNFKS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaLogPDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaLogPDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaLogPDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaLogPDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaLogPDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaLogPDFDDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaLogPDFDDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaLogPDFDDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaLogPDFDDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaLogPDFDDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NKD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaLogPDFNKD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaLogPDFNKD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaLogPDFNKD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaLogPDFNKD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaLogPDFNKD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NKD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NKS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaLogPDFNKS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaLogPDFNKS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaLogPDFNKS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaLogPDFNKS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaLogPDFNKS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NKS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGammaCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaCDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaCDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaCDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaCDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaCDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaCDFDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaCDFDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaCDFDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaCDFDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaCDFDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaCDFKD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaCDFKD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaCDFKD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaCDFKD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaCDFKD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaCDFKS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaCDFKS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaCDFKS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaCDFKS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaCDFKS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KR_ENABLED 1

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
    module procedure setGammaRandRNGD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaRandRNGD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaRandRNGD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaRandRNGD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaRandRNGD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
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
    module procedure setGammaRandRNGF_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaRandRNGF_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaRandRNGF_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaRandRNGF_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaRandRNGF_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
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
    module procedure setGammaRandRNGX_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaRandRNGX_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaRandRNGX_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaRandRNGX_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaRandRNGX_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
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
    module procedure setGammaRandRNGD_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaRandRNGD_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaRandRNGD_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaRandRNGD_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaRandRNGD_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
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
    module procedure setGammaRandRNGF_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaRandRNGF_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaRandRNGF_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaRandRNGF_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaRandRNGF_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
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
    module procedure setGammaRandRNGX_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaRandRNGX_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaRandRNGX_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaRandRNGX_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaRandRNGX_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGamma@routines.inc.F90"
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

#undef KR_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines