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
!>  This file contains procedure implementations of [pm_distGenGamma](@ref pm_distGenGamma).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distGenGamma) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_distGamma, only: setGammaRand
    use pm_mathGamma, only: setGammaInc!Low

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGenGammaLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGenGammaLogPDFNFKDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGenGammaLogPDFNFKDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGenGammaLogPDFNFKDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGenGammaLogPDFNFKDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGenGammaLogPDFNFKDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KOD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGenGammaLogPDFNFKOD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGenGammaLogPDFNFKOD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGenGammaLogPDFNFKOD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGenGammaLogPDFNFKOD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGenGammaLogPDFNFKOD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KOD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KOS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGenGammaLogPDFNFKOS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGenGammaLogPDFNFKOS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGenGammaLogPDFNFKOS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGenGammaLogPDFNFKOS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGenGammaLogPDFNFKOS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KOS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGenGammaLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGenGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGenGammaLogPDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGenGammaLogPDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGenGammaLogPDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGenGammaLogPDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGenGammaLogPDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGenGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGenGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DDDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGenGammaLogPDFDDDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaLogPDFDDDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaLogPDFDDDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaLogPDFDDDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaLogPDFDDDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DDDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NKDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGenGammaLogPDFNKDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaLogPDFNKDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaLogPDFNKDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaLogPDFNKDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaLogPDFNKDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NKDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NKOD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGenGammaLogPDFNKOD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaLogPDFNKOD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaLogPDFNKOD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaLogPDFNKOD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaLogPDFNKOD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NKOD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NKOS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGenGammaLogPDFNKOS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaLogPDFNKOS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaLogPDFNKOS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaLogPDFNKOS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaLogPDFNKOS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NKOS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGenGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGenGammaCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGenGammaCDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGenGammaCDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGenGammaCDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGenGammaCDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGenGammaCDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGenGammaCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGenGammaCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGenGammaCDFDDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaCDFDDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaCDFDDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaCDFDDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaCDFDDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGenGammaCDFKDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaCDFKDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaCDFKDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaCDFKDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaCDFKDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KDD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KOD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGenGammaCDFKOD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaCDFKOD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaCDFKOD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaCDFKOD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaCDFKOD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KOD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KOS_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGenGammaCDFKOS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaCDFKOS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaCDFKOS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaCDFKOS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaCDFKOS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KOS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGenGammaCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGenGammaRand_ENABLED 1

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
    module procedure setGenGammaRandRNGD_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaRandRNGD_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaRandRNGD_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaRandRNGD_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaRandRNGD_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
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
    module procedure setGenGammaRandRNGF_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaRandRNGF_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaRandRNGF_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaRandRNGF_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaRandRNGF_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
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
    module procedure setGenGammaRandRNGX_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaRandRNGX_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaRandRNGX_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaRandRNGX_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaRandRNGX_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
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
    module procedure setGenGammaRandRNGD_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaRandRNGD_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaRandRNGD_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaRandRNGD_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaRandRNGD_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
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
    module procedure setGenGammaRandRNGF_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaRandRNGF_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaRandRNGF_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaRandRNGF_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaRandRNGF_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
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
    module procedure setGenGammaRandRNGX_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenGammaRandRNGX_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenGammaRandRNGX_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenGammaRandRNGX_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenGammaRandRNGX_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenGamma@routines.inc.F90"
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

#undef setGenGammaRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines