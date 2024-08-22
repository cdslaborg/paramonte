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
!>  This file contains procedure implementations of [pm_distGenExpGamma](@ref pm_distGenExpGamma).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distGenExpGamma) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_mathGamma, only: setGammaInc!Low
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGenExpGammaLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGenExpGammaLogPDFNFKD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGenExpGammaLogPDFNFKD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGenExpGammaLogPDFNFKD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGenExpGammaLogPDFNFKD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGenExpGammaLogPDFNFKD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define KO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGenExpGammaLogPDFNFKO_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGenExpGammaLogPDFNFKO_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGenExpGammaLogPDFNFKO_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGenExpGammaLogPDFNFKO_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGenExpGammaLogPDFNFKO_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGenExpGammaLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGenExpGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGenExpGammaLogPDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGenExpGammaLogPDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGenExpGammaLogPDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGenExpGammaLogPDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGenExpGammaLogPDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGenExpGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGenExpGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DDDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGenExpGammaLogPDFDDDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenExpGammaLogPDFDDDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenExpGammaLogPDFDDDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenExpGammaLogPDFDDDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenExpGammaLogPDFDDDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
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
    module procedure setGenExpGammaLogPDFNKDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenExpGammaLogPDFNKDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenExpGammaLogPDFNKDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenExpGammaLogPDFNKDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenExpGammaLogPDFNKDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
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
    module procedure setGenExpGammaLogPDFNKOD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenExpGammaLogPDFNKOD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenExpGammaLogPDFNKOD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenExpGammaLogPDFNKOD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenExpGammaLogPDFNKOD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
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
    module procedure setGenExpGammaLogPDFNKOS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenExpGammaLogPDFNKOS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenExpGammaLogPDFNKOS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenExpGammaLogPDFNKOS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenExpGammaLogPDFNKOS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NKOS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGenExpGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGenExpGammaCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGenExpGammaCDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGenExpGammaCDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGenExpGammaCDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGenExpGammaCDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGenExpGammaCDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGenExpGammaCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGenExpGammaCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGenExpGammaCDFDDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenExpGammaCDFDDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenExpGammaCDFDDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenExpGammaCDFDDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenExpGammaCDFDDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
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
    module procedure setGenExpGammaCDFXKDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenExpGammaCDFXKDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenExpGammaCDFXKDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenExpGammaCDFXKDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenExpGammaCDFXKDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
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
    module procedure setGenExpGammaCDFKOD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenExpGammaCDFKOD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenExpGammaCDFKOD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenExpGammaCDFKOD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenExpGammaCDFKOD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
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
    module procedure setGenExpGammaCDFKOS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGenExpGammaCDFKOS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGenExpGammaCDFKOS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGenExpGammaCDFKOS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGenExpGammaCDFKOS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KOS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGenExpGammaCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines