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
!>  This file contains procedure implementations of [pm_distExpGamma](@ref pm_distExpGamma).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distExpGamma) routines ! LCOV_EXCL_LINE

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

#define getExpGammaLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getExpGammaLogPDFNF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExpGammaLogPDFNF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExpGammaLogPDFNF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExpGammaLogPDFNF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExpGammaLogPDFNF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpGammaLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getExpGammaLogPDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExpGammaLogPDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExpGammaLogPDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExpGammaLogPDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExpGammaLogPDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setExpGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DDD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExpGammaLogPDFDDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpGammaLogPDFDDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpGammaLogPDFDDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpGammaLogPDFDDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpGammaLogPDFDDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distExpGamma@routines.inc.F90"
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
    module procedure setExpGammaLogPDFNKD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpGammaLogPDFNKD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpGammaLogPDFNKD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpGammaLogPDFNKD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpGammaLogPDFNKD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distExpGamma@routines.inc.F90"
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
    module procedure setExpGammaLogPDFNKS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpGammaLogPDFNKS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpGammaLogPDFNKS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpGammaLogPDFNKS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpGammaLogPDFNKS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NKS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setExpGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpGammaCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getExpGammaCDF_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getExpGammaCDF_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getExpGammaCDF_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getExpGammaCDF_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getExpGammaCDF_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpGammaCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setExpGammaCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setExpGammaCDFDD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpGammaCDFDD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpGammaCDFDD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpGammaCDFDD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpGammaCDFDD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distExpGamma@routines.inc.F90"
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
    module procedure setExpGammaCDFKD_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpGammaCDFKD_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpGammaCDFKD_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpGammaCDFKD_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpGammaCDFKD_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distExpGamma@routines.inc.F90"
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
    module procedure setExpGammaCDFKS_RK5
        use pm_kind, only: RKG => RK5
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setExpGammaCDFKS_RK4
        use pm_kind, only: RKG => RK4
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setExpGammaCDFKS_RK3
        use pm_kind, only: RKG => RK3
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setExpGammaCDFKS_RK2
        use pm_kind, only: RKG => RK2
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setExpGammaCDFKS_RK1
        use pm_kind, only: RKG => RK1
#include "pm_distExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef KS_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setExpGammaCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines