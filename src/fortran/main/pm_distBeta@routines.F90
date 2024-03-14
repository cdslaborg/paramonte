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
!>  This file contains procedure implementations of [pm_distBeta](@ref pm_distBeta).
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

submodule (pm_distBeta) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_mathBeta, only: getLogBeta
    use pm_mathBeta, only: getBetaInc
    use pm_mathBeta, only: setBetaInc
    use pm_distUnif, only: setUnifRand
    use pm_distGamma, only: setGammaRand

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBetaPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBetaPDF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBetaPDF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBetaPDF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBetaPDF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBetaPDF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBetaPDF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBetaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DEFAULT_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBetaLogPDFD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBetaLogPDFD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBetaLogPDFD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBetaLogPDFD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBetaLogPDFD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DEFAULT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LOGBETA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBetaLogPDFL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBetaLogPDFL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBetaLogPDFL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBetaLogPDFL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBetaLogPDFL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LOGBETA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBetaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBetaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DEFAULT_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBetaLogPDFD_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaLogPDFD_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaLogPDFD_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaLogPDFD_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaLogPDFD_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef DEFAULT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define LOGBETA_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBetaLogPDFL_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaLogPDFL_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaLogPDFL_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaLogPDFL_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaLogPDFL_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef LOGBETA_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBetaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBetaCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBetaCDF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBetaCDF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBetaCDF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBetaCDF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBetaCDF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBetaCDF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBetaCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBetaCDF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaCDF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaCDF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaCDF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaCDF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBetaCDF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBetaRand_ENABLED 1

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
    module procedure setBetaRandRNGD_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaRandRNGD_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaRandRNGD_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaRandRNGD_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaRandRNGD_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
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
    module procedure setBetaRandRNGF_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaRandRNGF_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaRandRNGF_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaRandRNGF_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaRandRNGF_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
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
    module procedure setBetaRandRNGX_D0_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaRandRNGX_D0_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaRandRNGX_D0_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaRandRNGX_D0_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaRandRNGX_D0_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
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
    module procedure setBetaRandRNGD_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaRandRNGD_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaRandRNGD_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaRandRNGD_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaRandRNGD_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
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
    module procedure setBetaRandRNGF_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaRandRNGF_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaRandRNGF_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaRandRNGF_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaRandRNGF_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
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
    module procedure setBetaRandRNGX_D1_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaRandRNGX_D1_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaRandRNGX_D1_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaRandRNGX_D1_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaRandRNGX_D1_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distBeta@routines.inc.F90"
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

#undef setBetaRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines
