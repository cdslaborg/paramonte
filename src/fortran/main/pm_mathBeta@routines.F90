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
!>  This file contains procedure implementations of [pm_mathBeta](@ref pm_mathBeta).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathBeta) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_mathRoot, only: newton
    use pm_mathRoot, only: setRoot
    use pm_mathRoot, only: bisection
    use pm_mathLog1p, only: getLog1p
    use pm_math1mexp, only: get1mexp
    use pm_distNorm, only: setNormQuan
    use pm_distBeta, only: setBetaLogPDF
    use pm_quadPack, only: getQuadErr, weps
    use pm_distBeta, only: getBetaPDF
    use pm_quadPack, only: qrule => GK21
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogBeta_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogBeta_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogBeta_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogBeta_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogBeta_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogBeta_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogBeta_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBetaInc_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBetaInc_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBetaInc_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBetaInc_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBetaInc_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBetaInc_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBetaInc_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBetaInc_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBetaIncDef_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaIncDef_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaIncDef_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaIncDef_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaIncDef_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define GK21_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBetaIncGK21_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaIncGK21_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaIncGK21_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaIncGK21_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaIncGK21_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GK21_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBetaInc_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getBetaInv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getBetaInv_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getBetaInv_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getBetaInv_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getBetaInv_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getBetaInv_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getBetaInv_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setBetaInv_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setBetaInv_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setBetaInv_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setBetaInv_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setBetaInv_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setBetaInv_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathBeta@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setBetaInv_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines