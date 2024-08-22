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
!>  This file contains procedure implementations of [pm_mathGammaAM](@ref pm_mathGammaAM).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, July 22, 2024, 11:45 AM, NASA Goddard Space Flight Center, Washington, D.C.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathGammaAM) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
    use pm_option, only: getOption
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
    use pm_val2str, only: getStr
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGammaIncLowAM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaIncLowAM_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaIncLowAM_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaIncLowAM_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaIncLowAM_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaIncLowAM_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaIncLowAM_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGammaIncUppAM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaIncUppAM_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaIncUppAM_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaIncUppAM_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaIncUppAM_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaIncUppAM_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaIncUppAM_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaIncLowAM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaIncLowAM_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaIncLowAM_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaIncLowAM_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaIncLowAM_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaIncLowAM_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaIncLowAM_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaIncUppAM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaIncUppAM_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaIncUppAM_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaIncUppAM_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaIncUppAM_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaIncUppAM_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaIncUppAM_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaIncAM_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaIncAMDef_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaIncAMDef_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaIncAMDef_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaIncAMDef_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaIncAMDef_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define NXIK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaIncAMNXIK_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaIncAMNXIK_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaIncAMNXIK_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaIncAMNXIK_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaIncAMNXIK_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaAM@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef NXIK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaIncAM_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines