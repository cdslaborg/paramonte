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
!>  This file contains procedure implementations of [pm_mathGammaGil](@ref pm_mathGammaGil).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, July 22, 2024, 11:45 AM, NASA Goddard Space Flight Center, Washington, D.C.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathGammaGil) routines ! LCOV_EXCL_LINE

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

#define getGammaIncLowGil_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaIncLowGil_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaIncLowGil_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaIncLowGil_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaIncLowGil_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaIncLowGil_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaIncLowGil_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGammaIncUppGil_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaIncUppGil_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaIncUppGil_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaIncUppGil_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaIncUppGil_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaIncUppGil_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaIncUppGil_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaIncGil_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaIncGilDef_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaIncGilDef_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaIncGilDef_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaIncGilDef_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaIncGilDef_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaGil@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaIncGil_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines