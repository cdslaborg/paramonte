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
!>  This file contains procedure implementations of [pm_mathGamma](@ref pm_mathGamma).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathGamma) routines ! LCOV_EXCL_LINE

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

    use pm_mathGammaGil, only: setGammaIncGil

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGammaIncLow_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaIncLow_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaIncLow_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaIncLow_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaIncLow_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaIncLow_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaIncLow_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGammaIncUpp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaIncUpp_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaIncUpp_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaIncUpp_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaIncUpp_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaIncUpp_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaIncUpp_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaInc_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaInc_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaInc_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaInc_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaInc_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaInc_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaInc_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines