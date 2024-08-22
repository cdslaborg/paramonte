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
!>  This file contains procedure implementations of [pm_mathGammaNR](@ref pm_mathGammaNR).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathGammaNR) routines ! LCOV_EXCL_LINE

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

#define getGammaIncLowNR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaIncLowNR_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaIncLowNR_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaIncLowNR_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaIncLowNR_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaIncLowNR_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaIncLowNR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGammaIncUppNR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getGammaIncUppNR_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getGammaIncUppNR_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getGammaIncUppNR_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getGammaIncUppNR_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getGammaIncUppNR_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaIncUppNR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaIncLowNR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaIncLowNR_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaIncLowNR_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaIncLowNR_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaIncLowNR_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaIncLowNR_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaIncLowNR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaIncUppNR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaIncUppNR_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaIncUppNR_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaIncUppNR_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaIncUppNR_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaIncUppNR_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaIncUppNR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaIncLowSeriesNR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaIncLowSeriesNR_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaIncLowSeriesNR_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaIncLowSeriesNR_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaIncLowSeriesNR_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaIncLowSeriesNR_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaIncLowSeriesNR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaIncUppContFracNR_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setGammaIncUppContFracNR_RK5
        use pm_kind, only: RKG => RK5
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setGammaIncUppContFracNR_RK4
        use pm_kind, only: RKG => RK4
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setGammaIncUppContFracNR_RK3
        use pm_kind, only: RKG => RK3
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setGammaIncUppContFracNR_RK2
        use pm_kind, only: RKG => RK2
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setGammaIncUppContFracNR_RK1
        use pm_kind, only: RKG => RK1
#include "pm_mathGammaNR@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaIncUppContFracNR_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines