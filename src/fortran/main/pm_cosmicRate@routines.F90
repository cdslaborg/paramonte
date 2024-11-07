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
!>  This file contains procedure implementations of [pm_cosmicRate](@ref pm_cosmicRate).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Wednesday 5:43 PM, December 25, 2013, Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_cosmicRate) routines

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogRateDensity_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define H06_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogRateDensityH06_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogRateDensityH06_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogRateDensityH06_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogRateDensityH06_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogRateDensityH06_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef H06_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define L08_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogRateDensityL08_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogRateDensityL08_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogRateDensityL08_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogRateDensityL08_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogRateDensityL08_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef L08_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define B10_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogRateDensityB10_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogRateDensityB10_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogRateDensityB10_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogRateDensityB10_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogRateDensityB10_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef B10_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define M14_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogRateDensityM14_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogRateDensityM14_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogRateDensityM14_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogRateDensityM14_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogRateDensityM14_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef M14_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define P15_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogRateDensityP15_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogRateDensityP15_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogRateDensityP15_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogRateDensityP15_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogRateDensityP15_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef P15_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define M17_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogRateDensityM17_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogRateDensityM17_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogRateDensityM17_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogRateDensityM17_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogRateDensityM17_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef M17_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define F18_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogRateDensityF18_D0_RK5
        use pm_kind, only: RKG => RK5
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogRateDensityF18_D0_RK4
        use pm_kind, only: RKG => RK4
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogRateDensityF18_D0_RK3
        use pm_kind, only: RKG => RK3
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogRateDensityF18_D0_RK2
        use pm_kind, only: RKG => RK2
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogRateDensityF18_D0_RK1
        use pm_kind, only: RKG => RK1
#include "pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef F18_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogRateDensity_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines