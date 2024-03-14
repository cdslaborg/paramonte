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

!>  \brief This file contains the implementations of the tests of module [pm_cosmicRate](@ref pm_cosmicRate).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_cosmicRate) routines

    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
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
    module procedure test_getLogRateDensityH06_D0_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLogRateDensityH06_D0_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLogRateDensityH06_D0_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogRateDensityH06_D0_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogRateDensityH06_D0_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmicRate@routines.inc.F90"
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
    module procedure test_getLogRateDensityL08_D0_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLogRateDensityL08_D0_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLogRateDensityL08_D0_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogRateDensityL08_D0_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogRateDensityL08_D0_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmicRate@routines.inc.F90"
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
    module procedure test_getLogRateDensityB10_D0_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLogRateDensityB10_D0_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLogRateDensityB10_D0_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogRateDensityB10_D0_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogRateDensityB10_D0_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmicRate@routines.inc.F90"
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
    module procedure test_getLogRateDensityM14_D0_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLogRateDensityM14_D0_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLogRateDensityM14_D0_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogRateDensityM14_D0_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogRateDensityM14_D0_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmicRate@routines.inc.F90"
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
    module procedure test_getLogRateDensityP15_D0_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLogRateDensityP15_D0_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLogRateDensityP15_D0_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogRateDensityP15_D0_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogRateDensityP15_D0_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmicRate@routines.inc.F90"
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
    module procedure test_getLogRateDensityM17_D0_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLogRateDensityM17_D0_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLogRateDensityM17_D0_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogRateDensityM17_D0_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogRateDensityM17_D0_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmicRate@routines.inc.F90"
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
    module procedure test_getLogRateDensityF18_D0_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLogRateDensityF18_D0_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLogRateDensityF18_D0_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogRateDensityF18_D0_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmicRate@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogRateDensityF18_D0_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmicRate@routines.inc.F90"
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

end submodule routines ! LCOV_EXCL_LINE