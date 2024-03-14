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

!>  \brief This file contains the implementations of the tests of module [pm_cosmology](@ref pm_cosmology).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_cosmology) routines

    use pm_distUnif, only: getUnifRand
    use pm_val2str, only: getStr
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getSizeUnivNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getSizeUnivNormed_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getSizeUnivNormed_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getSizeUnivNormed_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getSizeUnivNormed_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getSizeUnivNormed_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getSizeUnivNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisLookbackNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getDisLookbackNormed_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getDisLookbackNormed_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getDisLookbackNormed_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getDisLookbackNormed_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getDisLookbackNormed_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisLookbackNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisComNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getDisComNormed_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getDisComNormed_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getDisComNormed_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getDisComNormed_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getDisComNormed_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisComNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisComTransNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getDisComTransNormed_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getDisComTransNormed_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getDisComTransNormed_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getDisComTransNormed_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getDisComTransNormed_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisComTransNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisAngNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getDisAngNormed_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getDisAngNormed_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getDisAngNormed_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getDisAngNormed_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getDisAngNormed_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisAngNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisLumNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getDisLumNormed_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getDisLumNormed_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getDisLumNormed_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getDisLumNormed_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getDisLumNormed_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisLumNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisComTransNormedWU10_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getDisComTransNormedWU10_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getDisComTransNormedWU10_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getDisComTransNormedWU10_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getDisComTransNormedWU10_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getDisComTransNormedWU10_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisComTransNormedWU10_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getHubbleParamNormedSq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getHubbleParamNormedSq_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getHubbleParamNormedSq_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getHubbleParamNormedSq_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getHubbleParamNormedSq_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getHubbleParamNormedSq_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getHubbleParamNormedSq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setVolComDiffNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setVolComDiffNormed_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setVolComDiffNormed_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setVolComDiffNormed_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setVolComDiffNormed_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setVolComDiffNormed_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVolComDiffNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getVolComDiffNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getVolComDiffNormed_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getVolComDiffNormed_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getVolComDiffNormed_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getVolComDiffNormed_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getVolComDiffNormed_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVolComDiffNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getVolComNormed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getVolComNormed_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getVolComNormed_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getVolComNormed_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getVolComNormed_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getVolComNormed_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_cosmology@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVolComNormed_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE