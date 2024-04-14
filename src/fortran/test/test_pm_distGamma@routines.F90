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

!>  \brief This file contains the implementations of the tests of module [pm_distGamma](@ref pm_distGamma).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distGamma) routines

    use pm_mathCompare, only: isClose
    use pm_distUnif, only: setUnifRand
    use pm_arraySpace, only: setLinSpace
    use pm_distGenGamma, only: getGenGammaLogPDF
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getGammaLogPDF_RK5_1
        use pm_kind, only: TKC => RK5
#include "test_pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getGammaLogPDF_RK4_1
        use pm_kind, only: TKC => RK4
#include "test_pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getGammaLogPDF_RK3_1
        use pm_kind, only: TKC => RK3
#include "test_pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getGammaLogPDF_RK2_1
        use pm_kind, only: TKC => RK2
#include "test_pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getGammaLogPDF_RK1_1
        use pm_kind, only: TKC => RK1
#include "test_pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setGammaLogPDF_RK5_1
        use pm_kind, only: TKC => RK5
#include "test_pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setGammaLogPDF_RK4_1
        use pm_kind, only: TKC => RK4
#include "test_pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setGammaLogPDF_RK3_1
        use pm_kind, only: TKC => RK3
#include "test_pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setGammaLogPDF_RK2_1
        use pm_kind, only: TKC => RK2
#include "test_pm_distGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setGammaLogPDF_RK1_1
        use pm_kind, only: TKC => RK1
#include "test_pm_distGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE