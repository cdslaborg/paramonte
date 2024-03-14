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

!>  \brief This file contains the implementations of the tests of module [pm_distGenExpGamma](@ref pm_distGenExpGamma).
!>
!>  \author
!>  \AmirShahmoradi

submodule (test_pm_distGenExpGamma) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGenExpGammaLogPDFNF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getGenExpGammaLogPDFNF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getGenExpGammaLogPDFNF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getGenExpGammaLogPDFNF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getGenExpGammaLogPDFNF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getGenExpGammaLogPDFNF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGenExpGammaLogPDFNF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getGenExpGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getGenExpGammaLogPDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getGenExpGammaLogPDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getGenExpGammaLogPDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getGenExpGammaLogPDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getGenExpGammaLogPDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getGenExpGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setGenExpGammaLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setGenExpGammaLogPDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setGenExpGammaLogPDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setGenExpGammaLogPDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setGenExpGammaLogPDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setGenExpGammaLogPDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distGenExpGamma@routines.inc.F90"
    end procedure
#endif

#undef setGenExpGammaLogPDF_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setGenExpGammaLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE