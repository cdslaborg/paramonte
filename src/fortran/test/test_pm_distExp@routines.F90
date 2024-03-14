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

!>  \brief This file contains the implementations of the tests of module [pm_distExp](@ref pm_distExp).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distExp) routines

    use pm_distUnif, only: getUnifRand
    use pm_arrayChoice, only: getChoice
    use pm_mathCompare, only: isClose
    use pm_kind, only: SK, IK
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getExpLogPDF_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getExpLogPDF_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getExpLogPDF_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getExpLogPDF_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getExpLogPDF_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setExpLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setExpLogPDF_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setExpLogPDF_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setExpLogPDF_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setExpLogPDF_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setExpLogPDF_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setExpLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getExpCDF_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getExpCDF_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getExpCDF_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getExpCDF_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getExpCDF_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setExpCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setExpCDF_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setExpCDF_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setExpCDF_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setExpCDF_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setExpCDF_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setExpCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getExpRand_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getExpRand_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getExpRand_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getExpRand_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getExpRand_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setExpRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setExpRand_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setExpRand_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setExpRand_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setExpRand_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setExpRand_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setExpRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE