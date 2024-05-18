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

!>  \brief This file contains the implementations of the tests of module [pm_distPareto](@ref pm_distPareto).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distPareto) routines

    use pm_option, only: getOption
    use pm_arraySort, only: setSorted
    use pm_quadPack, only: isFailedQuad
    use pm_arraySpace, only: getLinSpace
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_distPareto, only: getParetoLogPDF

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getParetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getParetoLogPDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getParetoLogPDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getParetoLogPDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getParetoLogPDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getParetoLogPDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getParetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setParetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setParetoLogPDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setParetoLogPDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setParetoLogPDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setParetoLogPDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setParetoLogPDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setParetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getParetoLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getParetoLogCDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getParetoLogCDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getParetoLogCDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getParetoLogCDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getParetoLogCDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getParetoLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setParetoLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setParetoLogCDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setParetoLogCDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setParetoLogCDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setParetoLogCDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setParetoLogCDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setParetoLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getParetoLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getParetoLogQuan_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getParetoLogQuan_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getParetoLogQuan_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getParetoLogQuan_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getParetoLogQuan_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getParetoLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setParetoLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setParetoLogQuan_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setParetoLogQuan_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setParetoLogQuan_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setParetoLogQuan_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setParetoLogQuan_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setParetoLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getParetoLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getParetoLogRand_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getParetoLogRand_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getParetoLogRand_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getParetoLogRand_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getParetoLogRand_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getParetoLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setParetoLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setParetoLogRand_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setParetoLogRand_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setParetoLogRand_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setParetoLogRand_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setParetoLogRand_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPareto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setParetoLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE