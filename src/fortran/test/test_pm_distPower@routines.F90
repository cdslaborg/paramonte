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

!>  \brief This file contains the implementations of the tests of module [pm_distPower](@ref pm_distPower).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distPower) routines

    use pm_option, only: getOption
    use pm_arraySort, only: setSorted
    use pm_arraySpace, only: getLinSpace
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowerLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPowerLogPDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPowerLogPDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPowerLogPDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPowerLogPDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPowerLogPDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowerLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowerLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPowerLogPDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPowerLogPDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPowerLogPDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPowerLogPDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPowerLogPDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowerLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowerLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPowerLogCDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPowerLogCDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPowerLogCDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPowerLogCDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPowerLogCDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowerLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowerLogCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPowerLogCDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPowerLogCDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPowerLogCDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPowerLogCDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPowerLogCDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowerLogCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowerLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPowerLogQuan_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPowerLogQuan_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPowerLogQuan_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPowerLogQuan_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPowerLogQuan_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowerLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowerLogQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPowerLogQuan_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPowerLogQuan_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPowerLogQuan_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPowerLogQuan_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPowerLogQuan_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowerLogQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowerLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPowerLogRand_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPowerLogRand_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPowerLogRand_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPowerLogRand_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPowerLogRand_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowerLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowerLogRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPowerLogRand_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPowerLogRand_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPowerLogRand_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPowerLogRand_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPowerLogRand_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPower@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowerLogRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE