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

!>  \brief This file contains the implementations of the tests of module [pm_distPoweto](@ref pm_distPoweto).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distPoweto) routines

    use pm_arraySort, only: setSorted
    use pm_quadPack, only: isFailedQuad
    use pm_arraySpace, only: getLinSpace
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPowetoLogPDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPowetoLogPDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPowetoLogPDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPowetoLogPDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPowetoLogPDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPowetoLogPDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPowetoLogPDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPowetoLogPDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPowetoLogPDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPowetoLogPDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPowetoCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPowetoCDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPowetoCDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPowetoCDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPowetoCDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPowetoCDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPowetoCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPowetoCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPowetoCDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPowetoCDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPowetoCDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPowetoCDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPowetoCDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPowetoCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
