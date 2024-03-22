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

!>  \brief This file contains the implementations of the tests of module [pm_distPiwiPoweto](@ref pm_distPiwiPoweto).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distPiwiPoweto) routines

    use pm_arraySort, only: setSorted
    use pm_quadPack, only: isFailedQuad
    use pm_arraySpace, only: getLinSpace
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPiwiPowetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPiwiPowetoLogPDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPiwiPowetoLogPDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPiwiPowetoLogPDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPiwiPowetoLogPDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPiwiPowetoLogPDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPiwiPowetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPiwiPowetoLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPiwiPowetoLogPDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPiwiPowetoLogPDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPiwiPowetoLogPDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPiwiPowetoLogPDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPiwiPowetoLogPDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPiwiPowetoLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getPiwiPowetoCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getPiwiPowetoCDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getPiwiPowetoCDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getPiwiPowetoCDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getPiwiPowetoCDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getPiwiPowetoCDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getPiwiPowetoCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setPiwiPowetoCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setPiwiPowetoCDF_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setPiwiPowetoCDF_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setPiwiPowetoCDF_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setPiwiPowetoCDF_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setPiwiPowetoCDF_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_distPiwiPoweto@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setPiwiPowetoCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
