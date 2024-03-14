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

!>  \brief This file contains the implementations of the tests of module [pm_mathLogAddExp](@ref pm_mathLogAddExp).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_mathLogAddExp) routines

    use pm_mathMinMax, only: setMinMax
    use pm_distUnif, only: setUnifRand
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogAddExp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getLogAddExp_CK5
        use pm_kind, only: CKC => CK5
#include "test_pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getLogAddExp_CK4
        use pm_kind, only: CKC => CK4
#include "test_pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getLogAddExp_CK3
        use pm_kind, only: CKC => CK3
#include "test_pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getLogAddExp_CK2
        use pm_kind, only: CKC => CK2
#include "test_pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getLogAddExp_CK1
        use pm_kind, only: CKC => CK1
#include "test_pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getLogAddExp_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getLogAddExp_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getLogAddExp_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getLogAddExp_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getLogAddExp_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_mathLogAddExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogAddExp_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE