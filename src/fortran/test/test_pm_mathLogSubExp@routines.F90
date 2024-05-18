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

!>  \brief This file contains the implementations of the tests of module [pm_mathLogSubExp](@ref pm_mathLogSubExp).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_mathLogSubExp) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_getLogSubExp_CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getLogSubExp_CK5_1
        use pm_kind, only: CKC => CK5
#include "test_pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getLogSubExp_CK4_1
        use pm_kind, only: CKC => CK4
#include "test_pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getLogSubExp_CK3_1
        use pm_kind, only: CKC => CK3
#include "test_pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getLogSubExp_CK2_1
        use pm_kind, only: CKC => CK2
#include "test_pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getLogSubExp_CK1_1
        use pm_kind, only: CKC => CK1
#include "test_pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#undef test_getLogSubExp_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_getLogSubExp_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getLogSubExp_RK5_1
        use pm_kind, only: RKC => RK5
#include "test_pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLogSubExp_RK4_1
        use pm_kind, only: RKC => RK4
#include "test_pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLogSubExp_RK3_1
        use pm_kind, only: RKC => RK3
#include "test_pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogSubExp_RK2_1
        use pm_kind, only: RKC => RK2
#include "test_pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogSubExp_RK1_1
        use pm_kind, only: RKC => RK1
#include "test_pm_mathLogSubExp@routines.inc.F90"
    end procedure
#endif

#undef test_getLogSubExp_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
