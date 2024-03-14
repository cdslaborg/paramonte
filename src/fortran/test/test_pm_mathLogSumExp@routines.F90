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

!>  \brief This file contains the implementations of the tests of module [pm_mathLogSumExp](@ref pm_mathLogSumExp).
!>
!>  \author
!>  \AmirShahmoradi

submodule (test_pm_mathLogSumExp) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_getLogSumExp_CK_ENABLED 1

#if CK3_ENABLED
    module procedure test_getLogSumExp_CK3_1
        use pm_kind, only: IK, CK => CK3
#include "test_pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getLogSumExp_CK2_1
        use pm_kind, only: IK, CK => CK2
#include "test_pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getLogSumExp_CK1_1
        use pm_kind, only: IK, CK => CK1
#include "test_pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#undef test_getLogSumExp_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_getLogSumExp_RK_ENABLED 1

#if RK3_ENABLED
    module procedure test_getLogSumExp_RK3_1
        use pm_kind, only: IK, RK => RK3
#include "test_pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogSumExp_RK2_1
        use pm_kind, only: IK, RK => RK2
#include "test_pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogSumExp_RK1_1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_mathLogSumExp@routines.inc.F90"
    end procedure
#endif

#undef test_getLogSumExp_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
