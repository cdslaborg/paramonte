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

!>  \brief
!>  This module contains implementations of the tests of the module [pm_mathFactorial](@ref pm_mathFactorial).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_mathFactorial) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getFactorial_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getFactorial_IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getFactorial_IK5_1
        use pm_kind, only: IKC => IK5
#include "test_pm_mathFactorial@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getFactorial_IK4_1
        use pm_kind, only: IKC => IK4
#include "test_pm_mathFactorial@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getFactorial_IK3_1
        use pm_kind, only: IKC => IK3
#include "test_pm_mathFactorial@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getFactorial_IK2_1
        use pm_kind, only: IKC => IK2
#include "test_pm_mathFactorial@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getFactorial_IK1_1
        use pm_kind, only: IKC => IK1
#include "test_pm_mathFactorial@routines.inc.F90"
    end procedure
#endif

#undef getFactorial_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getFactorial_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogFactorial_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLogFactorial_RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getLogFactorial_RK5_1
        use pm_kind, only: IK, RK => RK5
#include "test_pm_mathFactorial@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getLogFactorial_RK4_1
        use pm_kind, only: IK, RK => RK4
#include "test_pm_mathFactorial@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getLogFactorial_RK3_1
        use pm_kind, only: IK, RK => RK3
#include "test_pm_mathFactorial@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getLogFactorial_RK2_1
        use pm_kind, only: IK, RK => RK2
#include "test_pm_mathFactorial@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getLogFactorial_RK1_1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_mathFactorial@routines.inc.F90"
    end procedure
#endif

#undef getLogFactorial_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLogFactorial_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE
