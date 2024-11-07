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

!>  \brief This file contains the implementations of the tests of module [pm_matrixDet](@ref pm_matrixDet).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_matrixDet) routines

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_getDet_ENABLED 1

#if RK3_ENABLED
    module procedure test_getDet_RK3
        use pm_kind, only: IK, RK => RK3
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getDet_RK2
        use pm_kind, only: IK, RK => RK2
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getDet_RK1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef test_getDet_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_setDet_ENABLED 1

#if RK3_ENABLED
    module procedure test_setDet_RK3
        use pm_kind, only: IK, RK => RK3
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setDet_RK2
        use pm_kind, only: IK, RK => RK2
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setDet_RK1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef test_setDet_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_getMatDetSqrtLog_ENABLED 1

#if RK3_ENABLED
    module procedure test_getMatDetSqrtLog_RK3
        use pm_kind, only: IK, RK => RK3
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getMatDetSqrtLog_RK2
        use pm_kind, only: IK, RK => RK2
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getMatDetSqrtLog_RK1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef test_getMatDetSqrtLog_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_setMatDetSqrtLog_ENABLED 1

#if RK3_ENABLED
    module procedure test_setMatDetSqrtLog_RK3
        use pm_kind, only: IK, RK => RK3
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setMatDetSqrtLog_RK2
        use pm_kind, only: IK, RK => RK2
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setMatDetSqrtLog_RK1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef test_setMatDetSqrtLog_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_getMatDetSqrt_ENABLED 1

#if RK3_ENABLED
    module procedure test_getMatDetSqrt_RK3
        use pm_kind, only: IK, RK => RK3
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getMatDetSqrt_RK2
        use pm_kind, only: IK, RK => RK2
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getMatDetSqrt_RK1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef test_getMatDetSqrt_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define test_setMatDetSqrt_ENABLED 1

#if RK3_ENABLED
    module procedure test_setMatDetSqrt_RK3
        use pm_kind, only: IK, RK => RK3
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setMatDetSqrt_RK2
        use pm_kind, only: IK, RK => RK2
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setMatDetSqrt_RK1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_matrixDet@routines.inc.F90"
    end procedure
#endif

#undef test_setMatDetSqrt_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE