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

!>  \brief This file contains the implementations of the tests of module [pm_matrixInitDia](@ref pm_matrixInitDia).
!>
!>  \author
!>  \AmirShahmoradi

submodule (test_pm_matrixInitDia) routines

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatEye_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatEye_IK_ENABLED 1

#if IK4_ENABLED
    module procedure test_getMatEye_IK4_1
        use pm_kind, only: IK, IKC => IK4
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getMatEye_IK3_1
        use pm_kind, only: IK, IKC => IK3
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getMatEye_IK2_1
        use pm_kind, only: IK, IKC => IK2
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getMatEye_IK1_1
        use pm_kind, only: IK, IKC => IK1
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#undef getMatEye_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatEye_CK_ENABLED 1

#if CK3_ENABLED
    module procedure test_getMatEye_CK3
        use pm_kind, only: IK, CK => CK3
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getMatEye_CK2
        use pm_kind, only: IK, CK => CK2
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getMatEye_CK1
        use pm_kind, only: IK, CK => CK1
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#undef getMatEye_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMatEye_RK_ENABLED 1

#if RK3_ENABLED
    module procedure test_getMatEye_RK3
        use pm_kind, only: IK, RK => RK3
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getMatEye_RK2
        use pm_kind, only: IK, RK => RK2
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getMatEye_RK1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#undef getMatEye_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMatEye_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatEye_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatEye_IK_ENABLED 1

#if IK4_ENABLED
    module procedure test_setMatEye_IK4_1
        use pm_kind, only: IK, IKC => IK4
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setMatEye_IK3_1
        use pm_kind, only: IK, IKC => IK3
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setMatEye_IK2_1
        use pm_kind, only: IK, IKC => IK2
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setMatEye_IK1_1
        use pm_kind, only: IK, IKC => IK1
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#undef setMatEye_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatEye_CK_ENABLED 1

#if CK3_ENABLED
    module procedure test_setMatEye_CK3
        use pm_kind, only: IK, CK => CK3
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setMatEye_CK2
        use pm_kind, only: IK, CK => CK2
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setMatEye_CK1
        use pm_kind, only: IK, CK => CK1
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#undef setMatEye_CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMatEye_RK_ENABLED 1

#if RK3_ENABLED
    module procedure test_setMatEye_RK3
        use pm_kind, only: IK, RK => RK3
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setMatEye_RK2
        use pm_kind, only: IK, RK => RK2
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setMatEye_RK1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_matrixInitDia@routines.inc.F90"
    end procedure
#endif

#undef setMatEye_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMatEye_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE