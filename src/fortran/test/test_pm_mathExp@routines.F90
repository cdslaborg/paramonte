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
!>  This file contains the implementations of the tests of module [pm_mathExp](@ref pm_mathExp).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_mathExp) routines

    use pm_kind, only: LK
    use pm_val2str, only: getStr
    use pm_option, only: getOption
    use pm_distUnif, only: getUnifRand
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpNext_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getExpNext_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getExpNext_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getExpNext_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getExpNext_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getExpNext_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getExpNext_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getExpNext_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getExpNext_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getExpNext_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getExpNext_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpNext_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getExpPrev_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getExpPrev_IK5
        use pm_kind, only: IKG => IK5
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getExpPrev_IK4
        use pm_kind, only: IKG => IK4
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getExpPrev_IK3
        use pm_kind, only: IKG => IK3
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getExpPrev_IK2
        use pm_kind, only: IKG => IK2
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getExpPrev_IK1
        use pm_kind, only: IKG => IK1
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getExpPrev_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getExpPrev_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getExpPrev_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getExpPrev_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getExpPrev_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_mathExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getExpPrev_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE