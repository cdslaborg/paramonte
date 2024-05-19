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

!>  \brief This file contains the implementations of the tests of module [pm_mathCumSum](@ref pm_mathCumSum).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 12:01 AM, August 11, 2021, Dallas, TX

submodule (test_pm_mathCumSum) routines

    use pm_complexCompareAll, only: operator(<=)
    use pm_arrayReverse, only: setReversed
    use pm_arrayResize, only: setResized
    use pm_mathCumSum, only: getCumSum
    use pm_distUnif, only: getUnifRand
    use pm_val2str, only: getStr
    
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCumSum_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getCumSum_IK5
        use pm_kind, only: TKG => IK5
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getCumSum_IK4
        use pm_kind, only: TKG => IK4
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getCumSum_IK3
        use pm_kind, only: TKG => IK3
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getCumSum_IK2
        use pm_kind, only: TKG => IK2
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getCumSum_IK1
        use pm_kind, only: TKG => IK1
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getCumSum_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getCumSum_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getCumSum_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getCumSum_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getCumSum_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getCumSum_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getCumSum_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getCumSum_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getCumSum_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getCumSum_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCumSum_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCumSum_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setCumSum_IK5
        use pm_kind, only: TKG => IK5
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setCumSum_IK4
        use pm_kind, only: TKG => IK4
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setCumSum_IK3
        use pm_kind, only: TKG => IK3
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setCumSum_IK2
        use pm_kind, only: TKG => IK2
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setCumSum_IK1
        use pm_kind, only: TKG => IK1
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setCumSum_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setCumSum_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setCumSum_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setCumSum_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setCumSum_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setCumSum_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setCumSum_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setCumSum_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setCumSum_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setCumSum_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_mathCumSum@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCumSum_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE