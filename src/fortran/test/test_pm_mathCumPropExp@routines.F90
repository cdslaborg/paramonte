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

!>  \brief This file contains the implementations of the tests of module [pm_mathCumPropExp](@ref pm_mathCumPropExp).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 12:01 AM, August 11, 2021, Dallas, TX

submodule (test_pm_mathCumPropExp) routines

    use pm_arrayReverse, only: setReversed
    use pm_arrayResize, only: setResized
    use pm_mathCumSum, only: getCumSum
    use pm_distUnif, only: getUnifRand
    use pm_val2str, only: getStr

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCumPropExp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Sel_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getCumPropExpSel_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getCumPropExpSel_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getCumPropExpSel_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getCumPropExpSel_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getCumPropExpSel_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Sel_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Seq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getCumPropExpSeq_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getCumPropExpSeq_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getCumPropExpSeq_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getCumPropExpSeq_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getCumPropExpSeq_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Seq_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCumPropExp_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCumPropExp_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Sel_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setCumPropExpSel_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setCumPropExpSel_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setCumPropExpSel_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setCumPropExpSel_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setCumPropExpSel_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Sel_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Seq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setCumPropExpSeq_RK5
        use pm_kind, only: RKC => RK5
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setCumPropExpSeq_RK4
        use pm_kind, only: RKC => RK4
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setCumPropExpSeq_RK3
        use pm_kind, only: RKC => RK3
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setCumPropExpSeq_RK2
        use pm_kind, only: RKC => RK2
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setCumPropExpSeq_RK1
        use pm_kind, only: RKC => RK1
#include "test_pm_mathCumPropExp@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Seq_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCumPropExp_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE