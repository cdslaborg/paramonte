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
!>  This file contains the implementations of the tests of module [test_pm_mathRoot](@ref test_pm_mathRoot).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_mathRoot) routines

    use pm_arrayMembership, only: operator(.allinrange.)
    use pm_arrayRemove, only: getRemoved
    use pm_arrayChoice, only: getChoice
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type
    use pm_arrayRange, only: getRange
    use pm_option, only: getOption
    use pm_val2str, only: getStr

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRoot_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Def_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRootDef_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRootDef_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRootDef_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRootDef_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRootDef_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRoot_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRoot_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define False_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRootFalse_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRootFalse_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRootFalse_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRootFalse_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRootFalse_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef False_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Bisection_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRootBisection_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRootBisection_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRootBisection_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRootBisection_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRootBisection_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Bisection_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Secant_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRootSecant_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRootSecant_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRootSecant_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRootSecant_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRootSecant_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Secant_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Brent_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRootBrent_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRootBrent_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRootBrent_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRootBrent_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRootBrent_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Brent_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ridders_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRootRidders_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRootRidders_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRootRidders_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRootRidders_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRootRidders_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ridders_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TOMS748_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRootTOMS748_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRootTOMS748_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRootTOMS748_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRootTOMS748_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRootTOMS748_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TOMS748_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Newton_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRootNewton_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRootNewton_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRootNewton_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRootNewton_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRootNewton_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Newton_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Halley_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRootHalley_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRootHalley_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRootHalley_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRootHalley_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRootHalley_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Halley_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Schroder_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRootSchroder_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRootSchroder_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRootSchroder_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRootSchroder_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRootSchroder_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Schroder_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRoot_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRoot_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define False_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRootFalse_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRootFalse_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRootFalse_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRootFalse_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRootFalse_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef False_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Bisection_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRootBisection_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRootBisection_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRootBisection_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRootBisection_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRootBisection_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Bisection_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Secant_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRootSecant_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRootSecant_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRootSecant_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRootSecant_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRootSecant_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Secant_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Brent_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRootBrent_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRootBrent_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRootBrent_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRootBrent_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRootBrent_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Brent_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ridders_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRootRidders_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRootRidders_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRootRidders_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRootRidders_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRootRidders_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ridders_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TOMS748_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRootTOMS748_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRootTOMS748_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRootTOMS748_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRootTOMS748_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRootTOMS748_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TOMS748_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Newton_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRootNewton_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRootNewton_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRootNewton_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRootNewton_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRootNewton_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Newton_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Halley_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRootHalley_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRootHalley_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRootHalley_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRootHalley_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRootHalley_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Halley_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Schroder_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRootSchroder_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRootSchroder_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRootSchroder_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRootSchroder_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRootSchroder_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Schroder_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRoot_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE