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
!>  This file contains procedure implementations of [pm_mathRoot](@ref pm_mathRoot).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_mathRoot) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG)
#endif

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
    module procedure getRootDef_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRootDef_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRootDef_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRootDef_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRootDef_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Def_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define False_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getRootFalse_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRootFalse_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRootFalse_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRootFalse_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRootFalse_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
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
    module procedure getRootBisection_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRootBisection_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRootBisection_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRootBisection_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRootBisection_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
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
    module procedure getRootSecant_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRootSecant_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRootSecant_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRootSecant_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRootSecant_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
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
    module procedure getRootBrent_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRootBrent_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRootBrent_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRootBrent_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRootBrent_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
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
    module procedure getRootRidders_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRootRidders_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRootRidders_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRootRidders_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRootRidders_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
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
    module procedure getRootTOMS748_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRootTOMS748_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRootTOMS748_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRootTOMS748_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRootTOMS748_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
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
    module procedure getRootNewton_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRootNewton_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRootNewton_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRootNewton_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRootNewton_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
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
    module procedure getRootHalley_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRootHalley_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRootHalley_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRootHalley_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRootHalley_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
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
    module procedure getRootSchroder_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getRootSchroder_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getRootSchroder_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getRootSchroder_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getRootSchroder_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
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
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootFalseFixed_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootFalseFixed_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootFalseFixed_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootFalseFixed_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootFalseFixed_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Niter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootFalseNiter_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootFalseNiter_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootFalseNiter_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootFalseNiter_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootFalseNiter_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Niter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef False_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Bisection_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootBisectionFixed_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootBisectionFixed_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootBisectionFixed_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootBisectionFixed_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootBisectionFixed_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Niter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootBisectionNiter_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootBisectionNiter_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootBisectionNiter_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootBisectionNiter_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootBisectionNiter_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Niter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Bisection_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Secant_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootSecantFixed_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootSecantFixed_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootSecantFixed_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootSecantFixed_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootSecantFixed_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Niter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootSecantNiter_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootSecantNiter_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootSecantNiter_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootSecantNiter_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootSecantNiter_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Niter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Secant_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Brent_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootBrentFixed_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootBrentFixed_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootBrentFixed_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootBrentFixed_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootBrentFixed_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Niter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootBrentNiter_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootBrentNiter_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootBrentNiter_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootBrentNiter_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootBrentNiter_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Niter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Brent_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Ridders_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootRiddersFixed_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootRiddersFixed_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootRiddersFixed_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootRiddersFixed_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootRiddersFixed_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Niter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootRiddersNiter_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootRiddersNiter_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootRiddersNiter_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootRiddersNiter_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootRiddersNiter_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Niter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Ridders_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define TOMS748_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootTOMS748Fixed_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootTOMS748Fixed_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootTOMS748Fixed_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootTOMS748Fixed_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootTOMS748Fixed_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Niter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootTOMS748Niter_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootTOMS748Niter_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootTOMS748Niter_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootTOMS748Niter_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootTOMS748Niter_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Niter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef TOMS748_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Newton_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootNewtonFixed_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootNewtonFixed_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootNewtonFixed_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootNewtonFixed_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootNewtonFixed_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Niter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootNewtonNiter_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootNewtonNiter_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootNewtonNiter_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootNewtonNiter_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootNewtonNiter_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Niter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Newton_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Halley_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootHalleyFixed_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootHalleyFixed_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootHalleyFixed_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootHalleyFixed_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootHalleyFixed_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Niter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootHalleyNiter_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootHalleyNiter_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootHalleyNiter_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootHalleyNiter_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootHalleyNiter_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Niter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Halley_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Schroder_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Fixed_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootSchroderFixed_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootSchroderFixed_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootSchroderFixed_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootSchroderFixed_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootSchroderFixed_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Fixed_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Niter_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setRootSchroderNiter_RK5
        use pm_kind, only: RKC => RK5
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setRootSchroderNiter_RK4
        use pm_kind, only: RKC => RK4
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setRootSchroderNiter_RK3
        use pm_kind, only: RKC => RK3
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setRootSchroderNiter_RK2
        use pm_kind, only: RKC => RK2
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setRootSchroderNiter_RK1
        use pm_kind, only: RKC => RK1
#include "pm_mathRoot@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Niter_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Schroder_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRoot_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines