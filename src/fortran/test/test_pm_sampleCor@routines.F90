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
!>  This file contains procedure implementations of tests of [test_pm_sampleCor](@ref test_pm_sampleCor).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_sampleCor) routines ! LCOV_EXCL_LINE

    use pm_kind, only: TKR => RK
    use pm_sampleMean, only: getMean
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_arrayFill, only: getFilled
    use pm_matrixCopy, only: rdpack
    use pm_matrixCopy, only: transHerm
    use pm_matrixInit, only: getMatInit
    use pm_matrixInit, only: setMatInit
    use pm_matrixCopy, only: setMatCopy
    use pm_matrixSubset, only: uppLowDia
    use pm_matrixSubset, only: upp, low, dia
    use pm_matrixSubset, only: uppLowDia_type
    use pm_matrixSubset, only: getSubUnion
    use pm_matrixSubset, only: getSubComp
    use pm_matrixSubset, only: getSubSymm
    use pm_arrayChoice, only: getChoice
    use pm_arrayResize, only: setResized
    use pm_sampleShift, only: getShifted
    use pm_sampleCov, only: getCov, setCov
    use pm_complexCompareAll, only: operator(<)
    use pm_arrayRank, only: getRankFractional
    use pm_container, only: css_type
    !use pm_container, only: css_pdt
    use pm_complexAbs, only: abs
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCor_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getCor_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getCor_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getCor_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getCor_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getCor_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getCor_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getCor_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getCor_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getCor_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getCor_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCor_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCor_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setCor_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setCor_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setCor_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setCor_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setCor_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setCor_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setCor_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setCor_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setCor_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setCor_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCor_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getRho_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRho_D0_SK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRho_D0_SK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRho_D0_SK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRho_D0_SK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRho_D0_SK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XY_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRho_XY_SK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRho_XY_SK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRho_XY_SK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRho_XY_SK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRho_XY_SK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getRho_XY_IK5
        use pm_kind, only: TKG => IK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getRho_XY_IK4
        use pm_kind, only: TKG => IK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getRho_XY_IK3
        use pm_kind, only: TKG => IK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getRho_XY_IK2
        use pm_kind, only: TKG => IK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getRho_XY_IK1
        use pm_kind, only: TKG => IK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRho_XY_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRho_XY_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRho_XY_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRho_XY_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRho_XY_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRho_XY_PSSK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRho_XY_PSSK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRho_XY_PSSK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRho_XY_PSSK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRho_XY_PSSK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure test_getRho_XY_BSSK
        use pm_kind, only: TKG => SK
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XY_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRho_D2_SK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRho_D2_SK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRho_D2_SK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRho_D2_SK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRho_D2_SK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_getRho_D2_IK5
        use pm_kind, only: TKG => IK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_getRho_D2_IK4
        use pm_kind, only: TKG => IK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_getRho_D2_IK3
        use pm_kind, only: TKG => IK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_getRho_D2_IK2
        use pm_kind, only: TKG => IK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_getRho_D2_IK1
        use pm_kind, only: TKG => IK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getRho_D2_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getRho_D2_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getRho_D2_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getRho_D2_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getRho_D2_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_getRho_D2_PSSK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_getRho_D2_PSSK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_getRho_D2_PSSK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_getRho_D2_PSSK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_getRho_D2_PSSK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure test_getRho_D2_BSSK
        use pm_kind, only: TKG => SK
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getRho_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setRho_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRho_D0_SK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRho_D0_SK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRho_D0_SK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRho_D0_SK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRho_D0_SK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XY_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRho_XY_SK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRho_XY_SK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRho_XY_SK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRho_XY_SK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRho_XY_SK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setRho_XY_IK5
        use pm_kind, only: TKG => IK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setRho_XY_IK4
        use pm_kind, only: TKG => IK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setRho_XY_IK3
        use pm_kind, only: TKG => IK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setRho_XY_IK2
        use pm_kind, only: TKG => IK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setRho_XY_IK1
        use pm_kind, only: TKG => IK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRho_XY_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRho_XY_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRho_XY_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRho_XY_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRho_XY_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRho_XY_PSSK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRho_XY_PSSK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRho_XY_PSSK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRho_XY_PSSK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRho_XY_PSSK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure test_setRho_XY_BSSK
        use pm_kind, only: TKG => SK
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XY_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRho_D2_SK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRho_D2_SK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRho_D2_SK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRho_D2_SK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRho_D2_SK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setRho_D2_IK5
        use pm_kind, only: TKG => IK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setRho_D2_IK4
        use pm_kind, only: TKG => IK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setRho_D2_IK3
        use pm_kind, only: TKG => IK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setRho_D2_IK2
        use pm_kind, only: TKG => IK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setRho_D2_IK1
        use pm_kind, only: TKG => IK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setRho_D2_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setRho_D2_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setRho_D2_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setRho_D2_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setRho_D2_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setRho_D2_PSSK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setRho_D2_PSSK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setRho_D2_PSSK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setRho_D2_PSSK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setRho_D2_PSSK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure test_setRho_D2_BSSK
        use pm_kind, only: TKG => SK
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setRho_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCordance_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setCordance_D0_SK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setCordance_D0_SK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setCordance_D0_SK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setCordance_D0_SK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setCordance_D0_SK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setCordance_D1_SK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setCordance_D1_SK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setCordance_D1_SK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setCordance_D1_SK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setCordance_D1_SK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure test_setCordance_D1_IK5
        use pm_kind, only: TKG => IK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure test_setCordance_D1_IK4
        use pm_kind, only: TKG => IK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure test_setCordance_D1_IK3
        use pm_kind, only: TKG => IK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure test_setCordance_D1_IK2
        use pm_kind, only: TKG => IK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure test_setCordance_D1_IK1
        use pm_kind, only: TKG => IK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setCordance_D1_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setCordance_D1_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setCordance_D1_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setCordance_D1_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setCordance_D1_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED
#define PSSK_ENABLED 1

#if SK5_ENABLED
    module procedure test_setCordance_D1_PSSK5
        use pm_kind, only: TKG => SK5
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure test_setCordance_D1_PSSK4
        use pm_kind, only: TKG => SK4
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure test_setCordance_D1_PSSK3
        use pm_kind, only: TKG => SK3
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure test_setCordance_D1_PSSK2
        use pm_kind, only: TKG => SK2
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure test_setCordance_D1_PSSK1
        use pm_kind, only: TKG => SK1
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure
#endif

#undef PSSK_ENABLED
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define BSSK_ENABLED 1

    module procedure test_setCordance_D1_BSSK
        use pm_kind, only: TKG => SK
#include "test_pm_sampleCor@routines.inc.F90"
    end procedure

#undef BSSK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCordance_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines