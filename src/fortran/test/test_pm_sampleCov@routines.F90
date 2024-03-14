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
!>  This file contains procedure implementations of tests of [test_pm_sampleCov](@ref test_pm_sampleCov).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_sampleCov) routines ! LCOV_EXCL_LINE

    use pm_sampleMean, only: getMean
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_arrayFill, only: getFilled
    use pm_matrixCopy, only: rdpack
    use pm_matrixCopy, only: transHerm
    use pm_matrixInit, only: setMatInit
    use pm_matrixCopy, only: setMatCopy
    use pm_matrixSubset, only: uppLowDia
    use pm_matrixSubset, only: uppLowDia_type
    use pm_matrixSubset, only: getSubComp
    use pm_matrixSubset, only: getSubSymm
    use pm_arrayChoice, only: getChoice
    use pm_arrayResize, only: setResized
    use pm_sampleShift, only: getShifted
    use pm_complexCompareAll, only: operator(<)
    use pm_complexAbs, only: abs
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCov_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getCov_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getCov_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getCov_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getCov_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getCov_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getCov_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getCov_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getCov_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getCov_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getCov_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCov_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCov_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setCov_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setCov_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setCov_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setCov_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setCov_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setCov_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setCov_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setCov_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setCov_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setCov_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCov_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setCovMean_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setCovMean_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setCovMean_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setCovMean_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setCovMean_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setCovMean_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setCovMean_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setCovMean_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setCovMean_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setCovMean_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCovMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getCovMerged_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getCovMerged_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getCovMerged_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getCovMerged_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getCovMerged_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getCovMerged_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getCovMerged_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getCovMerged_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getCovMerged_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getCovMerged_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCovMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setCovMerged_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setCovMerged_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setCovMerged_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setCovMerged_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setCovMerged_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setCovMerged_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setCovMerged_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setCovMerged_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setCovMerged_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setCovMerged_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCovMeanMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setCovMeanMerged_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setCovMeanMerged_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setCovMeanMerged_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setCovMeanMerged_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setCovMeanMerged_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setCovMeanMerged_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setCovMeanMerged_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setCovMeanMerged_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setCovMeanMerged_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setCovMeanMerged_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_sampleCov@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCovMeanMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines