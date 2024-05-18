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
!>  This file contains procedure implementations of tests of [test_pm_sampleCCF](@ref test_pm_sampleCCF).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_sampleCCF) routines ! LCOV_EXCL_LINE

    use pm_io, only: getFormat
    use pm_fftpack, only: getFFTF
    use pm_fftpack, only: getFFTI
    use pm_fftpack, only: getFactorFFT
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
    use pm_arrayRange, only: getRange
    use pm_sampleMean, only: getMean
    use pm_sampleCov, only: getCov
    use pm_sampleCov, only: setCov
    use pm_sampleShift, only: getShifted
    use pm_complexCompareAll, only: operator(<)
    use pm_container, only: css_type, css_pdt
    use pm_arrayRank, only: getRankFractional
    use pm_complexAbs, only: abs
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCCF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getCCF_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getCCF_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getCCF_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getCCF_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getCCF_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getCCF_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getCCF_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getCCF_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getCCF_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getCCF_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCCF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCCF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setCCF_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setCCF_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setCCF_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setCCF_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setCCF_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setCCF_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setCCF_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setCCF_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setCCF_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setCCF_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCCF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getACF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getACF_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getACF_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getACF_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getACF_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getACF_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getACF_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getACF_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getACF_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getACF_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getACF_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getACF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setACF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setACF_CK5
        use pm_kind, only: TKC => CK5
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setACF_CK4
        use pm_kind, only: TKC => CK4
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setACF_CK3
        use pm_kind, only: TKC => CK3
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setACF_CK2
        use pm_kind, only: TKC => CK2
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setACF_CK1
        use pm_kind, only: TKC => CK1
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setACF_RK5
        use pm_kind, only: TKC => RK5
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setACF_RK4
        use pm_kind, only: TKC => RK4
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setACF_RK3
        use pm_kind, only: TKC => RK3
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setACF_RK2
        use pm_kind, only: TKC => RK2
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setACF_RK1
        use pm_kind, only: TKC => RK1
#include "test_pm_sampleCCF@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setACF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines