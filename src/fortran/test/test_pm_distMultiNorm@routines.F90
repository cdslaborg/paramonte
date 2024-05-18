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

!>  \brief This file contains the implementations of the tests of module [pm_distMultiNorm](@ref pm_distMultiNorm).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distMultiNorm) routines

    use pm_val2str, only: getStr
    use pm_matrix2vec, only: getVecDia
    use pm_sampleMean, only: getMean
    use pm_matrixChol, only: setMatChol, uppDia
    use pm_sampleCov, only: getCovUpp
    use pm_sampleShift, only: getShifted
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getMultiNormRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getMultiNormRand_RK5
        use pm_kind, only: RKC =>RK5
#include "test_pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getMultiNormRand_RK4
        use pm_kind, only: RKC =>RK4
#include "test_pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getMultiNormRand_RK3
        use pm_kind, only: RKC =>RK3
#include "test_pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getMultiNormRand_RK2
        use pm_kind, only: RKC =>RK2
#include "test_pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getMultiNormRand_RK1
        use pm_kind, only: RKC =>RK1
#include "test_pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getMultiNormRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setMultiNormRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setMultiNormRand_RK5
        use pm_kind, only: RKC =>RK5
#include "test_pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setMultiNormRand_RK4
        use pm_kind, only: RKC =>RK4
#include "test_pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setMultiNormRand_RK3
        use pm_kind, only: RKC =>RK3
#include "test_pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setMultiNormRand_RK2
        use pm_kind, only: RKC =>RK2
#include "test_pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setMultiNormRand_RK1
        use pm_kind, only: RKC =>RK1
#include "test_pm_distMultiNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setMultiNormRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE