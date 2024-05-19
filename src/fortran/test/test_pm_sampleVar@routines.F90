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
!>  This file contains procedure implementations of tests of [test_pm_sampleVar](@ref test_pm_sampleVar).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_sampleVar) routines ! LCOV_EXCL_LINE

    use pm_sampleMean, only: getMean
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_arrayResize, only: setResized
    use pm_arrayChoice, only: getChoice
    use pm_sampleShift, only: getShifted
    use pm_complexCompareAll, only: operator(<)
    use pm_complexAbs, only: abs
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getVarCorrection_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getVarCorrection_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getVarCorrection_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getVarCorrection_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getVarCorrection_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getVarCorrection_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVarCorrection_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getVar_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getVar_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getVar_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getVar_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getVar_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getVar_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getVar_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getVar_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getVar_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getVar_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getVar_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVar_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setVar_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setVar_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setVar_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setVar_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setVar_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setVar_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setVar_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setVar_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setVar_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setVar_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setVar_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVar_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setVarMean_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setVarMean_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setVarMean_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setVarMean_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setVarMean_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setVarMean_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setVarMean_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setVarMean_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setVarMean_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setVarMean_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setVarMean_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVarMean_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getVarMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_getVarMerged_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_getVarMerged_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_getVarMerged_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_getVarMerged_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_getVarMerged_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getVarMerged_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getVarMerged_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getVarMerged_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getVarMerged_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getVarMerged_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getVarMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setVarMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setVarMerged_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setVarMerged_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setVarMerged_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setVarMerged_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setVarMerged_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setVarMerged_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setVarMerged_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setVarMerged_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setVarMerged_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setVarMerged_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVarMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setVarMeanMerged_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure test_setVarMeanMerged_CK5
        use pm_kind, only: TKG => CK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure test_setVarMeanMerged_CK4
        use pm_kind, only: TKG => CK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure test_setVarMeanMerged_CK3
        use pm_kind, only: TKG => CK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure test_setVarMeanMerged_CK2
        use pm_kind, only: TKG => CK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure test_setVarMeanMerged_CK1
        use pm_kind, only: TKG => CK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setVarMeanMerged_RK5
        use pm_kind, only: TKG => RK5
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setVarMeanMerged_RK4
        use pm_kind, only: TKG => RK4
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setVarMeanMerged_RK3
        use pm_kind, only: TKG => RK3
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setVarMeanMerged_RK2
        use pm_kind, only: TKG => RK2
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setVarMeanMerged_RK1
        use pm_kind, only: TKG => RK1
#include "test_pm_sampleVar@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setVarMeanMerged_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines