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

!>  \brief This file contains the implementations of the tests of module [pm_distNorm](@ref pm_distNorm).
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (test_pm_distNorm) routines

    use pm_sampleMean, only: getMean
    use pm_sampleVar, only: getVar
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getNormCDF_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getNormCDF_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getNormCDF_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getNormCDF_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getNormCDF_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNormCDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setNormCDF_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setNormCDF_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setNormCDF_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setNormCDF_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setNormCDF_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNormCDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getNormLogPDF_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getNormLogPDF_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getNormLogPDF_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getNormLogPDF_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getNormLogPDF_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNormLogPDF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setNormLogPDF_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setNormLogPDF_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setNormLogPDF_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setNormLogPDF_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setNormLogPDF_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNormLogPDF_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getNormRand_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getNormRand_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getNormRand_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getNormRand_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getNormRand_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNormRand_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setNormRand_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setNormRand_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setNormRand_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setNormRand_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setNormRand_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNormRand_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getNormQuan_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getNormQuan_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getNormQuan_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getNormQuan_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getNormQuan_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setNormQuan_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_setNormQuan_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_setNormQuan_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_setNormQuan_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setNormQuan_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setNormQuan_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setNormQuan_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormEntropy_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getNormEntropy_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getNormEntropy_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getNormEntropy_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getNormEntropy_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getNormEntropy_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormEntropy_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormFisher_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getNormFisher_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getNormFisher_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getNormFisher_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getNormFisher_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getNormFisher_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormFisher_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getNormKLD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure test_getNormKLD_RK5
        use pm_kind, only: RKG => RK5
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure test_getNormKLD_RK4
        use pm_kind, only: RKG => RK4
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure test_getNormKLD_RK3
        use pm_kind, only: RKG => RK3
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_getNormKLD_RK2
        use pm_kind, only: RKG => RK2
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_getNormKLD_RK1
        use pm_kind, only: RKG => RK1
#include "test_pm_distNorm@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getNormKLD_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines ! LCOV_EXCL_LINE