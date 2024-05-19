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
!>  This file contains procedure implementations of [pm_sampleECDF](@ref pm_sampleECDF).
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

submodule (test_pm_sampleECDF) routines ! LCOV_EXCL_LINE

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setECDF_D1_RK_ENABLED 1

#if RK3_ENABLED
    module procedure test_setECDF_D1_RK3_RK3
        use pm_kind, only: IK, RK => RK3
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure test_setECDF_D1_RK2_RK2
        use pm_kind, only: IK, RK => RK2
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure test_setECDF_D1_RK1_RK1
        use pm_kind, only: IK, RK => RK1
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#undef setECDF_D1_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setECDF_D1_IK_ENABLED 1

#if IK4_ENABLED && RK3_ENABLED
    module procedure test_setECDF_D1_IK4_RK3
        use pm_kind, only: IK, IKG => IK4, RK => RK3
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && RK2_ENABLED
    module procedure test_setECDF_D1_IK4_RK2
        use pm_kind, only: IK, IKG => IK4, RK => RK2
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#if IK4_ENABLED && RK1_ENABLED
    module procedure test_setECDF_D1_IK4_RK1
        use pm_kind, only: IK, IKG => IK4, RK => RK1
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK3_ENABLED && RK3_ENABLED
    module procedure test_setECDF_D1_IK3_RK3
        use pm_kind, only: IK, IKG => IK3, RK => RK3
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && RK2_ENABLED
    module procedure test_setECDF_D1_IK3_RK2
        use pm_kind, only: IK, IKG => IK3, RK => RK2
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#if IK3_ENABLED && RK1_ENABLED
    module procedure test_setECDF_D1_IK3_RK1
        use pm_kind, only: IK, IKG => IK3, RK => RK1
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK2_ENABLED && RK3_ENABLED
    module procedure test_setECDF_D1_IK2_RK3
        use pm_kind, only: IK, IKG => IK2, RK => RK3
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && RK2_ENABLED
    module procedure test_setECDF_D1_IK2_RK2
        use pm_kind, only: IK, IKG => IK2, RK => RK2
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#if IK2_ENABLED && RK1_ENABLED
    module procedure test_setECDF_D1_IK2_RK1
        use pm_kind, only: IK, IKG => IK2, RK => RK1
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#if IK1_ENABLED && RK3_ENABLED
    module procedure test_setECDF_D1_IK1_RK3
        use pm_kind, only: IK, IKG => IK1, RK => RK3
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && RK2_ENABLED
    module procedure test_setECDF_D1_IK1_RK2
        use pm_kind, only: IK, IKG => IK1, RK => RK2
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#if IK1_ENABLED && RK1_ENABLED
    module procedure test_setECDF_D1_IK1_RK1
        use pm_kind, only: IK, IKG => IK1, RK => RK1
#include "test_pm_sampleECDF@routines@ECDF.inc.F90"
    end procedure
#endif

#undef setECDF_D1_IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines
