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
!>  This file contains procedure implementations of [pm_distanceEuclid](@ref pm_distanceEuclid).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Wednesday 5:43 PM, December 25, 2013, Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distanceEuclid) routines

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_ellipsoid, only: setVolUnitBall
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisEuclid_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define XYZ_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisEuclidXYZ_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisEuclidXYZ_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisEuclidXYZ_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisEuclidXYZ_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisEuclidXYZ_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef XYZ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDE_D1_XX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDE_D1_XX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDE_D1_XX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDE_D1_XX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDE_D1_XX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDE_D2_XX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDE_D2_XX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDE_D2_XX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDE_D2_XX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDE_D2_XX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDE_D1_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDE_D1_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDE_D1_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDE_D1_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDE_D1_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDE_D1_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDE_D1_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDE_D1_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDE_D1_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDE_D1_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDE_D2_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDE_D2_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDE_D2_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDE_D2_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDE_D2_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDE_D2_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDE_D2_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDE_D2_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDE_D2_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDE_D2_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisEuclid_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setDisEuclid_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MED_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MED_D0_D1_XX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MED_D0_D1_XX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MED_D0_D1_XX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MED_D0_D1_XX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MED_D0_D1_XX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MED_D1_D2_XX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MED_D1_D2_XX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MED_D1_D2_XX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MED_D1_D2_XX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MED_D1_D2_XX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MED_D0_D1_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MED_D0_D1_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MED_D0_D1_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MED_D0_D1_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MED_D0_D1_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MED_D1_D1_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MED_D1_D1_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MED_D1_D1_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MED_D1_D1_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MED_D1_D1_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MED_D1_D2_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MED_D1_D2_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MED_D1_D2_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MED_D1_D2_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MED_D1_D2_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D2_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MED_D2_D2_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MED_D2_D2_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MED_D2_D2_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MED_D2_D2_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MED_D2_D2_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D2_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MED_D2_D1_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MED_D2_D1_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MED_D2_D1_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MED_D2_D1_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MED_D2_D1_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MED_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEU_D0_D1_XX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEU_D0_D1_XX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEU_D0_D1_XX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEU_D0_D1_XX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEU_D0_D1_XX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEU_D1_D2_XX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEU_D1_D2_XX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEU_D1_D2_XX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEU_D1_D2_XX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEU_D1_D2_XX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEU_D0_D1_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEU_D0_D1_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEU_D0_D1_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEU_D0_D1_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEU_D0_D1_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEU_D1_D1_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEU_D1_D1_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEU_D1_D1_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEU_D1_D1_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEU_D1_D1_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEU_D1_D2_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEU_D1_D2_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEU_D1_D2_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEU_D1_D2_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEU_D1_D2_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D2_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEU_D2_D2_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEU_D2_D2_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEU_D2_D2_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEU_D2_D2_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEU_D2_D2_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D2_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEU_D2_D1_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEU_D2_D1_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEU_D2_D1_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEU_D2_D1_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEU_D2_D1_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEV_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEV_D0_D1_XX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEV_D0_D1_XX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEV_D0_D1_XX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEV_D0_D1_XX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEV_D0_D1_XX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEV_D1_D2_XX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEV_D1_D2_XX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEV_D1_D2_XX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEV_D1_D2_XX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEV_D1_D2_XX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEV_D0_D1_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEV_D0_D1_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEV_D0_D1_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEV_D0_D1_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEV_D0_D1_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEV_D1_D1_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEV_D1_D1_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEV_D1_D1_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEV_D1_D1_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEV_D1_D1_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEV_D1_D2_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEV_D1_D2_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEV_D1_D2_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEV_D1_D2_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEV_D1_D2_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D2_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEV_D2_D2_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEV_D2_D2_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEV_D2_D2_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEV_D2_D2_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEV_D2_D2_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D2_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEV_D2_D1_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEV_D2_D1_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEV_D2_D1_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEV_D2_D1_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEV_D2_D1_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEV_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEQ_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEQ_D0_D1_XX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEQ_D0_D1_XX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEQ_D0_D1_XX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEQ_D0_D1_XX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEQ_D0_D1_XX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_XX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEQ_D1_D2_XX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEQ_D1_D2_XX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEQ_D1_D2_XX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEQ_D1_D2_XX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEQ_D1_D2_XX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_XX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D0_D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEQ_D0_D1_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEQ_D0_D1_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEQ_D0_D1_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEQ_D0_D1_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEQ_D0_D1_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D0_D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D1_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEQ_D1_D1_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEQ_D1_D1_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEQ_D1_D1_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEQ_D1_D1_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEQ_D1_D1_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D1_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_D2_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEQ_D1_D2_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEQ_D1_D2_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEQ_D1_D2_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEQ_D1_D2_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEQ_D1_D2_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_D2_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D2_D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEQ_D2_D2_D2_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEQ_D2_D2_D2_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEQ_D2_D2_D2_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEQ_D2_D2_D2_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEQ_D2_D2_D2_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D2_D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_D1_D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDE_MEQ_D2_D1_D1_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDE_MEQ_D2_D1_D1_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDE_MEQ_D2_D1_D1_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDE_MEQ_D2_D1_D1_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDE_MEQ_D2_D1_D1_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_D1_D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEQ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setDisEuclid_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisMatEuclid_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define FUL_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDME_RDP_FUL_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDME_RDP_FUL_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDME_RDP_FUL_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDME_RDP_FUL_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDME_RDP_FUL_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef FUL_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDME_RDP_ULD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDME_RDP_ULD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDME_RDP_ULD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDME_RDP_ULD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDME_RDP_ULD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDME_RDP_ULX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDME_RDP_ULX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDME_RDP_ULX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDME_RDP_ULX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDME_RDP_ULX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisMatEuclid_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setDisMatEuclid_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MED_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDME_MED_RDP_ULD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDME_MED_RDP_ULD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDME_MED_RDP_ULD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDME_MED_RDP_ULD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDME_MED_RDP_ULD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDME_MED_RDP_ULX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDME_MED_RDP_ULX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDME_MED_RDP_ULX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDME_MED_RDP_ULX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDME_MED_RDP_ULX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MED_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEU_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDME_MEU_RDP_ULD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDME_MEU_RDP_ULD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDME_MEU_RDP_ULD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDME_MEU_RDP_ULD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDME_MEU_RDP_ULD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDME_MEU_RDP_ULX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDME_MEU_RDP_ULX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDME_MEU_RDP_ULX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDME_MEU_RDP_ULX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDME_MEU_RDP_ULX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEU_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define MEQ_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RDP_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULD_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDME_MEQ_RDP_ULD_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDME_MEQ_RDP_ULD_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDME_MEQ_RDP_ULD_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDME_MEQ_RDP_ULD_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDME_MEQ_RDP_ULD_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULD_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ULX_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setDME_MEQ_RDP_ULX_RK5
        use pm_kind, only: TKG => RK5
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setDME_MEQ_RDP_ULX_RK4
        use pm_kind, only: TKG => RK4
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setDME_MEQ_RDP_ULX_RK3
        use pm_kind, only: TKG => RK3
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setDME_MEQ_RDP_ULX_RK2
        use pm_kind, only: TKG => RK2
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setDME_MEQ_RDP_ULX_RK1
        use pm_kind, only: TKG => RK1
#include "pm_distanceEuclid@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ULX_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef RDP_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef MEQ_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setDisMatEuclid_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines