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
!>  This file contains procedure implementations of [pm_distanceHellinger](@ref pm_distanceHellinger).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_distanceHellinger) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    use pm_quadPack, only: isFailedQuad, weps
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getDisHellSq_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PMF_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisHellSqPMF_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceHellinger@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisHellSqPMF_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceHellinger@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisHellSqPMF_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceHellinger@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisHellSqPMF_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceHellinger@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisHellSqPMF_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceHellinger@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PMF_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define PRC_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getDisHellSqPRC_RK5
        use pm_kind, only: RKC => RK5
#include "pm_distanceHellinger@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getDisHellSqPRC_RK4
        use pm_kind, only: RKC => RK4
#include "pm_distanceHellinger@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getDisHellSqPRC_RK3
        use pm_kind, only: RKC => RK3
#include "pm_distanceHellinger@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getDisHellSqPRC_RK2
        use pm_kind, only: RKC => RK2
#include "pm_distanceHellinger@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getDisHellSqPRC_RK1
        use pm_kind, only: RKC => RK1
#include "pm_distanceHellinger@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PRC_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getDisHellSq_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines