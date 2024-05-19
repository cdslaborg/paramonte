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
!>  This file contains procedure implementations of [pm_sampleWeight](@ref pm_sampleWeight).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sampleWeight) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getReweight_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure getReweight_IK_IK5
        use pm_kind, only: IKG => IK5
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure getReweight_IK_IK4
        use pm_kind, only: IKG => IK4
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure getReweight_IK_IK3
        use pm_kind, only: IKG => IK3
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure getReweight_IK_IK2
        use pm_kind, only: IKG => IK2
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure getReweight_IK_IK1
        use pm_kind, only: IKG => IK1
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getReweight_RK_RK5
        use pm_kind, only: RKG => RK5
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getReweight_RK_RK4
        use pm_kind, only: RKG => RK4
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getReweight_RK_RK3
        use pm_kind, only: RKG => RK3
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getReweight_RK_RK2
        use pm_kind, only: RKG => RK2
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getReweight_RK_RK1
        use pm_kind, only: RKG => RK1
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_RK_ENABLED 1

#if RK5_ENABLED
    module procedure getReweight_IK_RK5
        use pm_kind, only: RKG => RK5
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getReweight_IK_RK4
        use pm_kind, only: RKG => RK4
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getReweight_IK_RK3
        use pm_kind, only: RKG => RK3
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getReweight_IK_RK2
        use pm_kind, only: RKG => RK2
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getReweight_IK_RK1
        use pm_kind, only: RKG => RK1
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#undef IK_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getReweight_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setReweight_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_ENABLED 1

#if IK5_ENABLED
    module procedure setReweight_IK_IK5
        use pm_kind, only: IKG => IK5
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if IK4_ENABLED
    module procedure setReweight_IK_IK4
        use pm_kind, only: IKG => IK4
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if IK3_ENABLED
    module procedure setReweight_IK_IK3
        use pm_kind, only: IKG => IK3
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if IK2_ENABLED
    module procedure setReweight_IK_IK2
        use pm_kind, only: IKG => IK2
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if IK1_ENABLED
    module procedure setReweight_IK_IK1
        use pm_kind, only: IKG => IK1
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#undef IK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReweight_RK_RK5
        use pm_kind, only: RKG => RK5
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReweight_RK_RK4
        use pm_kind, only: RKG => RK4
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReweight_RK_RK3
        use pm_kind, only: RKG => RK3
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReweight_RK_RK2
        use pm_kind, only: RKG => RK2
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReweight_RK_RK1
        use pm_kind, only: RKG => RK1
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define IK_RK_ENABLED 1

#if RK5_ENABLED
    module procedure setReweight_IK_RK5
        use pm_kind, only: IKG => IK, RKG => RK5
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setReweight_IK_RK4
        use pm_kind, only: IKG => IK, RKG => RK4
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setReweight_IK_RK3
        use pm_kind, only: IKG => IK, RKG => RK3
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setReweight_IK_RK2
        use pm_kind, only: IKG => IK, RKG => RK2
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setReweight_IK_RK1
        use pm_kind, only: IKG => IK, RKG => RK1
#include "pm_sampleWeight@routines.inc.F90"
    end procedure
#endif

#undef IK_RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setReweight_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines