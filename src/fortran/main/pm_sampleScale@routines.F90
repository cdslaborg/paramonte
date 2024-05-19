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
!>  This file contains procedure implementations of [pm_sampleScale](@ref pm_sampleScale).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_sampleScale) routines ! LCOV_EXCL_LINE

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

#define getScaled_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ARK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getScaled_ARK_ONO_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getScaled_ARK_ONO_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getScaled_ARK_ONO_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getScaled_ARK_ONO_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getScaled_ARK_ONO_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getScaled_ARK_ONO_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getScaled_ARK_ONO_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getScaled_ARK_ONO_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getScaled_ARK_ONO_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getScaled_ARK_ONO_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getScaled_ARK_ONO_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getScaled_ARK_ONO_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getScaled_ARK_ONO_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getScaled_ARK_ONO_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getScaled_ARK_ONO_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getScaled_ARK_ONO_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getScaled_ARK_ONO_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getScaled_ARK_ONO_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getScaled_ARK_ONO_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getScaled_ARK_ONO_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTH_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getScaled_ARK_OTH_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getScaled_ARK_OTH_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getScaled_ARK_OTH_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getScaled_ARK_OTH_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getScaled_ARK_OTH_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure getScaled_ARK_OTH_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure getScaled_ARK_OTH_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure getScaled_ARK_OTH_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure getScaled_ARK_OTH_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure getScaled_ARK_OTH_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ARK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getScaled_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getScaled_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ACK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getScaled_ACK_ONO_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getScaled_ACK_ONO_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getScaled_ACK_ONO_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getScaled_ACK_ONO_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getScaled_ACK_ONO_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ONO_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getScaled_ACK_ONO_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getScaled_ACK_ONO_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getScaled_ACK_ONO_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getScaled_ACK_ONO_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getScaled_ACK_ONO_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ONO_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define OTH_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure getScaled_ACK_OTH_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure getScaled_ACK_OTH_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure getScaled_ACK_OTH_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure getScaled_ACK_OTH_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure getScaled_ACK_OTH_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OTH_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ACK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getScaled_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setScaled_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ARK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setScaled_ARK_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setScaled_ARK_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setScaled_ARK_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setScaled_ARK_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setScaled_ARK_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setScaled_ARK_D1_RK5
        use pm_kind, only: RKG => RK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setScaled_ARK_D1_RK4
        use pm_kind, only: RKG => RK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setScaled_ARK_D1_RK3
        use pm_kind, only: RKG => RK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setScaled_ARK_D1_RK2
        use pm_kind, only: RKG => RK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setScaled_ARK_D1_RK1
        use pm_kind, only: RKG => RK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setScaled_ARK_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setScaled_ARK_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setScaled_ARK_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setScaled_ARK_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setScaled_ARK_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define RK_ENABLED 1

#if RK5_ENABLED
    module procedure setScaled_ARK_D2_RK5
        use pm_kind, only: RKG => RK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK4_ENABLED
    module procedure setScaled_ARK_D2_RK4
        use pm_kind, only: RKG => RK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK3_ENABLED
    module procedure setScaled_ARK_D2_RK3
        use pm_kind, only: RKG => RK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK2_ENABLED
    module procedure setScaled_ARK_D2_RK2
        use pm_kind, only: RKG => RK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if RK1_ENABLED
    module procedure setScaled_ARK_D2_RK1
        use pm_kind, only: RKG => RK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef RK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ARK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setScaled_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setScaled_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ACK_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D1_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setScaled_ACK_D1_CK5
        use pm_kind, only: CKG => CK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setScaled_ACK_D1_CK4
        use pm_kind, only: CKG => CK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setScaled_ACK_D1_CK3
        use pm_kind, only: CKG => CK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setScaled_ACK_D1_CK2
        use pm_kind, only: CKG => CK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setScaled_ACK_D1_CK1
        use pm_kind, only: CKG => CK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D1_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define D2_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CK_ENABLED 1

#if CK5_ENABLED
    module procedure setScaled_ACK_D2_CK5
        use pm_kind, only: CKG => CK5
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK4_ENABLED
    module procedure setScaled_ACK_D2_CK4
        use pm_kind, only: CKG => CK4
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK3_ENABLED
    module procedure setScaled_ACK_D2_CK3
        use pm_kind, only: CKG => CK3
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK2_ENABLED
    module procedure setScaled_ACK_D2_CK2
        use pm_kind, only: CKG => CK2
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#if CK1_ENABLED
    module procedure setScaled_ACK_D2_CK1
        use pm_kind, only: CKG => CK1
#include "pm_sampleScale@routines.inc.F90"
    end procedure
#endif

#undef CK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef D2_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef ACK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setScaled_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines